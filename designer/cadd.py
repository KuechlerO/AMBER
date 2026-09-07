"""CADD GRCh38 PHRED lookup via the BIH REST API, with on-disk cache."""

from __future__ import annotations

import json
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import requests
from django.conf import settings


Snv = tuple[str, int, str, str]


def _model_name() -> str:
    return str(getattr(settings, 'CADD_MODEL', 'GRCh38-v1.7'))


def _api_base() -> str:
    return str(getattr(settings, 'CADD_API_BASE', 'https://cadd.bihealth.org/api/v1.0')).rstrip('/')


def cadd_cache_dir() -> Path:
    raw = getattr(settings, 'CADD_CACHE_DIR', None)
    if raw:
        path = Path(raw)
    else:
        cache_root = Path(getattr(settings, 'CACHE_DIR', Path(tempfile.gettempdir()) / 'amber_django_cache'))
        path = cache_root / 'cadd'
    path.mkdir(parents=True, exist_ok=True)
    return path


def _position_cache_path(chrom: str, pos: int) -> Path:
    safe_chrom = str(chrom).replace('/', '_')
    return cadd_cache_dir() / _model_name() / f'{safe_chrom}_{int(pos)}.json'


def _read_position_cache(chrom: str, pos: int) -> list[dict[str, Any]] | None:
    path = _position_cache_path(chrom, pos)
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding='utf-8'))
    except (OSError, json.JSONDecodeError):
        return None
    if isinstance(payload, list):
        return payload
    return None


def _write_position_cache(chrom: str, pos: int, records: list[dict[str, Any]]) -> None:
    path = _position_cache_path(chrom, pos)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.json.tmp')
    tmp.write_text(json.dumps(records), encoding='utf-8')
    tmp.replace(path)


def parse_cadd_payload(payload: Any) -> list[dict[str, Any]]:
    """Normalize either list-of-dicts or header+rows CADD API responses."""
    if not payload:
        return []
    if isinstance(payload, list) and payload and isinstance(payload[0], dict):
        records = []
        for row in payload:
            rec = _record_from_mapping(row)
            if rec:
                records.append(rec)
        return records
    if isinstance(payload, list) and payload and isinstance(payload[0], list):
        header = [str(h) for h in payload[0]]
        records = []
        for raw in payload[1:]:
            if not isinstance(raw, list):
                continue
            row = {header[i]: raw[i] for i in range(min(len(header), len(raw)))}
            rec = _record_from_mapping(row)
            if rec:
                records.append(rec)
        return records
    return []


def _record_from_mapping(row: dict[str, Any]) -> dict[str, Any] | None:
    alt = str(row.get('Alt') or row.get('alt') or '').upper()
    ref = str(row.get('Ref') or row.get('ref') or '').upper()
    chrom = str(row.get('Chrom') or row.get('chrom') or '')
    pos_raw = row.get('Pos') if row.get('Pos') is not None else row.get('pos')
    phred_raw = row.get('PHRED') if row.get('PHRED') is not None else row.get('phred')
    raw_score = row.get('RawScore') if row.get('RawScore') is not None else row.get('rawscore')
    if not alt or pos_raw in (None, ''):
        return None
    try:
        pos = int(pos_raw)
    except (TypeError, ValueError):
        return None
    phred = _to_float(phred_raw)
    return {
        'chrom': chrom,
        'pos': pos,
        'ref': ref,
        'alt': alt,
        'phred': phred,
        'raw': _to_float(raw_score),
    }


def _to_float(value: Any) -> float | None:
    if value in (None, ''):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _fetch_position(chrom: str, pos: int, http_get=None, sleeper=None) -> list[dict[str, Any]]:
    cached = _read_position_cache(chrom, pos)
    if cached is not None:
        return cached
    getter = http_get or requests.get
    sleep = sleeper or time.sleep
    url = f'{_api_base()}/{_model_name()}/{chrom}:{int(pos)}'
    timeout = int(getattr(settings, 'CADD_TIMEOUT_SEC', 20))
    last_exc: Exception | None = None
    for attempt in range(3):
        try:
            response = getter(url, timeout=timeout)
            if response.status_code == 429:
                last_exc = RuntimeError('CADD API rate-limited (429)')
                sleep(1.5 * (attempt + 1))
                continue
            response.raise_for_status()
            records = parse_cadd_payload(response.json())
            _write_position_cache(chrom, pos, records)
            return records
        except requests.ConnectionError as exc:
            # Host unreachable (e.g. cadd.bihealth.org refused) — do not burn retries.
            last_exc = exc
            break
        except (requests.Timeout, requests.RequestException, ValueError) as exc:
            last_exc = exc
            if attempt + 1 < 3:
                sleep(1.5 * (attempt + 1))
    if last_exc:
        raise last_exc
    return []


def fetch_cadd_scores(
    snvs: list[Snv],
    *,
    http_get=None,
    sleeper=None,
) -> tuple[dict[Snv, dict[str, Any]], str | None]:
    """Batch-fetch CADD scores for unique genomic SNVs.

    Returns (snv -> {phred, raw}, warning). Failures leave scores missing
    and set a short warning; they do not raise.
    """
    score_map: dict[Snv, dict[str, Any]] = {}
    if not snvs:
        return score_map, None

    unique_positions = {(chrom, int(pos)) for chrom, pos, _ref, _alt in snvs}
    max_workers = max(1, int(getattr(settings, 'CADD_MAX_WORKERS', 4)))
    workers = min(max_workers, len(unique_positions))
    failures = 0
    by_position: dict[tuple[str, int], list[dict[str, Any]]] = {}

    def _one(item: tuple[str, int]) -> tuple[tuple[str, int], list[dict[str, Any]] | None]:
        chrom, pos = item
        try:
            return item, _fetch_position(chrom, pos, http_get=http_get, sleeper=sleeper)
        except Exception:
            return item, None

    if workers == 1:
        results = [_one(item) for item in unique_positions]
    else:
        results = []
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(_one, item) for item in unique_positions]
            for fut in as_completed(futures):
                results.append(fut.result())

    for item, records in results:
        if records is None:
            failures += 1
            continue
        by_position[item] = records

    for chrom, pos, ref, alt in snvs:
        records = by_position.get((chrom, int(pos))) or []
        match = next(
            (
                rec for rec in records
                if rec.get('alt') == alt.upper()
                and (not rec.get('ref') or rec.get('ref') == ref.upper())
            ),
            None,
        )
        if match is None:
            match = next((rec for rec in records if rec.get('alt') == alt.upper()), None)
        if match:
            score_map[(chrom, int(pos), ref.upper(), alt.upper())] = {
                'phred': match.get('phred'),
                'raw': match.get('raw'),
            }

    warning = None
    if failures:
        warning = (
            f'CADD scores unavailable for {failures} genomic position'
            f'{"s" if failures != 1 else ""} (API timeout or error). '
            'Guide ranking still uses AlphaMissense.'
        )
    return score_map, warning


def max_cadd_phred(snvs: list[Snv], score_map: dict[Snv, dict[str, Any]]) -> float | None:
    values = []
    for snv in snvs:
        rec = score_map.get(_norm_snv(snv))
        if rec and rec.get('phred') is not None:
            values.append(float(rec['phred']))
    return max(values) if values else None


def cadd_raw_for_best(snvs: list[Snv], score_map: dict[Snv, dict[str, Any]]) -> float | None:
    best_phred = None
    best_raw = None
    for snv in snvs:
        rec = score_map.get(_norm_snv(snv))
        if not rec or rec.get('phred') is None:
            continue
        phred = float(rec['phred'])
        if best_phred is None or phred > best_phred:
            best_phred = phred
            best_raw = rec.get('raw')
    return float(best_raw) if best_raw is not None else None


def _norm_snv(snv: Snv) -> Snv:
    chrom, pos, ref, alt = snv
    return (str(chrom), int(pos), str(ref).upper(), str(alt).upper())


def attach_cadd_to_outcomes(
    details: list[dict[str, Any]],
    score_map: dict[Snv, dict[str, Any]],
) -> list[dict[str, Any]]:
    """Set each outcome's `cadd` PHRED from its SNVs (max if multi-base)."""
    for detail in details:
        detail['cadd'] = max_cadd_phred(detail.get('snvs') or [], score_map)
    return details
