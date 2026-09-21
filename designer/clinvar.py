"""ClinVar protein-variant lookup via NCBI E-utilities, with on-disk cache."""

from __future__ import annotations

import json
import logging
import re
import tempfile
import time
from pathlib import Path
from typing import Any

import requests
from django.conf import settings

logger = logging.getLogger(__name__)

EUTILS_BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
_ONE_LETTER_RE = re.compile(r'^([A-Z*])(\d+)', re.IGNORECASE)
_THREE_LETTER_RE = re.compile(
    r'\(p\.([A-Za-z]{3}|\*)(\d+)',
)
_CACHE_VERSION = 1

# NCBI asks clients to identify themselves; override via settings/env.
_DEFAULT_TOOL = 'amber-crispr-tool'
_DEFAULT_EMAIL = 'amber-clinvar@localhost'


def clinvar_cache_dir() -> Path:
    raw = getattr(settings, 'CLINVAR_CACHE_DIR', None)
    if raw:
        path = Path(raw)
    else:
        cache_root = Path(
            getattr(settings, 'CACHE_DIR', Path(tempfile.gettempdir()) / 'amber_django_cache')
        )
        path = cache_root / 'clinvar'
    path.mkdir(parents=True, exist_ok=True)
    return path


def _cache_ttl_sec() -> int:
    return int(getattr(settings, 'CLINVAR_CACHE_TTL_SEC', 7 * 24 * 3600))


def _timeout_sec() -> float:
    return float(getattr(settings, 'CLINVAR_TIMEOUT_SEC', 20))


def _eutils_params() -> dict[str, str]:
    params = {
        'tool': str(getattr(settings, 'CLINVAR_EUTILS_TOOL', _DEFAULT_TOOL) or _DEFAULT_TOOL),
        'email': str(getattr(settings, 'CLINVAR_EUTILS_EMAIL', _DEFAULT_EMAIL) or _DEFAULT_EMAIL),
    }
    api_key = getattr(settings, 'NCBI_API_KEY', None) or ''
    api_key = str(api_key).strip().strip('"').strip("'")
    if api_key:
        params['api_key'] = api_key
    return params


def _cache_path(gene_symbol: str) -> Path:
    safe = re.sub(r'[^A-Za-z0-9_.-]+', '_', gene_symbol.upper())
    return clinvar_cache_dir() / f'{safe}.json'


def _read_cache(gene_symbol: str) -> list[dict[str, Any]] | None:
    path = _cache_path(gene_symbol)
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding='utf-8'))
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(payload, dict) or payload.get('version') != _CACHE_VERSION:
        return None
    fetched_at = payload.get('fetched_at')
    try:
        age = time.time() - float(fetched_at)
    except (TypeError, ValueError):
        return None
    if age > _cache_ttl_sec():
        return None
    variants = payload.get('variants')
    if not isinstance(variants, list):
        return None
    return variants


def _write_cache(gene_symbol: str, variants: list[dict[str, Any]]) -> None:
    path = _cache_path(gene_symbol)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        'version': _CACHE_VERSION,
        'gene_symbol': gene_symbol.upper(),
        'fetched_at': time.time(),
        'variants': variants,
    }
    tmp = path.with_suffix('.json.tmp')
    tmp.write_text(json.dumps(payload), encoding='utf-8')
    tmp.replace(path)


def _http_get(url: str, *, params: dict[str, Any] | None = None, retries: int = 2) -> dict[str, Any]:
    merged = dict(_eutils_params())
    if params:
        merged.update(params)
    timeout = _timeout_sec()
    last_exc: Exception | None = None
    for attempt in range(retries + 1):
        try:
            resp = requests.get(url, params=merged, timeout=timeout)
            resp.raise_for_status()
            return resp.json()
        except (requests.Timeout, requests.ConnectionError, requests.HTTPError, ValueError) as exc:
            last_exc = exc
            if attempt < retries:
                time.sleep(0.4 * (attempt + 1))
                continue
            raise
    raise RuntimeError(str(last_exc or 'ClinVar request failed'))


def parse_residue_from_protein_change(
    protein_change: str | None,
    *,
    protein_length: int | None = None,
) -> int | None:
    """Pick a 1-based residue from ClinVar protein_change (e.g. L260V, R514*)."""
    if not protein_change or not str(protein_change).strip():
        return None
    candidates: list[int] = []
    for part in str(protein_change).split(','):
        token = part.strip()
        if not token:
            continue
        match = _ONE_LETTER_RE.match(token)
        if not match:
            continue
        try:
            residue = int(match.group(2))
        except ValueError:
            continue
        if residue <= 0:
            continue
        if protein_length is not None and residue > int(protein_length):
            continue
        candidates.append(residue)
    if not candidates:
        return None
    return max(candidates)


def parse_residue_from_title(
    title: str | None,
    *,
    protein_length: int | None = None,
) -> int | None:
    """Fallback: parse p.Leu260Val / p.Arg514* from ClinVar title."""
    if not title:
        return None
    match = _THREE_LETTER_RE.search(title)
    if not match:
        return None
    try:
        residue = int(match.group(2))
    except ValueError:
        return None
    if residue <= 0:
        return None
    if protein_length is not None and residue > int(protein_length):
        return None
    return residue


def _normalize_classification(description: str | None) -> str:
    text = (description or '').strip()
    if not text:
        return 'Unknown'
    lower = text.lower()
    if 'conflict' in lower:
        return 'Conflicting classifications'
    if lower == 'pathogenic':
        return 'Pathogenic'
    if lower == 'likely pathogenic':
        return 'Likely pathogenic'
    if 'uncertain' in lower or lower in ('vus', 'variant of uncertain significance'):
        return 'Uncertain significance'
    if lower == 'likely benign':
        return 'Likely benign'
    if lower == 'benign':
        return 'Benign'
    return text


def parse_esummary_record(
    record: dict[str, Any],
    *,
    protein_length: int | None = None,
) -> dict[str, Any] | None:
    """Convert one ClinVar esummary JSON object into a protein-mapped variant."""
    if not isinstance(record, dict):
        return None
    protein_change = record.get('protein_change') or ''
    title = record.get('title') or ''
    residue = parse_residue_from_protein_change(protein_change, protein_length=protein_length)
    if residue is None:
        residue = parse_residue_from_title(title, protein_length=protein_length)
    if residue is None:
        return None

    germline = record.get('germline_classification') or {}
    if not isinstance(germline, dict):
        germline = {}
    classification = _normalize_classification(germline.get('description'))
    review_status = str(germline.get('review_status') or '').strip()

    mol = record.get('molecular_consequence_list') or []
    if isinstance(mol, list):
        molecular = [str(x) for x in mol if x]
    else:
        molecular = []

    uid = str(record.get('uid') or '').strip()
    accession = str(record.get('accession') or '').strip()
    return {
        'uid': uid,
        'accession': accession,
        'title': title,
        'protein_change': str(protein_change).strip(),
        'residue': int(residue),
        'classification': classification,
        'review_status': review_status,
        'molecular_consequence': molecular,
        'url': f'https://www.ncbi.nlm.nih.gov/clinvar/variation/{uid}/' if uid else '',
    }


def _esearch_uids(gene_symbol: str) -> list[str]:
    url = f'{EUTILS_BASE}/esearch.fcgi'
    retmax = 500
    retstart = 0
    ids: list[str] = []
    total: int | None = None
    while True:
        payload = _http_get(
            url,
            params={
                'db': 'clinvar',
                'term': f'{gene_symbol}[gene]',
                'retmode': 'json',
                'retmax': retmax,
                'retstart': retstart,
            },
        )
        result = payload.get('esearchresult') or {}
        batch = [str(x) for x in (result.get('idlist') or []) if x]
        if total is None:
            try:
                total = int(result.get('count') or 0)
            except (TypeError, ValueError):
                total = 0
        ids.extend(batch)
        retstart += len(batch)
        if not batch or (total is not None and retstart >= total):
            break
        # Be polite to NCBI when paging large genes.
        time.sleep(0.34)
    # Preserve order, drop duplicates.
    seen: set[str] = set()
    unique: list[str] = []
    for uid in ids:
        if uid in seen:
            continue
        seen.add(uid)
        unique.append(uid)
    return unique


def _esummary_records(uids: list[str]) -> list[dict[str, Any]]:
    if not uids:
        return []
    url = f'{EUTILS_BASE}/esummary.fcgi'
    records: list[dict[str, Any]] = []
    chunk_size = 200
    for i in range(0, len(uids), chunk_size):
        chunk = uids[i:i + chunk_size]
        payload = _http_get(
            url,
            params={
                'db': 'clinvar',
                'id': ','.join(chunk),
                'retmode': 'json',
            },
        )
        result = payload.get('result') or {}
        for uid in result.get('uids') or chunk:
            rec = result.get(str(uid))
            if isinstance(rec, dict):
                records.append(rec)
        if i + chunk_size < len(uids):
            time.sleep(0.34)
    return records


def fetch_clinvar_protein_variants(
    gene_symbol: str | None,
    *,
    protein_length: int | None = None,
    use_cache: bool = True,
) -> tuple[list[dict[str, Any]], str | None]:
    """
    Return (variants, error).

    Soft-fails: on missing gene or network/parse problems returns ([], message).
    """
    symbol = (gene_symbol or '').strip()
    if not symbol:
        return [], 'ClinVar lookup needs a gene symbol'

    if use_cache:
        cached = _read_cache(symbol)
        if cached is not None:
            # Re-filter by current protein_length in case length changed.
            if protein_length is not None:
                plen = int(protein_length)
                filtered = [v for v in cached if 1 <= int(v.get('residue') or 0) <= plen]
                return filtered, None
            return cached, None

    try:
        uids = _esearch_uids(symbol)
        raw_records = _esummary_records(uids)
        variants: list[dict[str, Any]] = []
        for rec in raw_records:
            parsed = parse_esummary_record(rec, protein_length=protein_length)
            if parsed:
                variants.append(parsed)
        variants.sort(key=lambda v: (int(v['residue']), v.get('classification') or '', v.get('uid') or ''))
        if use_cache:
            try:
                _write_cache(symbol, variants)
            except OSError as exc:
                logger.warning('ClinVar cache write failed for %s: %s', symbol, exc)
        return variants, None
    except Exception as exc:
        logger.warning('ClinVar fetch failed for %s: %s', symbol, exc)
        return [], f'ClinVar unavailable: {exc}'


_CLASSIFICATION_CSS = {
    'Pathogenic': 'clinvar-pathogenic',
    'Likely pathogenic': 'clinvar-likely-pathogenic',
    'Uncertain significance': 'clinvar-vus',
    'Likely benign': 'clinvar-likely-benign',
    'Benign': 'clinvar-benign',
    'Conflicting classifications': 'clinvar-conflict',
}


def classification_css_class(classification: str | None) -> str:
    key = (classification or '').strip()
    if key in _CLASSIFICATION_CSS:
        return _CLASSIFICATION_CSS[key]
    lower = key.lower()
    if 'conflict' in lower:
        return 'clinvar-conflict'
    if lower == 'pathogenic':
        return 'clinvar-pathogenic'
    if lower == 'likely pathogenic':
        return 'clinvar-likely-pathogenic'
    if 'uncertain' in lower:
        return 'clinvar-vus'
    if lower == 'likely benign':
        return 'clinvar-likely-benign'
    if lower == 'benign':
        return 'clinvar-benign'
    return 'clinvar-vus'


def index_clinvar_by_residue(variants: list[dict[str, Any]]) -> dict[int, list[dict[str, Any]]]:
    """Group parsed ClinVar variants by 1-based residue."""
    by_res: dict[int, list[dict[str, Any]]] = {}
    for raw in variants or []:
        try:
            residue = int(raw.get('residue'))
        except (TypeError, ValueError):
            continue
        entry = dict(raw)
        entry['css_class'] = classification_css_class(entry.get('classification'))
        by_res.setdefault(residue, []).append(entry)
    for residue in by_res:
        by_res[residue].sort(
            key=lambda v: (v.get('classification') or '', v.get('protein_change') or '', v.get('uid') or '')
        )
    return by_res


def format_clinvar_export(variants: list[dict[str, Any]] | None) -> str:
    """Semicolon-joined protein_change|classification|accession for CSV/Excel."""
    parts = []
    for v in variants or []:
        change = v.get('protein_change') or '—'
        classification = v.get('classification') or 'Unknown'
        accession = v.get('accession') or v.get('uid') or ''
        parts.append(f'{change}|{classification}|{accession}')
    return '; '.join(parts)


def attach_clinvar_to_rows(
    rows: list[dict[str, Any]] | None,
    by_residue: dict[int, list[dict[str, Any]]] | None,
) -> list[dict[str, Any]]:
    """Annotate each row with clinvar_at_position (list) and clinvar (export string)."""
    index = by_residue or {}
    annotated: list[dict[str, Any]] = []
    for row in rows or []:
        out = dict(row)
        try:
            pos = int(out.get('position'))
        except (TypeError, ValueError):
            pos = None
        variants = list(index.get(pos, [])) if pos is not None else []
        out['clinvar_at_position'] = variants
        out['clinvar'] = format_clinvar_export(variants)
        annotated.append(out)
    return annotated
