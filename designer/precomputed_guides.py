"""Load guide rows from precomputed Parquet shards (proteome batch output)."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import pandas as pd
from django.conf import settings

from .pipeline import (
    annotate_sgrna_binding,
    finalize_guide_scores,
    lookup_esm1b,
)

PRECOMPUTE_WINDOW_MIN = 1
PRECOMPUTE_WINDOW_MAX = 10
_SUPPORTED_PAMS = frozenset({'NGG', 'NG', 'NAG'})

_SHARD_RE = re.compile(
    r'^(?P<acc>[A-Z0-9]+)_sv(?P<sv>\d+)_tv(?P<tid>ENST\d+(?:\.\d+)?)\.parquet$'
)


def precomputed_guides_dir() -> Path:
    raw = getattr(settings, 'PRECOMPUTED_GUIDES_DIR', None)
    if raw:
        return Path(raw)
    return Path(settings.BASE_DIR) / 'data' / 'precomputed_guides'


def precomputed_lookup_enabled() -> bool:
    if getattr(settings, 'PRECOMPUTED_GUIDES_DISABLE', False):
        return False
    guides = precomputed_guides_dir() / 'guides'
    return guides.is_dir()


def _window_within_precomputed(window_min: int, window_max: int) -> bool:
    return (
        PRECOMPUTE_WINDOW_MIN <= window_min <= window_max <= PRECOMPUTE_WINDOW_MAX
    )


def _list_shards_for_accession(guides_dir: Path, accession: str) -> list[Path]:
    acc = accession.strip().upper()
    out = []
    for p in guides_dir.glob(f'{acc}_sv*_tv*.parquet'):
        if _SHARD_RE.match(p.name):
            out.append(p)
    return sorted(out, key=lambda p: (_shard_sequence_version(p.name), p.stat().st_mtime), reverse=True)


def _shard_sequence_version(name: str) -> int:
    m = _SHARD_RE.match(name)
    return int(m.group('sv')) if m else -1


def _pick_shard(guides_dir: Path, accession: str) -> Path | None:
    shards = _list_shards_for_accession(guides_dir, accession)
    return shards[0] if shards else None


def _parse_target_positions(raw) -> list[int]:
    if raw is None or (isinstance(raw, float) and pd.isna(raw)):
        return []
    if isinstance(raw, list):
        return [int(x) for x in raw]
    if isinstance(raw, str):
        text = raw.strip()
        if not text:
            return []
        try:
            parsed = json.loads(text)
            if isinstance(parsed, list):
                return [int(x) for x in parsed]
        except json.JSONDecodeError:
            pass
    return []


def _guide_intersects_window(positions: list[int], window_min: int, window_max: int) -> bool:
    if not positions:
        return False
    return any(window_min <= p <= window_max for p in positions)


def _parse_outcome_details(raw) -> list[dict[str, Any]]:
    if raw is None or (isinstance(raw, float) and pd.isna(raw)):
        return []
    if isinstance(raw, list):
        return [dict(x) for x in raw]
    if isinstance(raw, str):
        try:
            parsed = json.loads(raw)
            if isinstance(parsed, list):
                return [dict(x) for x in parsed]
        except json.JSONDecodeError:
            return []
    return []


def _prefers_center(positions) -> bool:
    if isinstance(positions, list):
        return any(p in (5, 6) for p in positions)
    return positions in (5, 6)


def _editors_for_request(editor: str) -> set[str]:
    if editor == 'BOTH':
        return {'ABE', 'CBE'}
    if editor in ('ABE', 'CBE'):
        return {editor}
    return set()


def _attach_esm1b(outcomes: list[dict], *, uniprot_id: str, position: int, wt_aa: str) -> None:
    for detail in outcomes:
        aa = detail.get('aa')
        if not aa:
            continue
        detail['esm1b'] = lookup_esm1b(uniprot_id, position, aa, wt_aa)


def _dataframe_to_pipeline_rows(
    df: pd.DataFrame,
    *,
    editor: str,
    pam_type: str,
    window_min: int,
    window_max: int,
    alpha_threshold: float,
    top_sgrnas: int,
    uniprot_id: str,
) -> list[dict[str, Any]]:
    editors = _editors_for_request(editor)
    pam = (pam_type or 'NGG').upper()
    if pam not in _SUPPORTED_PAMS:
        return []

    filtered = df[
        df['editor_used'].isin(editors) & (df['pam_type'].astype(str).str.upper() == pam)
    ].copy()
    if filtered.empty:
        return []

    rows: list[dict[str, Any]] = []
    for _, rec in filtered.iterrows():
        try:
            alpha = float(rec['alpha_score'])
        except (TypeError, ValueError):
            continue
        if alpha < alpha_threshold:
            continue
        target_positions = _parse_target_positions(rec.get('target_positions'))
        if not _guide_intersects_window(target_positions, window_min, window_max):
            continue
        position = int(rec['position'])
        wt_aa = rec.get('wt_aa')
        outcomes = _parse_outcome_details(rec.get('outcome_details'))
        # Drop outcomes below threshold for display consistency
        outcomes = [
            d for d in outcomes
            if d.get('score') is not None and float(d['score']) >= alpha_threshold
        ]
        if not outcomes:
            continue
        _attach_esm1b(outcomes, uniprot_id=uniprot_id, position=position, wt_aa=wt_aa)

        rows.append({
            'position': position,
            'wt_codon': rec.get('wt_codon'),
            'wt_aa': wt_aa,
            'mut_aa': rec.get('mut_aa'),
            'alpha_score': alpha,
            'avg_alpha_score': rec.get('avg_alpha_score'),
            'sgrna_seq': rec.get('sgrna_seq'),
            'protospacer': rec.get('protospacer'),
            'pam': rec.get('pam'),
            'strand': rec.get('strand'),
            'target_position': target_positions,
            'editor_used': rec.get('editor_used'),
            'num_as_window': rec.get('num_editable_window'),
            '_outcome_details': outcomes,
        })

    # top_sgrnas per (position, editor)
    try:
        n_top = int(top_sgrnas)
    except (TypeError, ValueError):
        n_top = 5

    grouped: dict[tuple[int, str], list[dict]] = {}
    for row in rows:
        key = (int(row['position']), str(row['editor_used']))
        grouped.setdefault(key, []).append(row)

    trimmed: list[dict[str, Any]] = []
    for key in sorted(grouped.keys()):
        sg_entries = grouped[key]
        sg_entries.sort(
            key=lambda g: (
                0 if _prefers_center(g.get('target_position')) else 1,
                g.get('num_as_window') if g.get('num_as_window') is not None else 999,
            )
        )
        trimmed.extend(sg_entries[:n_top])
    trimmed.sort(key=lambda r: r.get('position') if r.get('position') is not None else 999999)
    return trimmed


def _protein_length_from_cache(out_dir: Path, shard_path: Path) -> int:
    cache_name = shard_path.name.replace('.parquet', '.json')
    cache_path = out_dir / 'cds_cache' / cache_name
    if not cache_path.is_file():
        return 0
    try:
        with open(cache_path, encoding='utf-8') as f:
            data = json.load(f)
        cds = str(data.get('cds') or '')
        return len(cds) // 3 if cds else 0
    except (OSError, json.JSONDecodeError, TypeError):
        return 0


def _metadata_from_shard(df: pd.DataFrame, shard_path: Path, out_dir: Path) -> dict[str, Any]:
    if df.empty:
        return {}
    row = df.iloc[0]
    protein_length = _protein_length_from_cache(out_dir, shard_path)
    return {
        'uniprot_accession': str(row.get('uniprot_accession') or '').upper(),
        'gene_symbol': str(row.get('gene_symbol') or ''),
        'transcript_id': str(row.get('transcript_id') or ''),
        'protein_length': protein_length,
        'precomputed_shard': shard_path.name,
    }


def try_load_precomputed_guides(
    uniprot_id: str,
    editor: str,
    alpha_threshold,
    top_sgrnas,
    window_min,
    window_max,
    pam_type='NGG',
    *,
    include_cadd=False,
) -> dict[str, Any] | None:
    """Return run_pipeline-shaped result from Parquet, or None to fall back to live design."""
    if not precomputed_lookup_enabled():
        return None

    try:
        wmin = int(window_min)
        wmax = int(window_max)
    except (TypeError, ValueError):
        return None
    if not _window_within_precomputed(wmin, wmax):
        return None

    acc = (uniprot_id or '').strip().upper()
    if not acc or acc.startswith('ENST'):
        return None

    pam = (pam_type or 'NGG').upper()
    if pam not in _SUPPORTED_PAMS:
        return None

    try:
        threshold = float(alpha_threshold)
    except (TypeError, ValueError):
        threshold = 0.0

    out_dir = precomputed_guides_dir()
    guides_dir = out_dir / 'guides'
    shard = _pick_shard(guides_dir, acc)
    if shard is None:
        return None

    try:
        df = pd.read_parquet(shard)
    except Exception:
        return None
    if df.empty:
        return None

    meta = _metadata_from_shard(df, shard, out_dir)
    rows = _dataframe_to_pipeline_rows(
        df,
        editor=editor,
        pam_type=pam,
        window_min=wmin,
        window_max=wmax,
        alpha_threshold=threshold,
        top_sgrnas=top_sgrnas,
        uniprot_id=acc,
    )
    if not rows:
        return None

    finalize_guide_scores(rows, {}, keep_outcome_details=True)
    cadd_status = 'skipped'
    cadd_warning = None
    if include_cadd:
        cadd_status = 'pending'
        cadd_warning = (
            'CADD scores are loading in the background and will appear when ready. '
            'Precomputed guides may lack genomic SNV metadata for some outcomes.'
        )
    else:
        cadd_status = 'pending'
        cadd_warning = (
            'CADD scores are loading in the background and will appear when ready.'
        )

    return {
        'guide_rows': annotate_sgrna_binding(rows),
        'no_guide_positions': [],
        'protein_length': meta.get('protein_length') or 0,
        'uniprot_accession': meta.get('uniprot_accession') or acc,
        'gene_symbol': meta.get('gene_symbol') or '',
        'transcript_id': meta.get('transcript_id') or '',
        'cadd_warning': cadd_warning,
        'cadd_status': cadd_status,
        'precomputed_source': True,
        'precomputed_shard': meta.get('precomputed_shard'),
    }
