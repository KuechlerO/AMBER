"""Load Brandes/Ntranos ESM1b LLR matrices (one CSV per UniProt isoform)."""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import pandas as pd
from django.conf import settings


def esm1b_data_dir() -> Path:
    raw = getattr(settings, 'ESM1B_DATA_DIR', None)
    if raw:
        return Path(raw)
    base = Path(getattr(settings, 'DATA_DIR', Path('.')))
    return base / 'ESM1B' / 'content' / 'ALL_hum_isoforms_ESM1b_LLR'


def canonical_uniprot_id(uniprot_id: str) -> str:
    uid = (uniprot_id or '').strip()
    if '-' in uid:
        base, suffix = uid.rsplit('-', 1)
        if suffix.isdigit():
            return base
    return uid


def resolve_esm1b_path(uniprot_id: str, data_dir: Path | None = None) -> Path | None:
    """Exact isoform file first, then canonical accession without -N."""
    root = data_dir if data_dir is not None else esm1b_data_dir()
    uid = (uniprot_id or '').strip()
    if not uid:
        return None
    candidates = [uid]
    canon = canonical_uniprot_id(uid)
    if canon not in candidates:
        candidates.append(canon)
    for accession in candidates:
        path = root / f'{accession}_LLR.csv'
        if path.is_file():
            return path
    return None


def _parse_position_columns(df: pd.DataFrame) -> tuple[dict[int, str], dict[int, str]]:
    pos_to_wt: dict[int, str] = {}
    pos_to_col: dict[int, str] = {}
    for col in df.columns:
        parts = str(col).split()
        if len(parts) != 2 or not parts[1].isdigit():
            continue
        position = int(parts[1])
        pos_to_wt[position] = parts[0]
        pos_to_col[position] = str(col)
    return pos_to_wt, pos_to_col


@lru_cache(maxsize=32)
def _load_matrix(path_str: str) -> tuple[dict[int, str], dict[int, str], pd.DataFrame] | None:
    path = Path(path_str)
    if not path.is_file():
        return None
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.astype(str).str.strip()
    pos_to_wt, pos_to_col = _parse_position_columns(df)
    return pos_to_wt, pos_to_col, df


def clear_esm1b_cache() -> None:
    _load_matrix.cache_clear()


def lookup_esm1b(
    uniprot_id: str,
    position: int,
    mut_aa: str,
    wt_aa: str | None = None,
    *,
    data_dir: Path | None = None,
) -> float | None:
    """Return ESM1b LLR for a missense, or None if missing / WT mismatch."""
    path = resolve_esm1b_path(uniprot_id, data_dir=data_dir)
    if path is None:
        return None
    loaded = _load_matrix(str(path))
    if loaded is None:
        return None
    pos_to_wt, pos_to_col, df = loaded
    try:
        pos = int(position)
    except (TypeError, ValueError):
        return None
    col = pos_to_col.get(pos)
    if col is None:
        return None
    if wt_aa and pos_to_wt.get(pos) != wt_aa:
        return None
    mut = (mut_aa or '').strip()
    if mut not in df.index:
        return None
    value = df.loc[mut, col]
    if pd.isna(value):
        return None
    return float(value)
