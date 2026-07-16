"""Map UniProt accessions to NGG screen library genes (Supplementary Table 1)."""

from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path

import pandas as pd
from django.conf import settings

SCREEN_MARKERS = ['TNFa', 'IFNG', 'PD1', 'CD25']

_UNIPROT_RE = re.compile(r'^([OPQ][0-9][A-Z0-9]{3}[0-9])(?:-\d+)?$')


def _screen_data_root() -> Path:
    return Path(getattr(settings, 'SCREEN_DATA_DIR', settings.BASE_DIR / 'files-archive-dir'))


def _table1_path() -> Path:
    return _screen_data_root() / 'Supplementary Table 1.txt'


def normalize_uniprot_accession(uniprot_id: str) -> str:
    """Normalize user/API IDs for Table 1 matching (case, prefix, isoform suffix)."""
    uid = (uniprot_id or '').strip().upper()
    if not uid:
        return ''
    if ':' in uid:
        uid = uid.rsplit(':', 1)[-1]
    match = _UNIPROT_RE.match(uid)
    if match:
        return match.group(1)
    return uid


def resolve_lookup_uniprot_id(uniprot_id: str) -> str:
    """
    ID used for Supplementary Table 1 lookup.

    Resolves Ensembl transcript IDs to UniProt; normalizes isoforms (P60709-1 → P60709).
    """
    raw = (uniprot_id or '').strip()
    if not raw:
        return ''
    if raw.upper().startswith('ENST'):
        from designer.pipeline import ensembl_to_uniprot

        resolved = ensembl_to_uniprot(raw)
        if resolved:
            return normalize_uniprot_accession(resolved)
        return ''
    return normalize_uniprot_accession(raw)


def is_screen_library_data_readable() -> bool:
    """True if Supplementary Table 1 can be read (permissions / path)."""
    path = _table1_path()
    if not path.is_file():
        return False
    try:
        with path.open('rb') as fh:
            fh.read(1)
        return True
    except OSError:
        return False


@lru_cache(maxsize=1)
def _load_table1() -> pd.DataFrame:
    path = _table1_path()
    if not path.is_file():
        return pd.DataFrame(columns=['gene', 'uniprot_accession', 'uniprot_len', 'len'])
    try:
        return pd.read_csv(path, sep='\t', index_col=0)
    except OSError:
        return pd.DataFrame(columns=['gene', 'uniprot_accession', 'uniprot_len', 'len'])


def lookup_gene_by_uniprot(uniprot_id: str) -> dict | None:
    """
    Return gene metadata if uniprot_id is in Supplementary Table 1, else None.

    Keys: gene, uniprot_accession, uniprot_len (int), protein_len (int).
    """
    uid = resolve_lookup_uniprot_id(uniprot_id) or normalize_uniprot_accession(uniprot_id)
    if not uid:
        return None

    df = _load_table1()
    if df.empty or 'uniprot_accession' not in df.columns:
        return None

    accessions = (
        df['uniprot_accession']
        .astype(str)
        .str.strip()
        .apply(normalize_uniprot_accession)
    )
    matches = df[accessions == uid]
    if matches.empty:
        return None

    row = matches.iloc[0]
    gene = str(row['gene']) if 'gene' in row.index else str(matches.index[0])
    protein_len = int(row.get('uniprot_len', row.get('len', 0)))
    return {
        'gene': gene,
        'uniprot_accession': str(row['uniprot_accession']),
        'uniprot_len': protein_len,
        'protein_len': protein_len,
    }


def screen_plot_eligibility(uniprot_id: str) -> dict:
    """
    Whether the screen plot can be offered for this input.

    Returns keys: available (bool), gene, lookup_uniprot, reason ('not_in_library' | 'data_unavailable' | None).
    """
    lookup_uniprot = resolve_lookup_uniprot_id(uniprot_id)
    if not is_screen_library_data_readable():
        return {
            'available': False,
            'gene': None,
            'lookup_uniprot': lookup_uniprot,
            'reason': 'data_unavailable',
        }
    meta = lookup_gene_by_uniprot(uniprot_id)
    if meta is None:
        return {
            'available': False,
            'gene': None,
            'lookup_uniprot': lookup_uniprot,
            'reason': 'not_in_library',
        }
    return {
        'available': True,
        'gene': meta['gene'],
        'lookup_uniprot': lookup_uniprot or meta['uniprot_accession'],
        'reason': None,
        'meta': meta,
    }


def clear_gene_lookup_cache():
    _load_table1.cache_clear()
