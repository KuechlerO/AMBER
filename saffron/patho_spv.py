"""Gutierrez Guarnizo et al. (2023) pathogenic signal-peptide variant catalogue.

Source: NAR Genomics and Bioinformatics — PMC10583284
CSV: patho_SPVs_in_hs.csv (semicolon-delimited VEP + paper annotations).
"""

from __future__ import annotations

import csv
import threading
from functools import lru_cache
from pathlib import Path
from typing import Any

from django.conf import settings


_INDEX_LOCK = threading.Lock()
_INDEX: dict[tuple[str, int, str, str], dict[str, Any]] | None = None


def patho_spv_csv_path() -> Path:
    raw = getattr(settings, 'PATHO_SPV_CSV', None)
    if raw:
        return Path(raw)
    base = Path(getattr(settings, 'BASE_DIR', Path('.')))
    return base / 'files-archive-dir' / 'patho_spv_in_hs' / 'patho_SPVs_in_hs.csv'


def _norm_accession(swissprot: str) -> str:
    return (swissprot or '').split('.')[0].strip().upper()


def _norm_clin_sig(raw: str | None) -> str | None:
    text = (raw or '').strip()
    if not text or text == '-':
        return None
    return text


def _is_pathogenic_flag(raw: str | None) -> bool:
    return (raw or '').strip().lower() == 'yes'


def _load_index(path: Path) -> dict[tuple[str, int, str, str], dict[str, Any]]:
    index: dict[tuple[str, int, str, str], dict[str, Any]] = {}
    if not path.is_file():
        return index

    with path.open(newline='', encoding='utf-8', errors='replace') as fh:
        reader = csv.DictReader(fh, delimiter=';')
        for row in reader:
            acc = _norm_accession(row.get('SWISSPROT') or '')
            wt = (row.get('Reference') or '').strip().upper()
            mut = (row.get('Alternative') or '').strip().upper()
            try:
                pos = int(str(row.get('Protein_position') or '').strip())
            except ValueError:
                continue
            if not acc or not wt or not mut or len(wt) != 1 or len(mut) != 1:
                continue

            clin = _norm_clin_sig(row.get('CLIN_SIG'))
            patho = _is_pathogenic_flag(row.get('Potentially_pathogenic_variant'))
            key = (acc, pos, wt, mut)
            existing = index.get(key)
            if existing is None:
                index[key] = {
                    'clin_sig': clin,
                    'paper_pathogenic': patho,
                    'gene_symbol': (row.get('SYMBOL') or '').strip() or None,
                }
                continue

            # Merge duplicates: pathogenic if any Yes; union CLIN_SIG labels
            existing['paper_pathogenic'] = bool(existing['paper_pathogenic'] or patho)
            if clin:
                prev = existing.get('clin_sig')
                if not prev:
                    existing['clin_sig'] = clin
                elif clin not in prev.split(','):
                    # Keep readable; avoid huge concatenations
                    parts = [p.strip() for p in prev.split(',') if p.strip()]
                    if clin not in parts:
                        parts.append(clin)
                        existing['clin_sig'] = ','.join(parts[:6])
    return index


def get_patho_spv_index() -> dict[tuple[str, int, str, str], dict[str, Any]]:
    global _INDEX
    if _INDEX is not None:
        return _INDEX
    with _INDEX_LOCK:
        if _INDEX is None:
            _INDEX = _load_index(patho_spv_csv_path())
        return _INDEX


def clear_patho_spv_index() -> None:
    """Test helper."""
    global _INDEX
    with _INDEX_LOCK:
        _INDEX = None


def lookup_patho_spv(
    accession: str | None,
    position: int | None,
    wt_aa: str | None,
    mut_aa: str | None,
) -> dict[str, Any] | None:
    if not accession or position is None or not wt_aa or not mut_aa:
        return None
    key = (
        _norm_accession(accession),
        int(position),
        str(wt_aa).strip().upper(),
        str(mut_aa).strip().upper(),
    )
    hit = get_patho_spv_index().get(key)
    if not hit:
        return None
    return dict(hit)


def annotate_guide_rows(
    rows: list[dict[str, Any]],
    *,
    accession: str | None,
) -> list[dict[str, Any]]:
    """Attach clin_sig / paper_pathogenic fields from the catalogue."""
    if not accession or not rows:
        for r in rows:
            r.setdefault('clin_sig', None)
            r.setdefault('paper_pathogenic', False)
            r.setdefault('in_patho_catalogue', False)
        return rows

    # Ensure index is warm once per request
    get_patho_spv_index()
    for r in rows:
        hit = lookup_patho_spv(
            accession,
            r.get('position'),
            r.get('wt_aa'),
            r.get('mut_aa'),
        )
        if hit:
            r['clin_sig'] = hit.get('clin_sig')
            r['paper_pathogenic'] = bool(hit.get('paper_pathogenic'))
            r['in_patho_catalogue'] = True
        else:
            r['clin_sig'] = None
            r['paper_pathogenic'] = False
            r['in_patho_catalogue'] = False
    return rows


@lru_cache(maxsize=1)
def catalogue_citation() -> str:
    return (
        'Gutierrez Guarnizo et al., NAR Genom Bioinform (2023) — '
        'Pathogenic signal peptide variants in the human genome (PMC10583284)'
    )
