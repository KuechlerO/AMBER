"""Resolve human gene symbols to UniProt accessions via the UniProt REST API."""

from __future__ import annotations

import re
from typing import Any

import requests

from .pipeline import _http_get
from .services import looks_like_uniprot_accession, normalize_input_id

_GENE_SYMBOL_RE = re.compile(r'^[A-Z0-9][A-Z0-9_-]{0,31}$', re.IGNORECASE)


class GeneSearchError(ValueError):
    """User-facing gene search failure."""


def looks_like_gene_symbol(raw: str) -> bool:
    """True if input looks like a gene symbol (not UniProt / Ensembl)."""
    symbol = normalize_input_id(raw)
    if not symbol or len(symbol) < 2:
        return False
    if looks_like_uniprot_accession(symbol) or symbol.startswith('ENST'):
        return False
    return bool(_GENE_SYMBOL_RE.match(symbol))


def _protein_name(entry: dict[str, Any]) -> str:
    desc = entry.get('proteinDescription') or {}
    recommended = ((desc.get('recommendedName') or {}).get('fullName') or {}).get('value')
    if recommended:
        return str(recommended)
    for sub in desc.get('submissionNames') or []:
        value = ((sub.get('fullName') or {}).get('value'))
        if value:
            return str(value)
    return ''


def _gene_name(entry: dict[str, Any]) -> str:
    for gene in entry.get('genes') or []:
        name = (gene.get('geneName') or {}).get('value')
        if name:
            return str(name)
    return ''


def _is_reviewed(entry: dict[str, Any]) -> bool:
    entry_type = str(entry.get('entryType') or '').lower()
    # Prefer Swiss-Prot; avoid matching the substring "reviewed" inside "unreviewed".
    return 'swiss-prot' in entry_type or (
        'reviewed' in entry_type and 'unreviewed' not in entry_type
    )


def _parse_entry(entry: dict[str, Any]) -> dict[str, Any] | None:
    accession = str(entry.get('primaryAccession') or '').strip().upper()
    if not accession:
        return None
    length = (entry.get('sequence') or {}).get('length')
    try:
        length_int = int(length) if length is not None else None
    except (TypeError, ValueError):
        length_int = None
    return {
        'accession': accession,
        'gene': _gene_name(entry),
        'protein_name': _protein_name(entry),
        'length': length_int,
        'reviewed': _is_reviewed(entry),
        'am_available': False,
    }


def alphamissense_accessions_present(accessions: list[str]) -> set[str]:
    """Return the subset of accessions that have AlphaMissense rows in legacy_db."""
    ids = [str(a).strip().upper() for a in accessions if a]
    if not ids:
        return set()
    from .models import Alpha_missense

    return set(
        Alpha_missense.objects.using('legacy_db')
        .filter(uniprot_id__in=ids)
        .values_list('uniprot_id', flat=True)
        .distinct()
    )


def annotate_am_availability(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Attach am_available and re-sort: Swiss-Prot + AM first."""
    present = alphamissense_accessions_present([r.get('accession', '') for r in rows])
    for row in rows:
        acc = str(row.get('accession') or '').strip().upper()
        row['am_available'] = acc in present
    rows.sort(
        key=lambda row: (
            not row.get('reviewed', False),
            not row.get('am_available', False),
            row.get('accession') or '',
        )
    )
    return rows


def _search_query(query: str, *, size: int) -> list[dict[str, Any]]:
    url = 'https://rest.uniprot.org/uniprotkb/search'
    params = {
        'query': query,
        'fields': 'accession,gene_names,protein_name,length',
        'format': 'json',
        'size': size,
    }
    r = _http_get(url, params=params, timeout=30)
    r.raise_for_status()
    payload = r.json() or {}
    rows: list[dict[str, Any]] = []
    seen: set[str] = set()
    for entry in payload.get('results') or []:
        parsed = _parse_entry(entry)
        if not parsed or parsed['accession'] in seen:
            continue
        seen.add(parsed['accession'])
        rows.append(parsed)
    return rows


def search_uniprot_by_gene_symbol(
    symbol: str,
    *,
    organism_id: int = 9606,
    limit: int = 20,
) -> list[dict[str, Any]]:
    """
    Search UniProt for human (default) proteins matching a gene symbol.

    Swiss-Prot (reviewed) entries with AlphaMissense coverage are returned first.
    Raises GeneSearchError for invalid input; returns [] when UniProt has no matches.
    """
    raw = (symbol or '').strip()
    if not raw:
        raise GeneSearchError('Please enter a gene symbol.')
    normalized = normalize_input_id(raw)
    if looks_like_uniprot_accession(normalized):
        raise GeneSearchError(
            'That looks like a UniProt accession. Enter it directly, or use a gene symbol.'
        )
    if normalized.startswith('ENST'):
        raise GeneSearchError(
            'That looks like an Ensembl transcript ID. Enter it directly, or use a gene symbol.'
        )
    if not looks_like_gene_symbol(normalized):
        raise GeneSearchError('Invalid gene symbol.')

    size = max(1, min(int(limit), 50))
    exact_query = f'gene_exact:{normalized} AND organism_id:{organism_id}'
    try:
        rows = _search_query(exact_query, size=size)
        if not rows:
            # Broader gene: match when exact returns nothing (aliases / older names).
            rows = _search_query(
                f'gene:{normalized} AND organism_id:{organism_id}',
                size=size,
            )
    except requests.RequestException as exc:
        raise GeneSearchError(
            'Could not reach UniProt to look up this gene. Try again shortly.'
        ) from exc

    return annotate_am_availability(rows)
