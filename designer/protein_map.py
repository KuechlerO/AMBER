"""Protein map helpers: domain resolution for any UniProt accession."""

from __future__ import annotations

import json
import threading
from pathlib import Path
from typing import Any

from django.conf import settings

from .pipeline import fetch_uniprot_json
from .screen_plot.data_loader import get_screen_data_store
from .screen_plot.gene_lookup import lookup_gene_by_uniprot, normalize_uniprot_accession

_DOMAIN_TYPES = frozenset({
    'Domain',
    'Region',
    'Motif',
    'Repeat',
    'Zinc finger',
    'Topological domain',
    'Transmembrane',
    'Intramembrane',
    'Coiled coil',
    'Compositional bias',
})

_cache_lock = threading.Lock()
_domain_cache: dict[str, list[dict[str, Any]]] = {}


def _screen_data_root() -> Path:
    return Path(getattr(settings, 'SCREEN_DATA_DIR', settings.BASE_DIR / 'files-archive-dir'))


def domains_from_uniprot_features(uniprot_data: dict) -> list[dict[str, Any]]:
    """Convert UniProt JSON features into overview domain layout entries."""
    layout: list[dict[str, Any]] = []
    for feat in uniprot_data.get('features') or []:
        typ = feat.get('type') or ''
        if typ not in _DOMAIN_TYPES:
            continue
        loc = feat.get('location') or {}
        start = (loc.get('start') or {}).get('value')
        end = (loc.get('end') or {}).get('value')
        if start is None or end is None:
            continue
        try:
            start_i, end_i = int(start), int(end)
        except (TypeError, ValueError):
            continue
        if end_i < start_i:
            continue
        name = feat.get('description') or typ
        layout.append({
            'name': str(name),
            'start': start_i,
            'end': end_i,
            'avg_lfc': None,
            'domain_type': typ,
        })
    return layout


def domains_from_local_json(gene: str) -> list[dict[str, Any]] | None:
    """Return structural domains from local domain_data.json, or None if missing."""
    if not gene:
        return None
    store = get_screen_data_store()
    try:
        store.ensure_domains_only()
    except Exception:
        domain_path = _screen_data_root() / 'libraries' / 'domain_data.json'
        if not domain_path.is_file():
            return None
        with domain_path.open(encoding='utf-8') as fh:
            domains = json.load(fh)
        raw = domains.get(gene)
        if not raw:
            return None
        layout = []
        for domain in raw:
            layout.append({
                'name': str(domain['description']),
                'start': int(domain['location']['start']['value']),
                'end': int(domain['location']['end']['value']),
                'avg_lfc': None,
                'domain_type': domain.get('type', ''),
            })
        return layout
    layout = store.get_structural_domain_layout(gene)
    return layout if layout else None


def resolve_domain_layout(
    *,
    gene: str | None,
    uniprot_accession: str | None,
) -> tuple[list[dict[str, Any]], str]:
    """
    Prefer local screen-library domains; otherwise fetch UniProt features.

    Returns (layout, source) where source is 'local', 'uniprot', or 'none'.
    """
    if gene:
        local = domains_from_local_json(gene)
        if local:
            return local, 'local'

    accession = normalize_uniprot_accession(uniprot_accession or '')
    if not accession:
        return [], 'none'

    with _cache_lock:
        if accession in _domain_cache:
            return _domain_cache[accession], 'uniprot'

    try:
        data = fetch_uniprot_json(accession)
        layout = domains_from_uniprot_features(data)
    except Exception:
        return [], 'none'

    with _cache_lock:
        _domain_cache[accession] = layout
    return layout, ('uniprot' if layout else 'none')


def resolve_plot_gene_label(form_data: dict) -> str:
    """Best display/plot label: gene symbol, else accession, else input id."""
    symbol = (form_data.get('gene_symbol') or '').strip()
    if symbol:
        return symbol
    meta = lookup_gene_by_uniprot(
        form_data.get('uniprot_accession') or form_data.get('uniprot_id') or ''
    )
    if meta:
        return meta['gene']
    accession = (form_data.get('uniprot_accession') or '').strip()
    if accession:
        return accession
    return (form_data.get('uniprot_id') or 'protein').strip() or 'protein'


def alphafold_pdb_path(uniprot_accession: str) -> Path | None:
    """Local AF2 PDB if present for this accession."""
    accession = normalize_uniprot_accession(uniprot_accession)
    if not accession:
        return None
    path = _screen_data_root() / 'libraries' / 'AF_screen_proteins' / f'AF-{accession}-F1-model_v4.pdb'
    if path.is_file():
        return path
    return None


def alphafold_pdb_url(uniprot_accession: str, version: int | None = None) -> str:
    accession = normalize_uniprot_accession(uniprot_accession)
    ver = int(version) if version else 4
    return f'https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v{ver}.pdb'


def resolve_alphafold_pdb_url(uniprot_accession: str) -> str | None:
    """
    Resolve a working AFDB PDB URL for an accession.

    Prefers the API `latestVersion`, then falls back through recent versions.
    """
    import requests as http_requests

    accession = normalize_uniprot_accession(uniprot_accession)
    if not accession:
        return None

    versions: list[int] = []
    try:
        api = http_requests.get(
            f'https://alphafold.ebi.ac.uk/api/prediction/{accession}',
            timeout=20,
        )
        if api.status_code == 200:
            payload = api.json()
            if isinstance(payload, list) and payload:
                latest = payload[0].get('latestVersion')
                all_versions = payload[0].get('allVersions') or []
                if latest:
                    versions.append(int(latest))
                for v in all_versions:
                    try:
                        vi = int(v)
                    except (TypeError, ValueError):
                        continue
                    if vi not in versions:
                        versions.append(vi)
    except Exception:
        pass

    for v in (6, 5, 4, 3, 2, 1):
        if v not in versions:
            versions.append(v)

    for ver in versions:
        url = alphafold_pdb_url(accession, version=ver)
        try:
            head = http_requests.head(url, timeout=20, allow_redirects=True)
            if head.status_code == 200:
                return url
        except Exception:
            continue
    return None
