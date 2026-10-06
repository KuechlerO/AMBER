"""HTTP helpers + UniProt / Ensembl fetch with version identity."""

from __future__ import annotations

import hashlib
import json
import os
import time
from dataclasses import asdict, dataclass
from typing import Any

import requests


def http_get(url, *, headers=None, params=None, timeout=60, retries=None):
    retries = retries if retries is not None else int(os.environ.get('ENSEMBL_HTTP_RETRIES', '6'))
    last_exc: Exception | None = None
    for attempt in range(max(1, retries)):
        try:
            r = requests.get(url, headers=headers, params=params, timeout=timeout)
            if r.status_code in (429, 500, 502, 503, 504):
                last_exc = requests.HTTPError(f'{r.status_code} for {url}', response=r)
                time.sleep(min(30, 2.0 * (attempt + 1)))
                continue
            return r
        except (requests.Timeout, requests.ConnectionError) as exc:
            last_exc = exc
            if attempt + 1 >= retries:
                break
            time.sleep(min(30, 1.5 * (attempt + 1)))
    if last_exc is not None:
        raise last_exc
    raise RuntimeError(f'HTTP GET failed for {url}')


@dataclass
class SequenceIdentity:
    uniprot_accession: str
    uniprot_sequence_version: int
    uniprot_entry_version: int
    transcript_id: str  # versioned ENST…
    gene_symbol: str
    cds: str
    cds_sha256: str
    exon_regions: list[list[int]]  # JSON-serializable [[start, end], ...]

    def identity_key(self) -> str:
        return f'{self.uniprot_accession}|sv{self.uniprot_sequence_version}|{self.transcript_id}'

    def to_cache_dict(self) -> dict[str, Any]:
        d = asdict(self)
        return d

    @classmethod
    def from_cache_dict(cls, d: dict[str, Any]) -> 'SequenceIdentity':
        return cls(
            uniprot_accession=str(d['uniprot_accession']),
            uniprot_sequence_version=int(d['uniprot_sequence_version']),
            uniprot_entry_version=int(d['uniprot_entry_version']),
            transcript_id=str(d['transcript_id']),
            gene_symbol=str(d.get('gene_symbol') or ''),
            cds=str(d['cds']),
            cds_sha256=str(d['cds_sha256']),
            exon_regions=[list(x) for x in d['exon_regions']],
        )


def _extract_gene_symbol(data: dict) -> str:
    for gene in data.get('genes') or []:
        name = (gene.get('geneName') or {}).get('value')
        if name:
            return str(name)
    return ''


def _ensembl_transcript_ids_versioned(data: dict) -> list[str]:
    """ENST IDs from UniProt, keeping version suffix when present."""
    out: list[str] = []
    for ref in data.get('uniProtKBCrossReferences') or []:
        if ref.get('database') != 'Ensembl':
            continue
        tid = str(ref.get('id') or '').strip()
        if tid.startswith('ENST') and tid not in out:
            out.append(tid)
    return out


def _ensembl_gene_ids(data: dict) -> list[str]:
    out: list[str] = []

    def _add(raw: str) -> None:
        gid = str(raw or '').split('.')[0]
        if gid.startswith('ENSG') and gid not in out:
            out.append(gid)

    for ref in data.get('uniProtKBCrossReferences') or []:
        _add(ref.get('id') or '')
        if ref.get('database') == 'Ensembl':
            for prop in ref.get('properties') or []:
                if prop.get('key') == 'GeneId':
                    _add(prop.get('value') or '')
    return out


def _lookup_canonical_transcript(*, gene_symbol: str = '', gene_id: str = '') -> str | None:
    if gene_id:
        url = f'https://rest.ensembl.org/lookup/id/{gene_id}'
    elif gene_symbol:
        url = f'https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_symbol}'
    else:
        return None
    r = http_get(url, headers={'Content-Type': 'application/json'}, timeout=60)
    if r.status_code != 200:
        return None
    payload = r.json() or {}
    canonical = str(payload.get('canonical_transcript') or '')
    return canonical if canonical.startswith('ENST') else None


def resolve_transcript_id(data: dict) -> str:
    direct = _ensembl_transcript_ids_versioned(data)
    if direct:
        return direct[0]
    gene_symbol = _extract_gene_symbol(data)
    for gid in _ensembl_gene_ids(data):
        tid = _lookup_canonical_transcript(gene_id=gid)
        if tid:
            return tid
    if gene_symbol:
        tid = _lookup_canonical_transcript(gene_symbol=gene_symbol)
        if tid:
            return tid
    raise ValueError('Could not resolve Ensembl transcript for UniProt entry')


def ensure_versioned_transcript(transcript_id: str) -> str:
    """If ID has no version, resolve current version via Ensembl lookup."""
    tid = transcript_id.strip()
    if '.' in tid:
        return tid
    url = f'https://rest.ensembl.org/lookup/id/{tid}'
    r = http_get(url, headers={'Content-Type': 'application/json'}, timeout=60)
    if r.status_code != 200:
        return tid
    payload = r.json() or {}
    ver = payload.get('version')
    if ver is not None:
        return f'{tid.split(".")[0]}.{ver}'
    return tid


def fetch_uniprot_json(uniprot_id: str) -> dict:
    url = f'https://rest.uniprot.org/uniprotkb/{uniprot_id}.json'
    r = http_get(url, timeout=30)
    r.raise_for_status()
    return r.json()


def transcript_base_id(transcript_id: str) -> str:
    """ENST without version suffix (Ensembl sequence/map endpoints prefer this)."""
    return str(transcript_id).strip().split('.')[0]


def fetch_cds(transcript_id: str) -> str:
    base = transcript_base_id(transcript_id)
    url = f'https://rest.ensembl.org/sequence/id/{base}?type=cds'
    r = http_get(
        url,
        headers={'Content-Type': 'text/plain', 'Accept': 'text/plain'},
        timeout=60,
    )
    if r.status_code != 200:
        raise ValueError(f'CDS not found for {transcript_id} (queried {base}): HTTP {r.status_code}')
    return r.text.strip().upper()


def fetch_exon_regions(transcript_id: str, cds_length: int) -> list[list[int]]:
    base = transcript_base_id(transcript_id)
    url = f'https://rest.ensembl.org/map/cds/{base}/1..{cds_length}'
    r = http_get(url, headers={'Content-Type': 'application/json'}, timeout=60)
    r.raise_for_status()
    data = r.json()
    mappings = data.get('mappings') or []
    if not mappings:
        raise ValueError(f'No CDS mappings for {transcript_id}')
    strand = mappings[0]['strand']
    if strand == 1:
        mappings = sorted(mappings, key=lambda x: x['start'])
    else:
        mappings = sorted(mappings, key=lambda x: x['start'], reverse=True)
    exon_regions: list[list[int]] = []
    current = 0
    for block in mappings:
        genomic_len = abs(block['end'] - block['start']) + 1
        exon_regions.append([current, current + genomic_len])
        current += genomic_len
    return exon_regions


def fetch_uniprot_metadata(uniprot_accession: str) -> dict[str, Any]:
    """UniProt versions + transcript id + gene (no CDS yet)."""
    acc = uniprot_accession.strip().upper()
    data = fetch_uniprot_json(acc)
    audit = data.get('entryAudit') or {}
    transcript_id = ensure_versioned_transcript(resolve_transcript_id(data))
    return {
        'uniprot_accession': acc,
        'uniprot_sequence_version': int(audit.get('sequenceVersion') or 0),
        'uniprot_entry_version': int(audit.get('entryVersion') or 0),
        'gene_symbol': _extract_gene_symbol(data),
        'transcript_id': transcript_id,
    }


def complete_identity_with_cds(meta: dict[str, Any]) -> SequenceIdentity:
    """Fetch CDS + exons and build SequenceIdentity from metadata dict."""
    transcript_id = str(meta['transcript_id'])
    cds = fetch_cds(transcript_id)
    if not cds or len(cds) < 3:
        raise ValueError(f'Empty CDS for {meta["uniprot_accession"]} / {transcript_id}')
    exon_regions = fetch_exon_regions(transcript_id, len(cds))
    cds_sha = hashlib.sha256(cds.encode('ascii')).hexdigest()
    return SequenceIdentity(
        uniprot_accession=str(meta['uniprot_accession']),
        uniprot_sequence_version=int(meta['uniprot_sequence_version']),
        uniprot_entry_version=int(meta['uniprot_entry_version']),
        transcript_id=transcript_id,
        gene_symbol=str(meta.get('gene_symbol') or ''),
        cds=cds,
        cds_sha256=cds_sha,
        exon_regions=exon_regions,
    )


def load_sequence_identity(uniprot_accession: str) -> SequenceIdentity:
    """Fetch UniProt + Ensembl and build a versioned identity + CDS/exons."""
    return complete_identity_with_cds(fetch_uniprot_metadata(uniprot_accession))


def save_cds_cache(path, identity: SequenceIdentity) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.partial')
    payload = identity.to_cache_dict()
    with open(tmp, 'w', encoding='utf-8') as f:
        json.dump(payload, f)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def load_cds_cache(path) -> SequenceIdentity | None:
    if not path.is_file():
        return None
    try:
        with open(path, encoding='utf-8') as f:
            return SequenceIdentity.from_cache_dict(json.load(f))
    except (OSError, json.JSONDecodeError, KeyError, TypeError, ValueError):
        return None
