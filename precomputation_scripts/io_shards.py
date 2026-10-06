"""Versioned paths, atomic Parquet writes, append-only manifest."""

from __future__ import annotations

import json
import os
import re
import threading
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

from .uniprot_ensembl import SequenceIdentity

_MANIFEST_LOCK = threading.Lock()

# P04439_sv2_tvENST00000396221.7.parquet
_SHARD_RE = re.compile(
    r'^(?P<acc>[A-Z0-9]+)_sv(?P<sv>\d+)_tv(?P<tid>ENST\d+(?:\.\d+)?)\.parquet$'
)


def ensure_out_dirs(out_dir: Path) -> dict[str, Path]:
    out_dir = Path(out_dir)
    guides = out_dir / 'guides'
    cds_cache = out_dir / 'cds_cache'
    guides.mkdir(parents=True, exist_ok=True)
    cds_cache.mkdir(parents=True, exist_ok=True)
    return {
        'out': out_dir,
        'guides': guides,
        'cds_cache': cds_cache,
        'manifest': out_dir / 'manifest.jsonl',
    }


def safe_transcript_token(transcript_id: str) -> str:
    """Filename-safe transcript id (keep dots)."""
    return str(transcript_id).strip().replace('/', '_')


def shard_basename(accession: str, sequence_version: int, transcript_id: str) -> str:
    return (
        f'{accession.upper()}_sv{int(sequence_version)}'
        f'_tv{safe_transcript_token(transcript_id)}.parquet'
    )


def guide_shard_path(guides_dir: Path, identity: SequenceIdentity) -> Path:
    return Path(guides_dir) / shard_basename(
        identity.uniprot_accession,
        identity.uniprot_sequence_version,
        identity.transcript_id,
    )


def cds_cache_path(cds_dir: Path, identity: SequenceIdentity) -> Path:
    base = shard_basename(
        identity.uniprot_accession,
        identity.uniprot_sequence_version,
        identity.transcript_id,
    )
    return Path(cds_dir) / base.replace('.parquet', '.json')


def cds_cache_path_from_parts(cds_dir: Path, accession: str, seq_ver: int, transcript_id: str) -> Path:
    base = shard_basename(accession, seq_ver, transcript_id)
    return Path(cds_dir) / base.replace('.parquet', '.json')


def cleanup_partials(guides_dir: Path) -> int:
    """Remove orphaned *.parquet.partial files from interrupted runs."""
    n = 0
    for p in Path(guides_dir).glob('*.parquet.partial'):
        try:
            p.unlink()
            n += 1
        except OSError:
            pass
    return n


def scan_done_shards(guides_dir: Path) -> set[str]:
    """Return identity keys from final parquet filenames."""
    done: set[str] = set()
    for p in Path(guides_dir).glob('*.parquet'):
        if p.name.endswith('.partial'):
            continue
        m = _SHARD_RE.match(p.name)
        if not m:
            continue
        tid = m.group('tid')
        done.add(f"{m.group('acc')}|sv{m.group('sv')}|{tid}")
    return done


def parse_shard_name(name: str) -> dict[str, Any] | None:
    m = _SHARD_RE.match(name)
    if not m:
        return None
    return {
        'uniprot_accession': m.group('acc'),
        'uniprot_sequence_version': int(m.group('sv')),
        'transcript_id': m.group('tid'),
    }


def append_manifest(manifest_path: Path, record: dict[str, Any]) -> None:
    """Append one JSON line (thread-safe)."""
    record = dict(record)
    record.setdefault('ts', datetime.now(timezone.utc).isoformat())
    line = json.dumps(record, ensure_ascii=False, default=str) + '\n'
    with _MANIFEST_LOCK:
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        with open(manifest_path, 'a', encoding='utf-8') as f:
            f.write(line)
            f.flush()
            os.fsync(f.fileno())


def load_manifest_latest(manifest_path: Path) -> dict[str, dict[str, Any]]:
    """Accession → latest manifest record (by file order)."""
    latest: dict[str, dict[str, Any]] = {}
    if not manifest_path.is_file():
        return latest
    with open(manifest_path, encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            acc = str(rec.get('uniprot_accession') or '').upper()
            if acc:
                latest[acc] = rec
    return latest


def error_accessions(manifest_path: Path) -> list[str]:
    latest = load_manifest_latest(manifest_path)
    return sorted(a for a, r in latest.items() if r.get('status') == 'error')


def write_guides_parquet_atomic(path: Path, rows: list[dict[str, Any]]) -> int:
    """Write guide rows via .partial then os.replace. Returns row count."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = Path(str(path) + '.partial')
    if partial.exists():
        partial.unlink()
    df = pd.DataFrame(rows)
    # Lists → JSON strings for parquet stability
    if 'target_positions' in df.columns:
        df['target_positions'] = df['target_positions'].apply(
            lambda v: json.dumps(v) if not isinstance(v, str) else v
        )
    df.to_parquet(partial, index=False)
    # ensure durable
    with open(partial, 'rb') as f:
        os.fsync(f.fileno())
    os.replace(partial, path)
    return len(df)


def enrich_rows_with_identity(
    rows: list[dict[str, Any]],
    identity: SequenceIdentity,
    *,
    computed_at: str | None = None,
) -> list[dict[str, Any]]:
    ts = computed_at or datetime.now(timezone.utc).isoformat()
    out = []
    for r in rows:
        item = dict(r)
        item['uniprot_accession'] = identity.uniprot_accession
        item['uniprot_sequence_version'] = identity.uniprot_sequence_version
        item['uniprot_entry_version'] = identity.uniprot_entry_version
        item['transcript_id'] = identity.transcript_id
        item['cds_sha256'] = identity.cds_sha256
        item['gene_symbol'] = identity.gene_symbol
        item['computed_at'] = ts
        out.append(item)
    return out


def is_shard_complete(path: Path) -> bool:
    return Path(path).is_file() and Path(path).stat().st_size > 0
