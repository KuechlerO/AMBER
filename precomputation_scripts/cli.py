"""CLI for proteome guide precomputation (resumable)."""

from __future__ import annotations

import argparse
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

from . import WINDOW_MAX, WINDOW_MIN
from .am_db import fetch_variants, list_accessions
from .design import design_all_guides
from .uniprot_ensembl import (
    SequenceIdentity,
    complete_identity_with_cds,
    fetch_uniprot_metadata,
    load_cds_cache,
    save_cds_cache,
)
from .io_shards import (
    append_manifest,
    cds_cache_path_from_parts,
    cleanup_partials,
    enrich_rows_with_identity,
    ensure_out_dirs,
    error_accessions,
    guide_shard_path,
    is_shard_complete,
    write_guides_parquet_atomic,
)


def _resolve_identity(accession: str, cds_cache_dir: Path, force_refresh: bool) -> SequenceIdentity:
    """Prefer versioned CDS cache; only fetch Ensembl CDS/exons on miss or --force."""
    meta = fetch_uniprot_metadata(accession)
    cache_path = cds_cache_path_from_parts(
        cds_cache_dir,
        meta['uniprot_accession'],
        meta['uniprot_sequence_version'],
        meta['transcript_id'],
    )
    if not force_refresh:
        cached = load_cds_cache(cache_path)
        if (
            cached is not None
            and cached.uniprot_accession == meta['uniprot_accession']
            and cached.uniprot_sequence_version == meta['uniprot_sequence_version']
            and cached.transcript_id == meta['transcript_id']
            and cached.cds
            and cached.exon_regions
        ):
            # Keep entry version / gene from fresh UniProt metadata
            cached.uniprot_entry_version = meta['uniprot_entry_version']
            if meta.get('gene_symbol'):
                cached.gene_symbol = meta['gene_symbol']
            return cached
    identity = complete_identity_with_cds(meta)
    save_cds_cache(cache_path, identity)
    return identity


def process_one_accession(
    accession: str,
    out_dir: str,
    *,
    force: bool = False,
    dry_run: bool = False,
) -> dict:
    """Worker entry: design guides for one UniProt accession. Safe to call in a subprocess."""
    paths = ensure_out_dirs(Path(out_dir))
    acc = accession.strip().upper()
    try:
        identity = _resolve_identity(acc, paths['cds_cache'], force_refresh=force)
        shard = guide_shard_path(paths['guides'], identity)
        if is_shard_complete(shard) and not force:
            rec = {
                'uniprot_accession': acc,
                'uniprot_sequence_version': identity.uniprot_sequence_version,
                'uniprot_entry_version': identity.uniprot_entry_version,
                'transcript_id': identity.transcript_id,
                'cds_sha256': identity.cds_sha256,
                'gene_symbol': identity.gene_symbol,
                'status': 'skipped',
                'n_guides': None,
                'shard': str(shard.name),
                'reason': 'shard_exists',
            }
            append_manifest(paths['manifest'], rec)
            return rec

        if dry_run:
            rec = {
                'uniprot_accession': acc,
                'uniprot_sequence_version': identity.uniprot_sequence_version,
                'uniprot_entry_version': identity.uniprot_entry_version,
                'transcript_id': identity.transcript_id,
                'cds_sha256': identity.cds_sha256,
                'gene_symbol': identity.gene_symbol,
                'status': 'dry_run',
                'n_guides': None,
                'shard': str(shard.name),
            }
            append_manifest(paths['manifest'], rec)
            return rec

        am_df = fetch_variants(acc)
        if am_df.empty:
            raise ValueError(f'No AlphaMissense variants for {acc}')

        rows = design_all_guides(
            identity.cds,
            am_df,
            identity.exon_regions,
            window_min=WINDOW_MIN,
            window_max=WINDOW_MAX,
        )
        ts = datetime.now(timezone.utc).isoformat()
        rows = enrich_rows_with_identity(rows, identity, computed_at=ts)
        n = write_guides_parquet_atomic(shard, rows)
        rec = {
            'uniprot_accession': acc,
            'uniprot_sequence_version': identity.uniprot_sequence_version,
            'uniprot_entry_version': identity.uniprot_entry_version,
            'transcript_id': identity.transcript_id,
            'cds_sha256': identity.cds_sha256,
            'gene_symbol': identity.gene_symbol,
            'status': 'ok',
            'n_guides': n,
            'shard': str(shard.name),
            'window_min': WINDOW_MIN,
            'window_max': WINDOW_MAX,
        }
        append_manifest(paths['manifest'], rec)
        return rec
    except Exception as exc:  # noqa: BLE001
        rec = {
            'uniprot_accession': acc,
            'status': 'error',
            'error': str(exc),
            'traceback': traceback.format_exc()[-2000:],
        }
        try:
            append_manifest(paths['manifest'], rec)
        except OSError:
            pass
        return rec


def _build_worklist(args, paths) -> list[str]:
    if args.uniprot:
        ids = [x.strip().upper() for x in args.uniprot.replace(';', ',').split(',') if x.strip()]
        return ids
    if args.retry_errors:
        errs = error_accessions(paths['manifest'])
        if not errs:
            print('No error accessions in manifest.', flush=True)
        return errs
    print('Listing distinct UniProt IDs from AlphaMissense DB…', flush=True)
    return list_accessions()


def run(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            'Precompute AMBER guides for AlphaMissense proteins '
            f'(window {WINDOW_MIN}-{WINDOW_MAX}, ABE+CBE, PAM NGG/NG/NAG). '
            'Resumable: re-run the same command after interrupt.'
        )
    )
    root = Path(__file__).resolve().parents[1]
    parser.add_argument(
        '--out',
        type=Path,
        default=root / 'data' / 'precomputed_guides',
        help='Output directory (default: data/precomputed_guides)',
    )
    parser.add_argument('--workers', type=int, default=4, help='Parallel workers (default 4)')
    parser.add_argument(
        '--uniprot',
        type=str,
        default='',
        help='Comma-separated UniProt accessions (pilot). If omitted, all AM accessions.',
    )
    parser.add_argument(
        '--force',
        action='store_true',
        help='Rebuild even if versioned shard already exists',
    )
    parser.add_argument(
        '--retry-errors',
        action='store_true',
        help='Only re-process accessions whose latest manifest status is error',
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Resolve versions / paths only; do not design or write guides',
    )
    parser.add_argument(
        '--limit',
        type=int,
        default=0,
        help='Process at most N accessions (0 = all)',
    )
    args = parser.parse_args(argv)

    paths = ensure_out_dirs(args.out)
    n_partial = cleanup_partials(paths['guides'])
    if n_partial:
        print(f'Removed {n_partial} orphaned .parquet.partial file(s).', flush=True)

    work = _build_worklist(args, paths)
    if args.limit and args.limit > 0:
        work = work[: args.limit]
    if not work:
        print('Nothing to do.', flush=True)
        return 0

    print(
        f'Queue: {len(work)} accession(s); workers={args.workers}; out={paths["out"]}',
        flush=True,
    )

    ok = err = skipped = 0
    workers = max(1, int(args.workers))

    if workers == 1:
        results_iter = (
            process_one_accession(
                acc, str(paths['out']), force=args.force, dry_run=args.dry_run
            )
            for acc in work
        )
        for i, rec in enumerate(results_iter, 1):
            status = rec.get('status')
            if status == 'ok':
                ok += 1
            elif status == 'skipped':
                skipped += 1
            elif status == 'error':
                err += 1
            print(
                f'[{i}/{len(work)}] {rec.get("uniprot_accession")} → {status}'
                + (f' n={rec.get("n_guides")}' if rec.get('n_guides') is not None else '')
                + (f' err={rec.get("error")}' if status == 'error' else ''),
                flush=True,
            )
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futs = {
                pool.submit(
                    process_one_accession,
                    acc,
                    str(paths['out']),
                    force=args.force,
                    dry_run=args.dry_run,
                ): acc
                for acc in work
            }
            done_n = 0
            for fut in as_completed(futs):
                done_n += 1
                acc = futs[fut]
                try:
                    rec = fut.result()
                except Exception as exc:  # noqa: BLE001
                    rec = {'uniprot_accession': acc, 'status': 'error', 'error': str(exc)}
                    append_manifest(paths['manifest'], rec)
                status = rec.get('status')
                if status == 'ok':
                    ok += 1
                elif status == 'skipped':
                    skipped += 1
                elif status == 'error':
                    err += 1
                print(
                    f'[{done_n}/{len(work)}] {rec.get("uniprot_accession")} → {status}'
                    + (f' n={rec.get("n_guides")}' if rec.get('n_guides') is not None else '')
                    + (f' err={rec.get("error")}' if status == 'error' else ''),
                    flush=True,
                )

    print(f'Done. ok={ok} skipped={skipped} error={err} manifest={paths["manifest"]}', flush=True)
    return 1 if err and not ok and not skipped else 0


def main(argv: list[str] | None = None) -> None:
    sys.exit(run(argv))


if __name__ == '__main__':
    main()
