"""Quick validation of versioned shards for a few UniProt accessions."""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

from .io_shards import ensure_out_dirs, parse_shard_name


def check_shards(out_dir: Path, accessions: list[str] | None = None) -> int:
    paths = ensure_out_dirs(out_dir)
    guides = list(paths['guides'].glob('*.parquet'))
    if accessions:
        wanted = {a.upper() for a in accessions}
        guides = [p for p in guides if parse_shard_name(p.name) and parse_shard_name(p.name)['uniprot_accession'] in wanted]

    if not guides:
        print('No matching parquet shards found.', file=sys.stderr)
        return 1

    print(f'Checking {len(guides)} shard(s) under {paths["guides"]}')
    required = {
        'uniprot_accession',
        'uniprot_sequence_version',
        'uniprot_entry_version',
        'transcript_id',
        'cds_sha256',
        'editor_used',
        'pam_type',
        'sgrna_seq',
        'alpha_score',
        'target_positions',
        'outcome_details',
        'computed_at',
    }
    errors = 0
    for path in sorted(guides):
        meta = parse_shard_name(path.name)
        df = pd.read_parquet(path)
        missing = required - set(df.columns)
        counts = Counter(zip(df['editor_used'], df['pam_type'])) if len(df) else Counter()
        print(f'\n{path.name}')
        print(f'  filename meta: {meta}')
        print(f'  rows: {len(df)}')
        if len(df):
            print(f'  accession: {df["uniprot_accession"].iloc[0]}  '
                  f'sv={df["uniprot_sequence_version"].iloc[0]}  '
                  f'tid={df["transcript_id"].iloc[0]}')
            print(f'  cds_sha256: {str(df["cds_sha256"].iloc[0])[:16]}…')
            print(f'  editor×pam: {dict(counts)}')
        if missing:
            print(f'  MISSING columns: {sorted(missing)}')
            errors += 1
        if meta and len(df):
            if str(df['uniprot_accession'].iloc[0]).upper() != meta['uniprot_accession']:
                print('  ERROR: accession mismatch vs filename')
                errors += 1
            if int(df['uniprot_sequence_version'].iloc[0]) != meta['uniprot_sequence_version']:
                print('  ERROR: sequence_version mismatch vs filename')
                errors += 1
            if str(df['transcript_id'].iloc[0]) != meta['transcript_id']:
                print('  ERROR: transcript_id mismatch vs filename')
                errors += 1
        # spot-check outcome_details JSON
        if len(df) and 'outcome_details' in df.columns:
            sample = df['outcome_details'].iloc[0]
            try:
                json.loads(sample) if isinstance(sample, str) else sample
            except json.JSONDecodeError:
                print('  ERROR: outcome_details not valid JSON')
                errors += 1

    print(f'\nResult: {"OK" if errors == 0 else f"{errors} issue(s)"}')
    return 0 if errors == 0 else 1


def main(argv: list[str] | None = None) -> None:
    root = Path(__file__).resolve().parents[1]
    p = argparse.ArgumentParser(description='Validate precomputed guide shards')
    p.add_argument(
        '--out',
        type=Path,
        default=root / 'data' / 'precomputed_guides',
        help='Output directory used by the precompute CLI',
    )
    p.add_argument(
        '--uniprot',
        type=str,
        default='',
        help='Optional comma-separated accessions to filter',
    )
    args = p.parse_args(argv)
    accessions = None
    if args.uniprot:
        accessions = [x.strip() for x in args.uniprot.split(',') if x.strip()]
    sys.exit(check_shards(args.out, accessions))


if __name__ == '__main__':
    main()
