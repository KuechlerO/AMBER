# Precomputation scripts — proteome AMBER guides

Batch tools to precompute **all** base-editing guide RNAs for every UniProt accession
that has AlphaMissense scores in Postgres, and store them as **versioned Parquet shards**.

These scripts are **standalone** (no Django). They do **not** change the live AMBER web app.

| Parameter | Value |
|-----------|--------|
| Editing window | **1–10** (widest; filter narrower windows later) |
| AM threshold | **None** (all AM-mapped outcomes) |
| Editors | ABE + CBE |
| PAM | NGG, NG, NAG |
| ESM1b / CADD | Not computed |

---

## Prerequisites

1. Conda env with the project stack (e.g. `crispr_tool_env`):
   - `pandas`, `pyarrow` (or `fastparquet`), `biopython`, `requests`, `psycopg2`
2. Network access to **UniProt** and **Ensembl** REST APIs
3. Postgres AlphaMissense DB (same as AMBER `legacy_db`), via env vars:

| Variable | Default |
|----------|---------|
| `AMBER_DB_NAME` | `swpdb` |
| `AMBER_DB_USER` | `swp2026` |
| `AMBER_DB_PASSWORD` | `swp2026` |
| `AMBER_DB_HOST` | tried as `127.0.0.1` then `medgen_pgsql_db` (or set explicitly) |
| `AMBER_DB_PORT` | `5432` |
| `ENSEMBL_HTTP_RETRIES` | `6` |

4. Disk: budget about **40–50 GB** free for a full proteome run (~20.5k proteins).
5. Run commands from the **repository root** so `python -m precomputation_scripts` resolves.

```bash
cd /path/to/SWP2026-CRISPR-Tool
conda activate crispr_tool_env
```

---

## Quick start (pilot)

Process two proteins with two workers:

```bash
python -m precomputation_scripts --uniprot P04439,Q9BXL7 --workers 2
```

Validate shards:

```bash
python -m precomputation_scripts.pilot_check --uniprot P04439,Q9BXL7
```

---

## Full proteome run

```bash
python -m precomputation_scripts --workers 4 --out data/precomputed_guides
```

On a machine with **8 CPUs / 15 GB RAM**, prefer `--workers 4` (or at most `6`).  
Higher worker counts increase Ensembl/UniProt rate-limit risk more than they speed up design.

Rough expectations (first cold run):

- **Wall time:** ~1–2 days (API-bound); faster once `cds_cache/` is warm  
- **Storage:** ~25–50 GB Parquet + CDS cache  

Optional: `--limit 100` to process only the first N accessions while testing.

---

## Resume after interrupt (default)

The job is **resume-safe**. If it is killed (Ctrl+C, SSH drop, OOM):

1. Re-run **the same command**.
2. Accessions that already have a final versioned `*.parquet` are **skipped**.
3. Orphaned `*.parquet.partial` files from mid-write are deleted at startup and those proteins are redone.

You do **not** need a special flag to resume. Resume is the default.

```bash
# same as before — continues where it left off
python -m precomputation_scripts --workers 4 --out data/precomputed_guides
```

### Force rebuild

Rebuild the **current** UniProt/Ensembl identity even if the shard exists:

```bash
python -m precomputation_scripts --uniprot P04439 --force --workers 1
```

If UniProt’s **sequence version** or the Ensembl transcript version changed, a **new** filename is used automatically; the old shard is left in place (no silent overwrite of another version).

### Retry failures only

Re-queue accessions whose **latest** `manifest.jsonl` status is `error`:

```bash
python -m precomputation_scripts --workers 4 --retry-errors
```

---

## Output layout

```
data/precomputed_guides/          # default --out
  manifest.jsonl                  # append-only log (ok / skipped / error / dry_run)
  cds_cache/
    P04439_sv2_tvENST00000396221.7.json
  guides/
    P04439_sv2_tvENST00000396221.7.parquet
```

### Versioned filenames

```
{ACCESSION}_sv{sequenceVersion}_tv{ENSEMBL_TRANSCRIPT_ID}.parquet
```

Example: `P04439_sv2_tvENST00000396221.7.parquet`

Each guide row also stores:

- `uniprot_accession`, `uniprot_sequence_version`, `uniprot_entry_version`
- `transcript_id` (with version), `cds_sha256`, `gene_symbol`, `computed_at`
- design fields: `editor_used`, `pam_type`, `position`, scores, spacer, `target_positions`, `outcome_details`, …

---

## Inspecting progress

```bash
# How many finished shards?
ls data/precomputed_guides/guides/*.parquet 2>/dev/null | wc -l

# Tail the manifest
tail -n 20 data/precomputed_guides/manifest.jsonl

# Count errors in latest status per accession (simple scan)
python -m precomputation_scripts.pilot_check
```

---

## Dry run

Resolve UniProt/Ensembl versions and paths without designing guides:

```bash
python -m precomputation_scripts --uniprot P04439 --dry-run --workers 1
```

---

## Troubleshooting

| Symptom | What to try |
|---------|-------------|
| Cannot connect to DB | Set `AMBER_DB_HOST=127.0.0.1` (SSH tunnel) or `medgen_pgsql_db` (Docker network). Check password env vars. |
| Many `429` / timeouts | Lower `--workers` to 2; increase `ENSEMBL_HTTP_RETRIES`; re-run (resume + `--retry-errors`). |
| `No AlphaMissense variants` | Accession not in AM table, or ID typo. |
| `Could not resolve Ensembl transcript` | UniProt entry lacks usable Ensembl xrefs; logged as `error` in manifest — skip or fix manually. |
| `CDS not found` / Ensembl 400 | Scripts query Ensembl with the **unversioned** ENST but still **store** the versioned ID in filenames/rows. Re-run with `--retry-errors` after API blips. |
| Empty / corrupt shard after crash | Ensure no `*.partial` left; delete bad `.parquet` if needed and re-run with `--force` for that ID. |
| `ModuleNotFoundError: pyarrow` | `conda install pyarrow` or `pip install pyarrow` in the env. |

---

## Module map

| Module | Role |
|--------|------|
| `cli.py` | Main entry (`python -m precomputation_scripts`) |
| `design.py` | Guide search + AM outcomes |
| `uniprot_ensembl.py` | Versioned UniProt / Ensembl + CDS cache |
| `am_db.py` | List accessions / fetch AM variants |
| `io_shards.py` | Paths, atomic Parquet write, manifest |
| `pilot_check.py` | Validate shard columns and version consistency |
