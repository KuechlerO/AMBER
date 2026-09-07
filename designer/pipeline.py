import itertools
import time
from pathlib import Path
import pandas as pd
from Bio.Seq import Seq
from Bio.Data import CodonTable
import requests
from django.conf import settings
from django.db import connections

from .models import Alpha_missense
from .esm1b import lookup_esm1b
from .cadd import attach_cadd_to_outcomes, cadd_raw_for_best, fetch_cadd_scores
from .score_thresholds import (
    am_likely_pathogenic,
    cadd_significant,
    esm1b_damaging,
)

# ===========================
# simple help functions
# ===========================

def _http_get(url, *, headers=None, params=None, timeout=None, retries=None):
    """
    GET with retries for flaky upstreams (Ensembl via institutional proxy).
    CDS fetches often need >20s under load.
    """
    timeout = timeout if timeout is not None else int(
        getattr(settings, 'ENSEMBL_TIMEOUT_SEC', 60)
    )
    retries = retries if retries is not None else int(
        __import__('os').environ.get('ENSEMBL_HTTP_RETRIES', '3')
    )
    last_exc: Exception | None = None
    for attempt in range(max(1, retries)):
        try:
            r = requests.get(url, headers=headers, params=params, timeout=timeout)
            return r
        except (requests.Timeout, requests.ConnectionError) as exc:
            last_exc = exc
            if attempt + 1 >= retries:
                break
            time.sleep(1.5 * (attempt + 1))
    assert last_exc is not None
    raise last_exc


# ---------------------------
# Codon translation helpers
# ---------------------------
standard_table = CodonTable.unambiguous_dna_by_name['Standard']

def translate_codon(codon: str) -> str:
    """Translate a codon to amino acid; return 'X' on failure."""
    try:
        return str(Seq(codon).translate(table=standard_table))
    except Exception:
        return 'X'


# ---------------------------
# Gives back sgRNA changed base
# ---------------------------
def complementary_base(base):
    if(base == 'A'):
        return 'G'
    elif(base == 'T'):
        return 'C'
    elif(base == 'G'):
        return 'A'
    else: #base == "C"
        return 'T'


def wc_complement(base: str) -> str:
    """Watson–Crick complement (not the ABE/CBE substitution map)."""
    return {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}.get(str(base).upper(), base)


def normalize_cadd_chrom(name: str) -> str:
    text = str(name or '').strip()
    if text.lower().startswith('chr'):
        text = text[3:]
    return text or str(name)



# ===========================
#  UniProtID Input -> CDS, AlphaMissense results
# ===========================

# ---------------------------
# UniProtID/Ensembl -> CDS, transcriptID
# ---------------------------
def fetch_uniprot_json(uniprot_id):
    url = f'https://rest.uniprot.org/uniprotkb/{uniprot_id}.json'
    r = _http_get(url, timeout=30)
    r.raise_for_status()
    return r.json()


def extract_ensembl_transcript(data):
    for ref in data.get('uniProtKBCrossReferences', []):
        if ref['database'] == 'Ensembl':
            return ref['id'].split('.')[0]
    raise ValueError('Ensembl Transcipt not found!')


def fetch_cds(transcript_id):
    url = f'https://rest.ensembl.org/sequence/id/{transcript_id}?type=cds'
    r = _http_get(url, headers={'Content-Type': 'text/plain', 'Accept': 'text/plain'}, timeout=60)

    if r.status_code != 200:
        raise ValueError(f'CDS not found: {r.status_code}')

    return r.text.strip()

def get_ensembl_transcript_from_uniprot(uniprot_id):
    data = fetch_uniprot_json(uniprot_id)
    return extract_ensembl_transcript(data)

def get_cds_from_uniprot(uniprot_id):
    data = fetch_uniprot_json(uniprot_id)
    transcript_id = extract_ensembl_transcript(data)
    return fetch_cds(transcript_id)

def extract_gene_symbol(data) -> str | None:
    for gene in data.get('genes') or []:
        name = (gene.get('geneName') or {}).get('value')
        if name:
            return str(name)
    return None


def get_cds_and_transcript_from_uniprot(uniprot_id):
    data = fetch_uniprot_json(uniprot_id)
    transcript_id = extract_ensembl_transcript(data)
    cds = fetch_cds(transcript_id)
    return cds, transcript_id, extract_gene_symbol(data)

def ensembl_to_uniprot(ensembl_id: str) -> str | None:
    url = 'https://rest.uniprot.org/uniprotkb/search'

    params = {
        'query': ensembl_id,
        'format': 'json',
        'fields': 'accession'
    }

    try:
        r = _http_get(url, params=params, timeout=30)
        r.raise_for_status()
        data = r.json()

        results = data.get('results', [])
        if not results:
            return None

        print(results)
        # takes first UniProt ID
        return results[0]['primaryAccession']

    except requests.RequestException:
        return None

# ---------------------------
# UniProtID -> AlpaMissense results
# ---------------------------
def get_alphamissense_from_db(uniprot_nr, patho_threshold):

    results = Alpha_missense.objects.using('legacy_db').filter(
        uniprot_id=uniprot_nr
    ).values_list('protein_variant', 'am_pathogenicity')

    df = pd.DataFrame.from_records(
        results,
        columns=['protein_variant', 'am_pathogenicity']
    )

    df[['a.a.1', 'position', 'a.a.2']] = df['protein_variant'].str.extract(r'([A-Z])(\d+)([A-Z])')
    df = df.rename(columns={'am_pathogenicity': 'pathogenicity score'})

    return df

# ===========================
#  Exon/Intron filtering
# ===========================

# ---------------------------
# Ensemle-Mapping of exon boundaries
# ---------------------------
def parse_cds_mappings(mappings):
    """Build CDS exon regions plus genomic blocks from an Ensembl /map/cds payload."""
    if not mappings:
        return {
            'exon_regions': [],
            'exon_boundaries': [],
            'genomic_blocks': [],
        }

    strand = mappings[0]['strand']
    if strand == 1:
        mappings = sorted(mappings, key=lambda x: x['start'])
    else:
        mappings = sorted(mappings, key=lambda x: x['start'], reverse=True)

    exon_regions = []
    genomic_blocks = []
    current_cds_pos = 0

    for block in mappings:
        genomic_len = abs(block['end'] - block['start']) + 1
        cds_start = current_cds_pos
        cds_end = current_cds_pos + genomic_len
        exon_regions.append((cds_start, cds_end))
        genomic_blocks.append({
            'cds_start': cds_start,
            'cds_end': cds_end,
            'chrom': normalize_cadd_chrom(
                block.get('seq_region_name') or block.get('chromosome') or ''
            ),
            'start': int(block['start']),
            'end': int(block['end']),
            'strand': int(block['strand']),
        })
        current_cds_pos = cds_end

    exon_boundaries = [end for (_, end) in exon_regions[:-1]]
    return {
        'exon_regions': exon_regions,
        'exon_boundaries': exon_boundaries,
        'genomic_blocks': genomic_blocks,
    }


def get_exon_boundaries(transcript_id, cds_length):
    url = f'https://rest.ensembl.org/map/cds/{transcript_id}/1..{cds_length}'
    r = _http_get(url, headers={'Content-Type': 'application/json'}, timeout=60)
    r.raise_for_status()
    data = r.json()
    return parse_cds_mappings(data.get('mappings', []))


def cds_index_to_genomic(cds_pos, genomic_blocks):
    """Map a 0-based CDS index to (chrom, genomic_pos_1based, strand) or None."""
    try:
        idx = int(cds_pos)
    except (TypeError, ValueError):
        return None
    for block in genomic_blocks or []:
        if block['cds_start'] <= idx < block['cds_end']:
            offset = idx - block['cds_start']
            if int(block['strand']) == 1:
                genomic_pos = int(block['start']) + offset
            else:
                genomic_pos = int(block['end']) - offset
            return block['chrom'], genomic_pos, int(block['strand'])
    return None


def coding_edit_to_genomic_snv(cds_seq, cds_pos, alt_coding_base, genomic_blocks):
    """Coding-strand edit → genomic SNV (chrom, pos, ref, alt) for CADD."""
    mapped = cds_index_to_genomic(cds_pos, genomic_blocks)
    if mapped is None:
        return None
    chrom, genomic_pos, strand = mapped
    cds = str(cds_seq).upper()
    if cds_pos < 0 or cds_pos >= len(cds):
        return None
    ref_coding = cds[cds_pos]
    alt_coding = str(alt_coding_base).upper()
    if strand == 1:
        ref, alt = ref_coding, alt_coding
    else:
        ref, alt = wc_complement(ref_coding), wc_complement(alt_coding)
    if not chrom or ref == alt:
        return None
    return (chrom, int(genomic_pos), ref, alt)

# ---------------------------
# cuts guides over exon boundaries
# ---------------------------
def guide_is_within_single_exon(start, end, exon_regions):
    return any(start >= exon_start and end <= exon_end for exon_start, exon_end in exon_regions)



# ===========================
# sgRNA SEARCH
# ===========================

def find_sgrnas(cds_seq, codon_start, codon_end, base_sense, base_antisense, transcript_id, exon_regions, window_min, window_max, pam_type='NGG'):
    s = str(cds_seq).upper()
    L = len(s)
    removed_count = 0
    results = []

    # + strand
    for i in range(0, L - 23 + 1):
        protospacer = s[i:i + 20]
        pam = s[i + 20:i + 23]

        if len(pam) != 3:
            continue

        if pam_type == 'NGG':
            if not pam.endswith('GG'):
                continue
        elif pam_type == 'NG':
            if pam[1] != 'G':
                continue
        elif pam_type == 'NAG':
            if pam[1] != 'A' or pam[2] != 'G':
                continue

        guide_len = 23
        if not guide_is_within_single_exon(i, (i + 23), exon_regions):
            removed_count += 1
            continue

        window = protospacer[window_min-1:window_max]

        for j in range(20):
            pos_1b = j + 1
            # for editing window 4-8 -> 4 <= pos_1b_plus <= 8
            if window_min <= pos_1b <= window_max and protospacer[j] == base_sense:
                genome_idx = i + j
                if codon_start <= genome_idx < codon_end:
                    codon_idx = genome_idx - codon_start
                    key = ('+', i)
                    existing = next((r for r in results if r.get('_key') == key), None)

                    if existing is None:
                        entry = {
                            '_key': key,
                            'sgRNA_seq': protospacer + pam,
                            'protospacer_20nt': protospacer,
                            'pam': pam,
                            'strand': '+',
                            'target_A_pos': pos_1b,
                            'num_As_window': window.count(base_sense),
                            'codon_edit_indices': {codon_idx},
                            'edit_on_plus': {codon_idx: complementary_base(base_sense)},
                            '_target_positions': {pos_1b},
                        }
                        results.append(entry)
                    else:
                        existing['codon_edit_indices'].add(codon_idx)
                        existing['edit_on_plus'][codon_idx] = complementary_base(base_sense)
                        existing['_target_positions'].add(pos_1b)

    # - strand
    for i in range(3, L - 20):
        protospacer = s[i:i + 20]
        pam_upstream = s[i - 3:i]

        if len(pam_upstream) != 3:
            continue

        if pam_type == 'NGG':
            if not pam_upstream.startswith('CC'):
                continue
        elif pam_type == 'NG':
            if pam_upstream[1] != 'C':
                continue
        elif pam_type == 'NAG':
            if pam_upstream[0] != 'C' or pam_upstream[1] != 'T':
                continue

        if not guide_is_within_single_exon(i - 3, 23, exon_regions):
            removed_count += 1
            continue

        for j in range(20):
            pos_1b_plus = j + 1
            # for editing window 4-8 -> 13 <= pos_1b_plus <= 17
            if (20-window_max+1) <= pos_1b_plus <= (20-window_min+1) and protospacer[j] == base_antisense:
                genome_idx = i + j
                if codon_start <= genome_idx < codon_end:
                    codon_idx = genome_idx - codon_start
                    protospacer_rc = str(Seq(protospacer).reverse_complement())
                    pam_rc = str(Seq(pam_upstream).reverse_complement())
                    target_pos_rc = 20 - j
                    key = ('-', i)
                    existing = next((r for r in results if r.get('_key') == key), None)

                    if existing is None:
                        window_rc = protospacer_rc[window_min-1:window_max]
                        entry = {
                            '_key': key,
                            'sgRNA_seq': protospacer_rc + pam_rc,
                            'protospacer_20nt': protospacer_rc,
                            'pam': pam_rc,
                            'strand': '-',
                            'target_A_pos': target_pos_rc,
                            'num_As_window': window_rc.count(base_sense),
                            'codon_edit_indices': {codon_idx},
                            'edit_on_plus': {codon_idx: complementary_base(base_antisense)},
                            '_target_positions': {target_pos_rc},
                        }
                        results.append(entry)
                    else:
                        existing['codon_edit_indices'].add(codon_idx)
                        existing['edit_on_plus'][codon_idx] = complementary_base(base_antisense)
                        existing['_target_positions'].add(target_pos_rc)
    for r in results:
        r['target_A_pos'] = sorted(r['_target_positions'])
        del r['_target_positions']
        del r['_key']
        r['codon_edit_indices'] = sorted(r['codon_edit_indices'])

    unique = {}
    for r in results:
        key = r['sgRNA_seq']
        if key not in unique:
            unique[key] = r
        else:
            unique[key]['codon_edit_indices'] = sorted(
                set(unique[key]['codon_edit_indices']).union(r['codon_edit_indices'])
            )
            unique[key]['edit_on_plus'].update(r['edit_on_plus'])

    return list(unique.values())



# ===========================
# Enumerate codon outcomes
# ===========================
def codon_outcomes_for_sgRNA(wt_codon, codon_edit_indices, edit_on_plus, base_sense, base_antisense):
    outcomes = []
    valid_positions = []
    for pos in codon_edit_indices:
        desired_base = edit_on_plus.get(pos)
        if desired_base == complementary_base(base_sense) and wt_codon[pos] == base_sense:
            valid_positions.append(pos)
        elif desired_base == complementary_base(base_antisense) and wt_codon[pos] == base_antisense:
            valid_positions.append(pos)
    for r in range(1, len(valid_positions) + 1):
        for subset in itertools.combinations(valid_positions, r):
            m = list(wt_codon)
            for pos in subset:
                m[pos] = edit_on_plus[pos]
            new_codon = ''.join(m)
            outcomes.append({
                'codon': new_codon,
                'aa': translate_codon(new_codon),
                'edited_indices': list(subset),
            })
    seen = set()
    uniq = []
    for o in outcomes:
        if o['codon'] not in seen:
            uniq.append(o)
            seen.add(o['codon'])
    return uniq


# ===========================
# Candidate generation
# ===========================
def _outcome_snvs(cds_seq, codon_start, edited_indices, edit_on_plus, genomic_blocks):
    snvs = []
    if not genomic_blocks:
        return snvs
    for idx in edited_indices or []:
        alt = edit_on_plus.get(idx)
        if not alt:
            continue
        snv = coding_edit_to_genomic_snv(cds_seq, codon_start + idx, alt, genomic_blocks)
        if snv:
            snvs.append(snv)
    return snvs


def generate_candidates(
    cds_seq,
    patho_df,
    a_or_c,
    patho_threshold,
    transcript_id,
    window_min,
    window_max,
    pam_type='NGG',
    *,
    exon_data=None,
    uniprot_id=None,
    esm_lookup=None,
):
    cds = Seq(str(cds_seq).upper())

    # ABE: A -> G und T -> C
    # CBE: G -> A und C -> T
    if a_or_c == 'ABE':
        b1 = 'A'
        b2 = 'T'
    elif a_or_c == 'CBE':
        b1 = 'C'
        b2 = 'G'
    else:
        raise ValueError(f'Unsupported editor mode: {a_or_c}')

    if exon_data is None:
        exon_data = get_exon_boundaries(transcript_id, len(cds))
    exon_regions = exon_data['exon_regions']
    genomic_blocks = exon_data.get('genomic_blocks') or []

    agg = (
        patho_df.groupby(['position', 'a.a.1'])
        .agg({'a.a.2': list, 'pathogenicity score': list})
        .reset_index()
        .rename(columns={'a.a.2': 'pathogenic_AA', 'pathogenicity score': 'pathogenicity_scores'})
    )

    out_rows = []

    for _, row in agg.iterrows():
        pos_aa = int(row['position']) - 1
        wt_aa = row['a.a.1']
        pathogenic_aas = list(row['pathogenic_AA'])
        pathogenic_scores = list(row['pathogenicity_scores'])
        pat_map = {aa: sc for aa, sc in zip(pathogenic_aas, pathogenic_scores)}

        codon_start = pos_aa * 3
        codon_end = codon_start + 3
        codon_seq = str(cds[codon_start:codon_end])

        if b1 not in codon_seq and b2 not in codon_seq:
            continue

        sgrnas = find_sgrnas(
            cds,
            codon_start,
            codon_end,
            b1,
            b2,
            transcript_id,
            exon_regions,
            window_min,
            window_max,
            pam_type=pam_type
        )

        sg_cols = []
        any_pathogenic = False
        best_aa = None
        best_score = None
        avg_alpha_score = None

        for sg in sgrnas:

            outcomes = codon_outcomes_for_sgRNA(
                codon_seq,
                sg['codon_edit_indices'],
                sg['edit_on_plus'],
                b1,
                b2
            )

            annotated = []
            sg_is_pathogenic = False
            local_best_score = None
            local_best_aa = None

            for o in outcomes:
                aa = o['aa']
                cod = o['codon']
                score = pat_map.get(aa, None)

                is_pat = aa in pat_map
                esm_score = None
                if esm_lookup is not None and uniprot_id and aa:
                    esm_score = esm_lookup(uniprot_id, pos_aa + 1, aa, wt_aa)
                snvs = _outcome_snvs(
                    cds,
                    codon_start,
                    o.get('edited_indices') or [],
                    sg['edit_on_plus'],
                    genomic_blocks,
                )

                annotated.append({
                    'codon': cod,
                    'aa': aa,
                    'score': score,
                    'is_pathogenic': is_pat,
                    'esm1b': esm_score,
                    'snvs': snvs,
                })

                if is_pat:
                    sg_is_pathogenic = True
                    any_pathogenic = True

                    if score is not None and (best_score is None or score > best_score):
                        best_score = score
                        best_aa = aa

                    if score is not None and (local_best_score is None or score > local_best_score):
                        local_best_score = score
                        local_best_aa = aa

            # filter only sgRNAs with min one pathogenic score
            if not sg_is_pathogenic:
                continue

            annotated.sort(
                key=lambda x: (
                    not x['is_pathogenic'],
                    -(x['score'] if x['score'] is not None else -1),
                    x['codon']
                )
            )

            all_scores = [x['score'] for x in annotated if x['score'] is not None]
            avg_alpha_score = sum(all_scores) / len(all_scores) if all_scores else None

            sg_cols.append({
                'seq': sg['sgRNA_seq'],
                'prot': sg['protospacer_20nt'],
                'pam': sg['pam'],
                'strand': sg['strand'],
                'tpos': sg['target_A_pos'],
                'nAwin': sg['num_As_window'],
                'outcome_details': annotated,
                'local_best_aa': local_best_aa,
                'local_best_score': local_best_score,
                'avg_alpha_score': avg_alpha_score,
            })

        row_dict = {
            'position': pos_aa + 1,
            'WT_codon': codon_seq,
            'WT_AA': wt_aa,
            'pathogenic_AA': ','.join(pathogenic_aas),
            f'{a_or_c}_is_pathogenic': bool(any_pathogenic),
            f'{a_or_c}_mutated_AA': best_aa,
            f'{a_or_c}_mutated_AA_score': best_score,
            'sgRNA_found': bool(len(sg_cols) > 0),
            'avg_alpha_score': avg_alpha_score,
        }

        seen_seqs, idx = set(), 1
        for sg in sg_cols:
            if sg['seq'] in seen_seqs:
                continue
            seen_seqs.add(sg['seq'])

            row_dict[f'sgRNA_{idx}_seq'] = sg['seq']
            row_dict[f'sgRNA_{idx}_protospacer'] = sg['prot']
            row_dict[f'sgRNA_{idx}_pam'] = sg['pam']
            row_dict[f'sgRNA_{idx}_strand'] = sg['strand']
            row_dict[f'sgRNA_{idx}_TargetApos'] = sg['tpos']
            row_dict[f'sgRNA_{idx}_numAsWindow'] = sg['nAwin']
            row_dict[f'sgRNA_{idx}_outcome_details'] = sg.get('outcome_details') or []
            row_dict[f'sgRNA_{idx}_local_mut_aa'] = sg.get('local_best_aa')
            row_dict[f'sgRNA_{idx}_local_am_score'] = sg.get('local_best_score')
            row_dict[f'sgRNA_{idx}_avg_alpha_score'] = sg.get('avg_alpha_score')
            idx += 1

        out_rows.append(row_dict)

    df = pd.DataFrame(out_rows)

    if not df.empty:
        df = df.sort_values(
            by=[f'{a_or_c}_is_pathogenic', 'sgRNA_found'],
            ascending=[False, False],
            kind='stable'
        ).reset_index(drop=True)

        df = df.where(pd.notnull(df), None)

    return df



def kill_doubles(df):
    """Legacy: dedupe wide dataframe rows by first sgRNA column (prefer apply_duplicate_mode)."""
    if df is None or df.empty:
        return df

    seq_col = None
    for col in df.columns:
        if col.endswith('_seq') and 'sgRNA' in col:
            seq_col = col
            break

    if seq_col is None:
        raise ValueError('No sgRNA sequence column found')

    df = df.drop_duplicates(subset=[seq_col], keep='first')

    return df


def normalize_duplicate_mode(mode):
    """Map form values (and legacy hide) to unique | best | group."""
    if mode in ('unique', 'best', 'group'):
        return mode
    if mode == 'hide':
        return 'best'
    return 'best'


def annotate_sgrna_binding(rows):
    """Add binding_count and binding_positions from the full candidate list."""
    from collections import Counter, defaultdict

    seq_counts = Counter(r.get('sgrna_seq') for r in rows if r.get('sgrna_seq'))
    seq_positions = defaultdict(set)
    for row in rows:
        seq = row.get('sgrna_seq')
        pos = row.get('position')
        if seq and pos is not None:
            seq_positions[seq].add(pos)

    for row in rows:
        seq = row.get('sgrna_seq')
        row['binding_positions'] = sorted(seq_positions.get(seq, []))
        row['binding_count'] = seq_counts.get(seq, 0)
    return rows


def _row_rank_for_best(row):
    """Higher alpha_score wins; tie-break by lower protein position."""
    score = row.get('alpha_score')
    try:
        score_val = float(score) if score is not None else -1.0
    except (TypeError, ValueError):
        score_val = -1.0
    pos = row.get('position')
    pos_val = int(pos) if pos is not None else 999999
    return (score_val, -pos_val)


def apply_duplicate_mode(rows, duplicate_mode):
    """
    Filter/sort flattened guide rows.

    unique — only guides that map to a single position.
    best   — one row per guide at the best position (highest AlphaMissense; tie → lower position).
    group  — all rows; same guide sequence grouped together.
    """
    mode = normalize_duplicate_mode(duplicate_mode)
    rows = annotate_sgrna_binding(list(rows))

    if mode == 'unique':
        filtered = [r for r in rows if r.get('binding_count', 0) == 1]
        filtered.sort(key=lambda r: r.get('position') if r.get('position') is not None else 999999)
        return filtered

    if mode == 'best':
        best_by_seq = {}
        for row in rows:
            seq = row.get('sgrna_seq')
            if not seq:
                continue
            if seq not in best_by_seq or _row_rank_for_best(row) > _row_rank_for_best(best_by_seq[seq]):
                best_by_seq[seq] = row
        result = list(best_by_seq.values())
        result.sort(key=lambda r: r.get('position') if r.get('position') is not None else 999999)
        return result

    if mode == 'group':
        return sort_rows_by_sgrna_group(rows)

    return rows


def sort_rows_by_sgrna_group(rows):
    """Place rows with the same guide RNA sequence adjacent (stable within groups)."""
    seq_first_rank = {}
    for i, row in enumerate(rows):
        seq = row.get('sgrna_seq')
        if seq and seq not in seq_first_rank:
            seq_first_rank[seq] = i

    return [
        row for _, row in sorted(
            enumerate(rows),
            key=lambda item: (
                seq_first_rank.get(item[1].get('sgrna_seq'), item[0]),
                item[0],
            ),
        )
    ]


def _none_if_na(value):
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    return value


def _list_or_empty(value):
    cleaned = _none_if_na(value)
    return cleaned if isinstance(cleaned, list) else []


def _format_score_line(codon, aa, score, *, digits=3, star=False):
    if score is None:
        return f"{codon}({aa})"
    marker = '*' if star else ''
    return f"{codon}({aa})[{score:.{digits}f}]{marker}"


def finalize_guide_scores(rows, score_map):
    """Attach CADD to each outcome and build per-guide lists + max/avg scalars."""
    for row in rows:
        details = []
        for item in row.pop('_outcome_details', None) or []:
            details.append(dict(item))
        attach_cadd_to_outcomes(details, score_map)

        am_scores = [d['score'] for d in details if d.get('score') is not None]
        esm_scores = [d['esm1b'] for d in details if d.get('esm1b') is not None]
        cadd_scores = [d['cadd'] for d in details if d.get('cadd') is not None]

        am_best = None
        for detail in details:
            score = detail.get('score')
            if score is None:
                continue
            if am_best is None or score > am_best['score']:
                am_best = detail
        if am_best is not None:
            row['alpha_score'] = am_best['score']
            row['mut_aa'] = am_best['aa']

        row['outcomes'] = [
            _format_score_line(
                d.get('codon'), d.get('aa'), d.get('score'),
                digits=3, star=am_likely_pathogenic(d.get('score')),
            )
            for d in details
        ]
        row['esm1b_outcomes'] = [
            _format_score_line(
                d.get('codon'), d.get('aa'), d.get('esm1b'),
                digits=3, star=esm1b_damaging(d.get('esm1b')),
            )
            for d in details
        ]
        row['cadd_outcomes'] = [
            _format_score_line(
                d.get('codon'), d.get('aa'), d.get('cadd'),
                digits=2, star=cadd_significant(d.get('cadd')),
            )
            for d in details
        ]
        row['avg_alpha_score'] = sum(am_scores) / len(am_scores) if am_scores else None
        row['esm1b_score'] = min(esm_scores) if esm_scores else None
        row['avg_esm1b_score'] = sum(esm_scores) / len(esm_scores) if esm_scores else None
        row['cadd_phred'] = max(cadd_scores) if cadd_scores else None
        row['avg_cadd_phred'] = sum(cadd_scores) / len(cadd_scores) if cadd_scores else None
        row['alpha_score_significant'] = am_likely_pathogenic(row.get('alpha_score'))
        row['esm1b_score_significant'] = esm1b_damaging(row.get('esm1b_score'))
        row['cadd_phred_significant'] = cadd_significant(row.get('cadd_phred'))

        max_cadd_snvs = []
        best_cadd = None
        for detail in details:
            cadd = detail.get('cadd')
            if cadd is None:
                continue
            if best_cadd is None or cadd > best_cadd:
                best_cadd = cadd
                max_cadd_snvs = detail.get('snvs') or []
        row['cadd_raw'] = cadd_raw_for_best(max_cadd_snvs, score_map)
    return rows


def flatten_candidate_dataframe(results_df, editor_mode, top_sgrnas):
    try:
        n = int(top_sgrnas)
    except (TypeError, ValueError):
        n = 5

    if results_df is None or results_df.empty:
        return []

    rows = []
    score_col = f'{editor_mode}_mutated_AA_score'
    aa_col = f'{editor_mode}_mutated_AA'

    for _, row in results_df.iterrows():
        sg_entries = []
        idx = 1

        while f'sgRNA_{idx}_seq' in results_df.columns:
            seq = row.get(f'sgRNA_{idx}_seq')
            if seq:
                local_aa = _none_if_na(row.get(f'sgRNA_{idx}_local_mut_aa'))
                local_am = _none_if_na(row.get(f'sgRNA_{idx}_local_am_score'))
                sg_entries.append({
                    'position': row.get('position'),
                    'wt_codon': row.get('WT_codon'),
                    'wt_aa': row.get('WT_AA'),
                    'mut_aa': local_aa if local_aa is not None else row.get(aa_col),
                    'alpha_score': local_am if local_am is not None else row.get(score_col),
                    'sgrna_seq': row.get(f'sgRNA_{idx}_seq'),
                    'protospacer': row.get(f'sgRNA_{idx}_protospacer'),
                    'pam': row.get(f'sgRNA_{idx}_pam'),
                    'strand': row.get(f'sgRNA_{idx}_strand'),
                    'target_position': row.get(f'sgRNA_{idx}_TargetApos'),
                    'editor_used': editor_mode,
                    'num_as_window': row.get(f'sgRNA_{idx}_numAsWindow'),
                    'avg_alpha_score': _none_if_na(row.get(f'sgRNA_{idx}_avg_alpha_score')),
                    '_outcome_details': _list_or_empty(row.get(f'sgRNA_{idx}_outcome_details')),
                })
            idx += 1

        sg_entries.sort(
            key=lambda g: (
                0 if g['target_position'] in (5, 6) else 1,
                g['num_as_window'] if g['num_as_window'] is not None else 999,
            )
        )

        rows.extend(sg_entries[:n])

    return rows


# ===========================
# run whole pipeline
# ===========================
def run_pipeline(uniprot_id, editor, alpha_threshold, top_sgrnas, window_min, window_max, pam_type='NGG'):

    transcript_id = ''
    cds = ''
    gene_symbol = None
    input_was_ensembl = uniprot_id.startswith('ENST')

    if input_was_ensembl:
        transcript_id = uniprot_id
        uniprot_id = ensembl_to_uniprot(transcript_id)
        if not uniprot_id:
            raise ValueError(f'Could not map Ensembl ID {transcript_id} to a UniProt accession.')
        cds = fetch_cds(transcript_id)
        try:
            gene_symbol = extract_gene_symbol(fetch_uniprot_json(uniprot_id))
        except Exception:
            gene_symbol = None
    else:
        cds, transcript_id, gene_symbol = get_cds_and_transcript_from_uniprot(uniprot_id)
        print(transcript_id)

    threshold = float(alpha_threshold)
    patho_df = get_alphamissense_from_db(uniprot_id, threshold)
    exon_data = get_exon_boundaries(transcript_id, len(cds))

    all_rows = []
    gen_kwargs = dict(
        exon_data=exon_data,
        uniprot_id=uniprot_id,
        esm_lookup=lookup_esm1b,
    )

    if editor in ('ABE', 'BOTH'):
        df_abe = generate_candidates(
            cds, patho_df, 'ABE', threshold, transcript_id, window_min, window_max,
            pam_type=pam_type, **gen_kwargs,
        )
        all_rows.extend(flatten_candidate_dataframe(df_abe, 'ABE', top_sgrnas))

    if editor in ('CBE', 'BOTH'):
        df_cbe = generate_candidates(
            cds, patho_df, 'CBE', threshold, transcript_id, window_min, window_max,
            pam_type=pam_type, **gen_kwargs,
        )
        all_rows.extend(flatten_candidate_dataframe(df_cbe, 'CBE', top_sgrnas))

    # Filter out incomplete rows
    all_rows = [
        r for r in all_rows
        if r.get('sgrna_seq')
        and pd.notna(r.get('pam'))
        and pd.notna(r.get('strand'))
        and pd.notna(r.get('alpha_score'))
    ]

    cadd_snvs = []
    for row in all_rows:
        for detail in row.get('_outcome_details') or []:
            cadd_snvs.extend(detail.get('snvs') or [])
    cadd_map, cadd_warning = fetch_cadd_scores(cadd_snvs)
    finalize_guide_scores(all_rows, cadd_map)

    all_rows.sort(
        key=lambda r: (
            r['position'] if r['position'] is not None else 999999,
        )
    )

    protein_length = len(cds) // 3
    return {
        'guide_rows': annotate_sgrna_binding(all_rows),
        'no_guide_positions': [],
        'protein_length': protein_length,
        'uniprot_accession': uniprot_id,
        'gene_symbol': gene_symbol or '',
        'transcript_id': transcript_id,
        'cadd_warning': cadd_warning,
    }
