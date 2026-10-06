"""Guide RNA design for precompute (window 1–10, ABE/CBE × NGG/NG/NAG).

Logic matches AMBER designer/pipeline.py (PAM, exon filter, pathogenic outcomes).
Minus-strand exon span uses the full 23-nt guide (CDS start i-3 … i+20).
"""

from __future__ import annotations

import itertools
import json
from typing import Any

import pandas as pd
from Bio.Data import CodonTable
from Bio.Seq import Seq

from . import EDITORS, PAM_TYPES, WINDOW_MAX, WINDOW_MIN

standard_table = CodonTable.unambiguous_dna_by_name['Standard']


def translate_codon(codon: str) -> str:
    try:
        return str(Seq(codon).translate(table=standard_table))
    except Exception:
        return 'X'


def complementary_base(base: str) -> str:
    """Base-editor substitution map (not Watson–Crick)."""
    b = str(base).upper()
    if b == 'A':
        return 'G'
    if b == 'T':
        return 'C'
    if b == 'G':
        return 'A'
    return 'T'


def guide_is_within_single_exon(start, end, exon_regions) -> bool:
    return any(start >= exon_start and end <= exon_end for exon_start, exon_end in exon_regions)


def _pam_ok_plus(pam: str, pam_type: str) -> bool:
    if len(pam) != 3:
        return False
    if pam_type == 'NGG':
        return pam.endswith('GG')
    if pam_type == 'NG':
        return pam[1] == 'G'
    if pam_type == 'NAG':
        return pam[1] == 'A' and pam[2] == 'G'
    return False


def _pam_ok_minus(pam_upstream: str, pam_type: str) -> bool:
    if len(pam_upstream) != 3:
        return False
    if pam_type == 'NGG':
        return pam_upstream.startswith('CC')
    if pam_type == 'NG':
        return pam_upstream[1] == 'C'
    if pam_type == 'NAG':
        return pam_upstream[0] == 'C' and pam_upstream[1] == 'T'
    return False


def index_guides_by_codon(
    cds_seq,
    base_sense,
    base_antisense,
    exon_regions,
    window_min=WINDOW_MIN,
    window_max=WINDOW_MAX,
    pam_type='NGG',
):
    """Map 0-indexed codon → list of guide dicts that edit that codon."""
    s = str(cds_seq).upper()
    L = len(s)
    by_codon: dict[int, dict] = {}

    for i in range(0, L - 23 + 1):
        protospacer = s[i:i + 20]
        pam = s[i + 20:i + 23]
        if not _pam_ok_plus(pam, pam_type):
            continue
        if not guide_is_within_single_exon(i, i + 23, exon_regions):
            continue
        window = protospacer[window_min - 1:window_max]
        for j in range(20):
            pos_1b = j + 1
            if not (window_min <= pos_1b <= window_max and protospacer[j] == base_sense):
                continue
            genome_idx = i + j
            codon_idx = genome_idx % 3
            key = ('+', i)
            codon = genome_idx // 3
            bucket = by_codon.setdefault(codon, {})
            existing = bucket.get(key)
            if existing is None:
                bucket[key] = {
                    'sgRNA_seq': protospacer + pam,
                    'protospacer_20nt': protospacer,
                    'pam': pam,
                    'strand': '+',
                    'num_As_window': window.count(base_sense),
                    'codon_edit_indices': {codon_idx},
                    'edit_on_plus': {codon_idx: complementary_base(base_sense)},
                    '_target_positions': {pos_1b},
                }
            else:
                existing['codon_edit_indices'].add(codon_idx)
                existing['edit_on_plus'][codon_idx] = complementary_base(base_sense)
                existing['_target_positions'].add(pos_1b)

    for i in range(3, L - 20):
        protospacer = s[i:i + 20]
        pam_upstream = s[i - 3:i]
        if not _pam_ok_minus(pam_upstream, pam_type):
            continue
        if not guide_is_within_single_exon(i - 3, i + 20, exon_regions):
            continue
        for j in range(20):
            pos_1b_plus = j + 1
            if not (
                (20 - window_max + 1) <= pos_1b_plus <= (20 - window_min + 1)
                and protospacer[j] == base_antisense
            ):
                continue
            genome_idx = i + j
            codon_idx = genome_idx % 3
            protospacer_rc = str(Seq(protospacer).reverse_complement())
            pam_rc = str(Seq(pam_upstream).reverse_complement())
            target_pos_rc = 20 - j
            window_rc = protospacer_rc[window_min - 1:window_max]
            key = ('-', i)
            codon = genome_idx // 3
            bucket = by_codon.setdefault(codon, {})
            existing = bucket.get(key)
            if existing is None:
                bucket[key] = {
                    'sgRNA_seq': protospacer_rc + pam_rc,
                    'protospacer_20nt': protospacer_rc,
                    'pam': pam_rc,
                    'strand': '-',
                    'num_As_window': window_rc.count(base_sense),
                    'codon_edit_indices': {codon_idx},
                    'edit_on_plus': {codon_idx: complementary_base(base_antisense)},
                    '_target_positions': {target_pos_rc},
                }
            else:
                existing['codon_edit_indices'].add(codon_idx)
                existing['edit_on_plus'][codon_idx] = complementary_base(base_antisense)
                existing['_target_positions'].add(target_pos_rc)

    out: dict[int, list[dict[str, Any]]] = {}
    for codon, keyed in by_codon.items():
        unique: dict[str, dict[str, Any]] = {}
        for r in keyed.values():
            r['target_A_pos'] = sorted(r.pop('_target_positions'))
            r['codon_edit_indices'] = sorted(r['codon_edit_indices'])
            seq = r['sgRNA_seq']
            if seq not in unique:
                unique[seq] = r
            else:
                unique[seq]['codon_edit_indices'] = sorted(
                    set(unique[seq]['codon_edit_indices']).union(r['codon_edit_indices'])
                )
                unique[seq]['edit_on_plus'].update(r['edit_on_plus'])
                unique[seq]['target_A_pos'] = sorted(
                    set(unique[seq]['target_A_pos']).union(r['target_A_pos'])
                )
        out[codon] = list(unique.values())
    return out


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


def _score_guides_at_residue(cds, pos_aa, wt_aa, pat_map, sgrnas, b1, b2, editor, pam_type):
    codon_start = pos_aa * 3
    codon_seq = str(cds[codon_start:codon_start + 3])
    if len(codon_seq) < 3:
        return []
    if b1 not in codon_seq and b2 not in codon_seq:
        return []

    guide_rows = []
    for sg in sgrnas:
        outcomes = codon_outcomes_for_sgRNA(
            codon_seq, sg['codon_edit_indices'], sg['edit_on_plus'], b1, b2
        )
        annotated = []
        local_best = None
        local_aa = None
        has_am = False
        for o in outcomes:
            score = pat_map.get(o['aa'])
            is_pat = o['aa'] in pat_map
            detail = {
                'codon': o['codon'],
                'aa': o['aa'],
                'score': None if score is None else float(score),
                'is_pathogenic': is_pat,
            }
            annotated.append(detail)
            if is_pat:
                has_am = True
                if score is not None and (local_best is None or score > local_best):
                    local_best = float(score)
                    local_aa = o['aa']
        if not has_am:
            continue
        am_scores = [d['score'] for d in annotated if d.get('score') is not None]
        avg_alpha = sum(am_scores) / len(am_scores) if am_scores else None
        target_positions = list(sg.get('target_A_pos') or [])
        guide_rows.append({
            'position': pos_aa + 1,
            'wt_codon': codon_seq,
            'wt_aa': wt_aa,
            'mut_aa': local_aa,
            'alpha_score': local_best,
            'avg_alpha_score': avg_alpha,
            'sgrna_seq': sg['sgRNA_seq'],
            'protospacer': sg['protospacer_20nt'],
            'pam': sg['pam'],
            'strand': sg['strand'],
            'target_positions': target_positions,
            'num_editable_window': sg.get('num_As_window'),
            'editor_used': editor,
            'pam_type': pam_type,
            'outcome_details': json.dumps(annotated, separators=(',', ':')),
        })
    return guide_rows


def design_all_guides(
    cds: str,
    patho_df: pd.DataFrame,
    exon_regions: list,
    *,
    window_min: int = WINDOW_MIN,
    window_max: int = WINDOW_MAX,
    editors: tuple[str, ...] = EDITORS,
    pam_types: tuple[str, ...] = PAM_TYPES,
) -> list[dict[str, Any]]:
    """Design all ABE/CBE × PAM guides; no AM score threshold, no top-N cut."""
    cds_seq = Seq(str(cds).upper())
    regions = [tuple(r) for r in exon_regions]

    if patho_df is None or patho_df.empty:
        return []

    work = patho_df.copy()
    work['position'] = work['position'].astype(int)
    grouped = (
        work.groupby(['position', 'a.a.1'], as_index=False)
        .agg({'a.a.2': list, 'pathogenicity score': list})
    )

    # Pre-index per (editor, pam)
    indexes: dict[tuple[str, str], dict[int, list]] = {}
    for editor in editors:
        if editor == 'ABE':
            b1, b2 = 'A', 'T'
        elif editor == 'CBE':
            b1, b2 = 'C', 'G'
        else:
            continue
        for pam_type in pam_types:
            indexes[(editor, pam_type)] = index_guides_by_codon(
                cds_seq, b1, b2, regions, window_min, window_max, pam_type
            )

    all_rows: list[dict[str, Any]] = []
    for _, row in grouped.iterrows():
        pos_aa = int(row['position']) - 1
        wt_aa = row['a.a.1']
        pat_map = {aa: sc for aa, sc in zip(row['a.a.2'], row['pathogenicity score'])}
        for editor, pam_type in indexes:
            if editor == 'ABE':
                b1, b2 = 'A', 'T'
            else:
                b1, b2 = 'C', 'G'
            sgrnas = indexes[(editor, pam_type)].get(pos_aa, [])
            all_rows.extend(
                _score_guides_at_residue(
                    cds_seq, pos_aa, wt_aa, pat_map, sgrnas, b1, b2, editor, pam_type
                )
            )
    return all_rows
