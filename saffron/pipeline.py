"""SAFFRON analysis pipeline: SignalP + SP-region base-editing guides."""

from __future__ import annotations

import re
import uuid
from typing import Any

from Bio.Seq import Seq

from designer.pipeline import (
    codon_outcomes_for_sgRNA,
    ensembl_to_uniprot,
    extract_gene_symbol,
    fetch_cds,
    fetch_uniprot_json,
    find_sgrnas,
    get_cds_and_transcript_from_uniprot,
    get_exon_boundaries,
    translate_codon,
)
from designer.services import normalize_input_id

from .signalp_client import (
    SignalPPrediction,
    compute_deltas,
    prediction_to_dict,
    run_signalp,
    signalp_status,
)
from .patho_spv import annotate_guide_rows


def _with_backend(payload: dict[str, Any]) -> dict[str, Any]:
    payload = dict(payload)
    payload['signalp_backend'] = signalp_status()['backend']
    return payload

MUTATION_RE = re.compile(r'^([A-Z])(\d+)([A-Z])$')
AA_ALPHABET = set('ACDEFGHIKLMNPQRSTVWYUBZOX')


class SaffronUserError(ValueError):
    """User-facing validation / analysis error."""


def parse_mutation(raw: str | None) -> tuple[str, int, str] | None:
    text = (raw or '').strip().upper().replace(' ', '')
    if not text:
        return None
    m = MUTATION_RE.match(text)
    if not m:
        raise SaffronUserError('Mutation must look like A23V (WT amino acid, position, mutant amino acid).')
    return m.group(1), int(m.group(2)), m.group(3)


def normalize_aa_sequence(raw: str) -> str:
    seq = ''.join(c for c in (raw or '').upper() if c.isalpha())
    if len(seq) < 10:
        raise SaffronUserError('Amino-acid sequence must be at least 10 residues.')
    if len(seq) > 10000:
        raise SaffronUserError('Amino-acid sequence exceeds the 10,000 residue limit.')
    return seq


def extract_uniprot_sequence(data: dict) -> str:
    seq = ((data.get('sequence') or {}).get('value') or '').strip().upper()
    if not seq:
        raise SaffronUserError('UniProt entry has no protein sequence.')
    return ''.join(c for c in seq if c.isalpha())


def apply_mutation(sequence: str, wt: str, pos: int, mut: str) -> str:
    if pos < 1 or pos > len(sequence):
        raise SaffronUserError(f'Mutation position {pos} is outside the protein length ({len(sequence)}).')
    observed = sequence[pos - 1]
    if observed != wt:
        raise SaffronUserError(
            f'Mutation {wt}{pos}{mut} does not match the sequence (found {observed} at position {pos}).'
        )
    if wt == mut:
        raise SaffronUserError('Mutation must change the amino acid.')
    return sequence[: pos - 1] + mut + sequence[pos:]


def _protein_from_cds(cds: str) -> str:
    return str(Seq(str(cds).upper()).translate(to_stop=False)).replace('*', '')


def resolve_uniprot_inputs(uniprot_id: str) -> dict[str, Any]:
    uid = normalize_input_id(uniprot_id)
    gene_symbol = None
    transcript_id = ''
    cds = ''
    accession = uid

    if uid.startswith('ENST'):
        transcript_id = uid
        accession = ensembl_to_uniprot(transcript_id)
        if not accession:
            raise SaffronUserError(f'Could not map Ensembl ID {transcript_id} to a UniProt accession.')
        cds = fetch_cds(transcript_id)
        data = fetch_uniprot_json(accession)
        gene_symbol = extract_gene_symbol(data)
        aa_seq = extract_uniprot_sequence(data)
    else:
        data = fetch_uniprot_json(uid)
        accession = uid
        gene_symbol = extract_gene_symbol(data)
        aa_seq = extract_uniprot_sequence(data)
        cds, transcript_id, _ = get_cds_and_transcript_from_uniprot(uid)

    # Prefer UniProt AA sequence; fall back to CDS translation if needed
    if not aa_seq and cds:
        aa_seq = _protein_from_cds(cds)

    return {
        'accession': accession,
        'gene_symbol': gene_symbol or '',
        'transcript_id': transcript_id,
        'cds': cds,
        'aa_sequence': aa_seq,
        'protein_length': len(aa_seq),
    }


def _editor_bases(editor: str) -> tuple[str, str]:
    if editor == 'ABE':
        return 'A', 'T'
    if editor == 'CBE':
        return 'C', 'G'
    raise SaffronUserError(f'Unsupported editor: {editor}')


def _missense_outcomes(wt_codon: str, wt_aa: str, sg: dict, base_sense: str, base_antisense: str) -> list[dict]:
    outcomes = codon_outcomes_for_sgRNA(
        wt_codon,
        sg['codon_edit_indices'],
        sg['edit_on_plus'],
        base_sense,
        base_antisense,
    )
    return [
        o for o in outcomes
        if o.get('aa') and o['aa'] != wt_aa and o['aa'] not in ('*', 'X')
    ]


def design_sp_guides(
    cds: str,
    transcript_id: str,
    sp_start: int,
    sp_end: int,
    *,
    editor: str,
    window_min: int,
    window_max: int,
    pam_type: str,
    top_sgrnas: int,
) -> list[dict[str, Any]]:
    """Enumerate ABE/CBE guides that introduce missense edits inside [sp_start, sp_end] (1-based)."""
    cds_seq = Seq(str(cds).upper())
    protein_len = len(cds_seq) // 3
    start = max(1, int(sp_start))
    end = min(protein_len, int(sp_end))
    if start > end:
        return []

    exon_regions = get_exon_boundaries(transcript_id, len(cds_seq))['exon_regions']
    editors = ['ABE', 'CBE'] if editor == 'BOTH' else [editor]
    try:
        top_n = max(1, int(top_sgrnas))
    except (TypeError, ValueError):
        top_n = 5

    rows: list[dict[str, Any]] = []
    for ed in editors:
        b1, b2 = _editor_bases(ed)
        for pos in range(start, end + 1):
            codon_start = (pos - 1) * 3
            codon_end = codon_start + 3
            wt_codon = str(cds_seq[codon_start:codon_end])
            if len(wt_codon) != 3:
                continue
            wt_aa = translate_codon(wt_codon)
            if b1 not in wt_codon and b2 not in wt_codon:
                continue

            sgrnas = find_sgrnas(
                cds_seq,
                codon_start,
                codon_end,
                b1,
                b2,
                transcript_id,
                exon_regions,
                window_min,
                window_max,
                pam_type=pam_type,
            )

            scored: list[dict] = []
            for sg in sgrnas:
                missense = _missense_outcomes(wt_codon, wt_aa, sg, b1, b2)
                if not missense:
                    continue
                # Prefer single-base outcomes that change AA; take first unique mut AA entries
                for o in missense:
                    scored.append({
                        'position': pos,
                        'wt_codon': wt_codon,
                        'wt_aa': wt_aa,
                        'mut_aa': o['aa'],
                        'mut_codon': o['codon'],
                        'sgrna_seq': sg['sgRNA_seq'],
                        'protospacer': sg['protospacer_20nt'],
                        'pam': sg['pam'],
                        'strand': sg['strand'],
                        'target_position': sg.get('target_A_pos'),
                        'editor_used': ed,
                        'num_as_window': sg.get('num_As_window'),
                        'outcomes': [f"{oo['codon']}→{oo['aa']}" for oo in missense],
                    })

            # Prefer guides with target positions 5/6, then fewer editable bases in window
            scored.sort(
                key=lambda g: (
                    0 if (
                        isinstance(g['target_position'], list)
                        and any(t in (5, 6) for t in g['target_position'])
                    ) or g['target_position'] in (5, 6) else 1,
                    g['num_as_window'] if g['num_as_window'] is not None else 999,
                    g['mut_aa'],
                )
            )
            # Deduplicate by (guide, mut_aa) then take top_n guides (by guide seq first occurrence)
            seen_guides: set[str] = set()
            kept = []
            for g in scored:
                key = g['sgrna_seq']
                if key in seen_guides:
                    # still allow different mut_aa from same guide if not yet kept for that pair
                    if any(k['sgrna_seq'] == key and k['mut_aa'] == g['mut_aa'] for k in kept):
                        continue
                else:
                    if len(seen_guides) >= top_n:
                        continue
                    seen_guides.add(key)
                kept.append(g)
            rows.extend(kept)

    rows.sort(key=lambda r: (r['position'], r['editor_used'], r['mut_aa'], r['sgrna_seq']))
    return rows


def _mutant_protein(wt_seq: str, position: int, mut_aa: str) -> str:
    return wt_seq[: position - 1] + mut_aa + wt_seq[position:]


def _attach_signalp_to_rows(
    rows: list[dict],
    wt_seq: str,
    wt_pred: SignalPPrediction,
    mut_preds: dict[str, SignalPPrediction],
) -> list[dict]:
    enriched = []
    for row in rows:
        mid = row['mutant_id']
        mut_pred = mut_preds.get(mid)
        if mut_pred is None:
            continue
        deltas = compute_deltas(wt_pred, mut_pred)
        enriched.append({
            **row,
            'sp_prediction': mut_pred.prediction,
            'sp_prob': mut_pred.sp_prob,
            'cs_before': mut_pred.cs_before,
            'cs_after': mut_pred.cs_after,
            'cs_prob': mut_pred.cs_prob,
            'plot_path': mut_pred.plot_path,
            **deltas,
            'highlighted': bool(row.get('highlighted')),
        })
    return enriched


def run_sequence_analysis(
    aa_sequence: str,
    *,
    organism: str = 'eukarya',
    mutation: str | None = None,
) -> dict[str, Any]:
    seq = normalize_aa_sequence(aa_sequence)
    mut = parse_mutation(mutation)
    job_id = f'seq_{uuid.uuid4().hex[:12]}'

    sequences = {'WT': seq}
    focus = None
    if mut:
        wt, pos, mut_aa = mut
        mut_seq = apply_mutation(seq, wt, pos, mut_aa)
        sequences['FOCUS'] = mut_seq
        focus = {'wt_aa': wt, 'position': pos, 'mut_aa': mut_aa, 'label': f'{wt}{pos}{mut_aa}'}

    preds = run_signalp(sequences, organism=organism, job_id=job_id, want_plots=True)
    wt_pred = preds.get('WT')
    if not wt_pred:
        raise SaffronUserError('SignalP did not return a WT prediction.')

    focus_row = None
    if focus and 'FOCUS' in preds:
        deltas = compute_deltas(wt_pred, preds['FOCUS'])
        focus_row = {
            'position': focus['position'],
            'wt_aa': focus['wt_aa'],
            'mut_aa': focus['mut_aa'],
            'label': focus['label'],
            'sp_prediction': preds['FOCUS'].prediction,
            'sp_prob': preds['FOCUS'].sp_prob,
            'cs_before': preds['FOCUS'].cs_before,
            'cs_prob': preds['FOCUS'].cs_prob,
            'plot_path': preds['FOCUS'].plot_path,
            **deltas,
        }

    return _with_backend({
        'mode': 'sequence',
        'aa_sequence': seq,
        'protein_length': len(seq),
        'organism': organism,
        'wt_signalp': prediction_to_dict(wt_pred),
        'focus_mutation': focus_row,
        'guide_rows': [],
        'guides_available': False,
        'message': 'Guide design requires a UniProt or Ensembl ID (CDS needed).',
        'job_id': job_id,
    })


def run_uniprot_analysis(
    uniprot_id: str,
    *,
    editor: str = 'BOTH',
    pam_type: str = 'NGG',
    window_min: int = 4,
    window_max: int = 8,
    top_sgrnas: int = 5,
    organism: str = 'eukarya',
    mutation: str | None = None,
) -> dict[str, Any]:
    resolved = resolve_uniprot_inputs(uniprot_id)
    seq = resolved['aa_sequence']
    cds = resolved['cds']
    transcript_id = resolved['transcript_id']
    mut = parse_mutation(mutation)
    job_id = f'{resolved["accession"]}_{uuid.uuid4().hex[:10]}'

    wt_preds = run_signalp({'WT': seq}, organism=organism, job_id=f'{job_id}_wt', want_plots=True)
    wt_pred = wt_preds.get('WT')
    if not wt_pred:
        raise SaffronUserError('SignalP did not return a WT prediction.')

    if wt_pred.prediction == 'OTHER' or not wt_pred.sp_end:
        return _with_backend({
            'mode': 'uniprot',
            **{k: resolved[k] for k in ('accession', 'gene_symbol', 'transcript_id', 'protein_length')},
            'organism': organism,
            'wt_signalp': prediction_to_dict(wt_pred),
            'guide_rows': [],
            'guides_available': True,
            'no_sp': True,
            'message': 'SignalP did not predict a signal peptide for this protein. Guide design in an SP region is not available.',
            'focus_mutation': None,
            'job_id': job_id,
            'form_meta': {
                'editor': editor,
                'pam_type': pam_type,
                'window_min': window_min,
                'window_max': window_max,
                'top_sgrnas': top_sgrnas,
                'mutation': mutation or '',
            },
        })

    sp_start = wt_pred.sp_start or 1
    sp_end = wt_pred.sp_end

    guide_rows = design_sp_guides(
        cds,
        transcript_id,
        sp_start,
        sp_end,
        editor=editor,
        window_min=window_min,
        window_max=window_max,
        pam_type=pam_type,
        top_sgrnas=top_sgrnas,
    )

    # Unique mutant proteins
    mutant_seqs: dict[str, str] = {}
    for row in guide_rows:
        mid = f"m{row['position']}_{row['wt_aa']}{row['position']}{row['mut_aa']}"
        row['mutant_id'] = mid
        row['highlighted'] = bool(
            mut and mut[1] == row['position'] and mut[0] == row['wt_aa'] and mut[2] == row['mut_aa']
        )
        if mid not in mutant_seqs:
            mutant_seqs[mid] = _mutant_protein(seq, row['position'], row['mut_aa'])

    # Optional focus mutation not reachable by guides: still score it
    focus_extra = None
    if mut:
        wt_aa, pos, mut_aa = mut
        apply_mutation(seq, wt_aa, pos, mut_aa)  # validate
        mid = f'focus_{wt_aa}{pos}{mut_aa}'
        if mid not in mutant_seqs and f'm{pos}_{wt_aa}{pos}{mut_aa}' not in mutant_seqs:
            mutant_seqs[mid] = _mutant_protein(seq, pos, mut_aa)
            focus_extra = mid

    mut_preds: dict[str, SignalPPrediction] = {}
    if mutant_seqs:
        mut_preds = run_signalp(mutant_seqs, organism=organism, job_id=f'{job_id}_mut', want_plots=True)

    enriched = _attach_signalp_to_rows(guide_rows, seq, wt_pred, mut_preds)
    annotate_guide_rows(enriched, accession=resolved['accession'])

    focus_row = None
    if mut:
        wt_aa, pos, mut_aa = mut
        label = f'{wt_aa}{pos}{mut_aa}'
        # Prefer matching guide row
        matches = [r for r in enriched if r['position'] == pos and r['wt_aa'] == wt_aa and r['mut_aa'] == mut_aa]
        if matches:
            focus_row = {**matches[0], 'label': label}
        elif focus_extra and focus_extra in mut_preds:
            deltas = compute_deltas(wt_pred, mut_preds[focus_extra])
            focus_row = {
                'position': pos,
                'wt_aa': wt_aa,
                'mut_aa': mut_aa,
                'label': label,
                'sp_prediction': mut_preds[focus_extra].prediction,
                'sp_prob': mut_preds[focus_extra].sp_prob,
                'cs_before': mut_preds[focus_extra].cs_before,
                'cs_prob': mut_preds[focus_extra].cs_prob,
                'plot_path': mut_preds[focus_extra].plot_path,
                'sgrna_seq': None,
                **deltas,
            }
            annotate_guide_rows([focus_row], accession=resolved['accession'])

    # Sort: paper-pathogenic first, then highlighted, then |delta|, then position
    def sort_key(r):
        d = r.get('delta_wt_class_prob')
        return (
            0 if r.get('paper_pathogenic') else 1,
            0 if r.get('highlighted') else 1,
            0 if d is None else -abs(d),
            r['position'],
            r.get('editor_used') or '',
        )

    enriched.sort(key=sort_key)

    return _with_backend({
        'mode': 'uniprot',
        'accession': resolved['accession'],
        'gene_symbol': resolved['gene_symbol'],
        'transcript_id': transcript_id,
        'protein_length': resolved['protein_length'],
        'aa_sequence': seq,
        'organism': organism,
        'wt_signalp': prediction_to_dict(wt_pred),
        'guide_rows': enriched,
        'guides_available': True,
        'no_sp': False,
        'message': '',
        'focus_mutation': focus_row,
        'job_id': job_id,
        'sp_span': {'start': sp_start, 'end': sp_end},
        'form_meta': {
            'editor': editor,
            'pam_type': pam_type,
            'window_min': window_min,
            'window_max': window_max,
            'top_sgrnas': top_sgrnas,
            'mutation': mutation or '',
        },
    })
