"""Input validation and entry points for SAFFRON analyses."""

from __future__ import annotations

import requests

from designer.services import (
    UserInputError,
    is_valid_ensembl_id,
    is_valid_uniprot_id,
    normalize_input_id,
)

from .pipeline import (
    SaffronUserError,
    run_sequence_analysis,
    run_uniprot_analysis,
)


def run_saffron_analysis(form: dict) -> dict:
    organism = (form.get('organism') or 'eukarya').strip().lower()
    if organism not in ('eukarya', 'other'):
        raise UserInputError('Organism must be eukarya or other.')

    input_mode = (form.get('input_mode') or 'uniprot').strip().lower()
    mutation = (form.get('mutation') or '').strip() or None

    if input_mode == 'sequence':
        seq = form.get('aa_sequence') or ''
        try:
            return run_sequence_analysis(seq, organism=organism, mutation=mutation)
        except SaffronUserError as exc:
            raise UserInputError(str(exc)) from exc

    # UniProt / Ensembl path
    uid = normalize_input_id(form.get('uniprot_id') or '')
    if not uid:
        raise UserInputError('Please provide a UniProt or Ensembl ID.')
    if not (is_valid_uniprot_id(uid) or is_valid_ensembl_id(uid)):
        raise UserInputError('Invalid UniProt ID / Ensembl ID provided.')

    editor = (form.get('editor') or 'BOTH').strip().upper()
    if editor not in ('ABE', 'CBE', 'BOTH'):
        raise UserInputError('Invalid editor selected.')

    pam_type = (form.get('pam_type') or 'NGG').strip().upper()
    if pam_type not in ('NGG', 'NG', 'NAG'):
        raise UserInputError('Invalid PAM type.')

    try:
        window_min = int(form.get('window_min') or 4)
        window_max = int(form.get('window_max') or 8)
    except (TypeError, ValueError) as exc:
        raise UserInputError('Invalid editing window values.') from exc
    if not (1 <= window_min <= 20 and 1 <= window_max <= 20) or window_min > window_max:
        raise UserInputError('Editing window must be between 1 and 20 with min ≤ max.')

    try:
        top_sgrnas = int(form.get('top_sgrnas') or 5)
    except (TypeError, ValueError) as exc:
        raise UserInputError('Invalid number of top guides.') from exc

    try:
        return run_uniprot_analysis(
            uid,
            editor=editor,
            pam_type=pam_type,
            window_min=window_min,
            window_max=window_max,
            top_sgrnas=top_sgrnas,
            organism=organism,
            mutation=mutation,
        )
    except SaffronUserError as exc:
        raise UserInputError(str(exc)) from exc
    except requests.RequestException as exc:
        raise UserInputError(f'Could not retrieve sequence data: {exc}') from exc
    except FileNotFoundError as exc:
        raise UserInputError(str(exc)) from exc
    except ValueError as exc:
        raise UserInputError(str(exc)) from exc
