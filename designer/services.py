# services.py
import requests
from .pipeline import run_pipeline


# For error analysis, to distinguish between user and server errors:
class UserInputError(Exception):
    """Error due to invalid user input"""
    pass


def normalize_input_id(raw_id: str) -> str:
    """Strip surrounding whitespace and remove internal spaces from IDs."""
    return ''.join((raw_id or '').split())


# tests if UniProt ID is valid
def is_valid_uniprot_id(uniprot_id: str) -> bool:
    uid = normalize_input_id(uniprot_id)
    if not uid or uid[0] not in ['P', 'O', 'Q']:
        return False

    url = f"https://rest.uniprot.org/uniprotkb/{uid}.json"
    response = requests.get(url)
    return response.status_code == 200

# tests if Ensembl ID is valid
def is_valid_ensembl_id(ensembl_id: str) -> bool:
    eid = normalize_input_id(ensembl_id)
    if not eid.startswith('ENST'):
        return False

    url = f"https://rest.ensembl.org/sequence/id/{eid}?type=cds"

    try:
        response = requests.get(
            url,
            headers={'Content-Type': 'text/plain', 'Accept': 'text/plain'},
            timeout=60,
        )
        return response.status_code == 200
    except requests.RequestException:
        return False

# analysis pipeline

def run_analysis(uniprot_id, editor, alpha_threshold, top_sgrnas, window_min, window_max, duplicate_mode='hide', pam_type='NGG'):
    uniprot_id = normalize_input_id(uniprot_id)
    print('RUN ANALYSIS')
    print('UniProt ID:', uniprot_id)
    print('Editor:', editor)
    print('Alpha threshold:', alpha_threshold)
    print('Top sgRNAs per position:', top_sgrnas)
    print('Editing window positions:', window_min, '-', window_max)
    print('Duplicate mode:', duplicate_mode)

    try:
        window_min = int(window_min)
        window_max = int(window_max)
    except (TypeError, ValueError):
        raise UserInputError('Invalid editing window values.')

    if not (1 <= window_min <= 20 and 1 <= window_max <= 20):
        print('Editing window must be between 1 and 20.')
        raise UserInputError('Editing window must be between 1 and 20.')

    if window_min > window_max:
        print('window_min cannot be greater than window_max.')
        raise UserInputError('window_min cannot be greater than window_max.')

    if not uniprot_id or len(uniprot_id) < 5:
        print('Invalid UniProt ID/ Ensembl ID input.')
        raise UserInputError('Invalid UniProt ID/ Ensembl ID provided.')

    if (is_valid_uniprot_id(uniprot_id) == False and is_valid_ensembl_id(uniprot_id) == False):
        print('Invalid UniProt ID/ Ensembl ID provided.')
        raise UserInputError('Invalid UniProt ID/ Ensembl provided.')

    if editor not in ['ABE', 'CBE', 'BOTH']:
        print('Invalid editor selected.')
        raise UserInputError('Invalid editor selected.')

    try:
        float(alpha_threshold)
    except (TypeError, ValueError):
        print('Invalid AlphaMissense threshold.')
        raise UserInputError('Invalid AlphaMissense threshold.')

    try:
        int(top_sgrnas)
    except (TypeError, ValueError):
        print('Invalid number of top sgRNAs.')
        raise UserInputError('Invalid number of top sgRNAs.')

    try:
        return run_pipeline(
            uniprot_id=uniprot_id,
            editor=editor,
            alpha_threshold=alpha_threshold,
            top_sgrnas=top_sgrnas,
            window_min=window_min,
            window_max=window_max,
            pam_type=pam_type,
        )

    except FileNotFoundError as e:
        raise UserInputError(str(e))
    except ValueError as e:
        raise UserInputError(str(e))
    except requests.RequestException as e:
        raise UserInputError(f"Could not retrieve sequence data: {e}")
