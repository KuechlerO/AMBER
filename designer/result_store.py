"""File-backed cache storage for analysis results (not session DB)."""

from django.core.cache import cache

CACHE_PREFIX = 'crispr_analysis:'
CACHE_TIMEOUT = 60 * 60 * 24  # 24 hours


def session_key_from_request(request) -> str:
    if not request.session.session_key:
        request.session.save()
    return request.session.session_key


def analysis_cache_key(session_key: str) -> str:
    return f'{CACHE_PREFIX}{session_key}'


def _session_cache_key(request) -> str:
    return analysis_cache_key(session_key_from_request(request))


def save_analysis_results_by_key(
    session_key: str,
    *,
    full_rows,
    form_data,
    filtered_rows=None,
    no_guide_positions=None,
    cadd_status=None,
) -> None:
    existing = load_analysis_results_by_key(session_key)
    payload = {
        'full_rows': full_rows,
        'no_guide_positions': no_guide_positions or [],
        'filtered_rows': filtered_rows,
        'raw_rows': full_rows,  # legacy alias
        'form_data': form_data,
        'cadd_status': (
            cadd_status
            if cadd_status is not None
            else existing.get('cadd_status')
            or form_data.get('cadd_status')
            or 'done'
        ),
    }
    cache.set(analysis_cache_key(session_key), payload, CACHE_TIMEOUT)


def load_analysis_results_by_key(session_key: str) -> dict:
    return cache.get(analysis_cache_key(session_key)) or {}


def save_analysis_results(
    request,
    *,
    full_rows,
    form_data,
    filtered_rows=None,
    no_guide_positions=None,
    cadd_status=None,
) -> None:
    save_analysis_results_by_key(
        session_key_from_request(request),
        full_rows=full_rows,
        form_data=form_data,
        filtered_rows=filtered_rows,
        no_guide_positions=no_guide_positions,
        cadd_status=cadd_status,
    )


def load_analysis_results(request) -> dict:
    return load_analysis_results_by_key(session_key_from_request(request))


def get_full_rows(request):
    data = load_analysis_results(request)
    return data.get('full_rows') or data.get('raw_rows', [])


def get_filtered_rows(request):
    """Rows currently shown on the results table (score filter + duplicate mode)."""
    data = load_analysis_results(request)
    return data.get('filtered_rows') or data.get('full_rows') or data.get('raw_rows', [])


def get_raw_rows(request):
    return get_full_rows(request)


def get_form_data(request):
    data = load_analysis_results(request)
    return data.get('form_data', {})


def get_no_guide_positions(request):
    data = load_analysis_results(request)
    return data.get('no_guide_positions', [])


def get_cadd_status(request) -> str:
    data = load_analysis_results(request)
    return data.get('cadd_status') or data.get('form_data', {}).get('cadd_status') or 'done'
