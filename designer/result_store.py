"""File-backed cache storage for analysis results (not session DB)."""

from django.core.cache import cache

CACHE_PREFIX = 'crispr_analysis:'
CACHE_TIMEOUT = 60 * 60 * 24  # 24 hours


def _session_cache_key(request) -> str:
    if not request.session.session_key:
        request.session.save()
    return f'{CACHE_PREFIX}{request.session.session_key}'


def save_analysis_results(
    request,
    *,
    full_rows,
    form_data,
    filtered_rows=None,
    no_guide_positions=None,
) -> None:
    cache.set(
        _session_cache_key(request),
        {
            'full_rows': full_rows,
            'no_guide_positions': no_guide_positions or [],
            'filtered_rows': filtered_rows,
            'raw_rows': full_rows,  # legacy alias
            'form_data': form_data,
        },
        CACHE_TIMEOUT,
    )


def load_analysis_results(request) -> dict:
    return cache.get(_session_cache_key(request)) or {}


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
