"""Session-scoped cache for SAFFRON analysis results."""

from django.core.cache import cache

CACHE_PREFIX = 'saffron_analysis:'
CACHE_TIMEOUT = 60 * 60 * 24  # 24 hours


def session_key_from_request(request) -> str:
    if not request.session.session_key:
        request.session.save()
    return request.session.session_key


def analysis_cache_key(session_key: str) -> str:
    return f'{CACHE_PREFIX}{session_key}'


def save_analysis_results_by_key(session_key: str, *, payload: dict) -> None:
    cache.set(analysis_cache_key(session_key), payload, CACHE_TIMEOUT)


def load_analysis_results_by_key(session_key: str) -> dict:
    return cache.get(analysis_cache_key(session_key)) or {}


def _session_cache_key(request) -> str:
    return analysis_cache_key(session_key_from_request(request))


def save_analysis_results(request, *, payload: dict) -> None:
    save_analysis_results_by_key(session_key_from_request(request), payload=payload)


def load_analysis_results(request) -> dict:
    return load_analysis_results_by_key(session_key_from_request(request))
