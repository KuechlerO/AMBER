"""Session-scoped cache for SAFFRON analysis results."""

from django.core.cache import cache

CACHE_PREFIX = 'saffron_analysis:'
CACHE_TIMEOUT = 60 * 60 * 24  # 24 hours


def _session_cache_key(request) -> str:
    if not request.session.session_key:
        request.session.save()
    return f'{CACHE_PREFIX}{request.session.session_key}'


def save_analysis_results(request, *, payload: dict) -> None:
    cache.set(_session_cache_key(request), payload, CACHE_TIMEOUT)


def load_analysis_results(request) -> dict:
    return cache.get(_session_cache_key(request)) or {}
