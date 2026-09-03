"""Background SAFFRON analysis jobs (avoids gateway timeout on long SignalP runs)."""

from __future__ import annotations

import threading
from typing import Any

from django.core.cache import cache

from designer.services import UserInputError

from . import result_store
from .services import run_saffron_analysis

JOB_PREFIX = 'saffron_job:'
JOB_TIMEOUT = 60 * 60  # 1 hour

_lock = threading.Lock()


def _job_cache_key(session_key: str) -> str:
    return f'{JOB_PREFIX}{session_key}'


def start_analysis(session_key: str, form: dict[str, Any]) -> bool:
    """Start analysis in a daemon thread. Returns True if a new job was started."""
    key = _job_cache_key(session_key)
    with _lock:
        state = cache.get(key)
        if state and state.get('status') == 'running':
            return False
        cache.set(key, {'status': 'running', 'error': None}, JOB_TIMEOUT)

    def _work() -> None:
        try:
            payload = run_saffron_analysis(form)
            payload['form_data'] = form
            result_store.save_analysis_results_by_key(session_key, payload=payload)
            cache.set(key, {'status': 'done', 'error': None}, JOB_TIMEOUT)
        except UserInputError as exc:
            cache.set(key, {'status': 'error', 'error': str(exc)}, JOB_TIMEOUT)
        except Exception as exc:  # noqa: BLE001
            cache.set(key, {'status': 'error', 'error': f'Analysis failed: {exc}'}, JOB_TIMEOUT)

    thread = threading.Thread(target=_work, daemon=True, name=f'saffron-{session_key[:8]}')
    thread.start()
    return True


def get_status(session_key: str) -> dict[str, Any]:
    """Return job status for polling: idle | running | done | error."""
    key = _job_cache_key(session_key)
    state = cache.get(key)
    if state:
        return dict(state)
    if result_store.load_analysis_results_by_key(session_key):
        return {'status': 'done', 'error': None}
    return {'status': 'idle', 'error': None}


def clear_job_state(session_key: str) -> None:
    """Test helper."""
    cache.delete(_job_cache_key(session_key))
