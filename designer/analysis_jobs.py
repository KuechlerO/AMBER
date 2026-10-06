"""Background AMBER analysis jobs (avoids gateway timeouts; defers CADD)."""

from __future__ import annotations

import threading
from typing import Any

from django.core.cache import cache

from . import result_store
from .pipeline import enrich_guide_rows_with_cadd, normalize_duplicate_mode
from .services import UserInputError, normalize_input_id, run_analysis

JOB_PREFIX = 'amber_job:'
JOB_TIMEOUT = 60 * 60  # 1 hour
CADD_PREFIX = 'amber_cadd:'

_lock = threading.Lock()


def _job_cache_key(session_key: str) -> str:
    return f'{JOB_PREFIX}{session_key}'


def _cadd_cache_key(session_key: str) -> str:
    return f'{CADD_PREFIX}{session_key}'


def _form_adapter(form: dict[str, Any]):
    """Minimal QueryDict-like wrapper so view helpers accept plain dicts."""

    class _Adapter:
        def get(self, key, default=None):
            return form.get(key, default)

        def getlist(self, key):
            value = form.get(key, [])
            if value is None:
                return []
            if isinstance(value, (list, tuple)):
                return list(value)
            return [value]

    return _Adapter()


def _build_form_data(form: dict[str, Any], analysis_result: dict[str, Any]) -> dict[str, Any]:
    from .views import build_form_data, resolve_score_filter, resolve_selected_columns

    adapter = _form_adapter(form)
    editor = form.get('editor') or 'BOTH'
    alpha_threshold = form.get('alpha_threshold', '0.5')
    score_filter_mode, display_score_cutoff = resolve_score_filter(
        adapter, default_cutoff=alpha_threshold
    )
    selected_columns = form.get('selected_columns')
    if not selected_columns:
        selected_columns = resolve_selected_columns(adapter, editor)
    uniprot_id = form.get('uniprot_id') or ''
    return build_form_data(
        uniprot_id=uniprot_id,
        editor=editor,
        editor_display='ABE & CBE' if editor == 'BOTH' else editor,
        alpha_threshold=alpha_threshold,
        top_sgrnas=form.get('top_sgrnas', '5'),
        gene_name=analysis_result.get('gene_symbol') or uniprot_id,
        gene_symbol=analysis_result.get('gene_symbol') or '',
        uniprot_accession=analysis_result.get('uniprot_accession') or uniprot_id,
        protein_length=analysis_result.get('protein_length') or 0,
        transcript_id=analysis_result.get('transcript_id') or '',
        selected_columns=selected_columns,
        score_filter_mode=score_filter_mode,
        display_score_cutoff=display_score_cutoff,
        filter_threshold='true' if score_filter_mode == 'above' else 'false',
        window_min=form.get('window_min', '4'),
        window_max=form.get('window_max', '8'),
        duplicate_mode=normalize_duplicate_mode(form.get('duplicate_mode', 'best')),
        pam_type=form.get('pam_type', 'NGG'),
        cadd_warning=analysis_result.get('cadd_warning'),
        cadd_status=analysis_result.get('cadd_status', 'pending'),
    )


def _display_rows(full_rows, form_data):
    from .views import _rows_for_display

    return _rows_for_display(full_rows, form_data)


def _start_cadd_enrichment(session_key: str) -> None:
    key = _cadd_cache_key(session_key)
    with _lock:
        state = cache.get(key)
        if state and state.get('status') == 'running':
            return
        cache.set(key, {'status': 'running', 'error': None}, JOB_TIMEOUT)

    def _work() -> None:
        try:
            data = result_store.load_analysis_results_by_key(session_key)
            rows = list(data.get('full_rows') or data.get('raw_rows') or [])
            form_data = dict(data.get('form_data') or {})
            if not rows:
                cache.set(key, {'status': 'done', 'error': None}, JOB_TIMEOUT)
                return
            enriched, cadd_warning = enrich_guide_rows_with_cadd(rows)
            form_data['cadd_warning'] = cadd_warning
            form_data['cadd_status'] = 'done'
            filtered = _display_rows(enriched, form_data)
            result_store.save_analysis_results_by_key(
                session_key,
                full_rows=enriched,
                form_data=form_data,
                filtered_rows=filtered,
                no_guide_positions=data.get('no_guide_positions') or [],
                cadd_status='done',
            )
            cache.set(key, {'status': 'done', 'error': None}, JOB_TIMEOUT)
        except Exception as exc:  # noqa: BLE001
            data = result_store.load_analysis_results_by_key(session_key)
            form_data = dict(data.get('form_data') or {})
            form_data['cadd_warning'] = f'CADD enrichment failed: {exc}'
            form_data['cadd_status'] = 'error'
            if data:
                result_store.save_analysis_results_by_key(
                    session_key,
                    full_rows=data.get('full_rows') or [],
                    form_data=form_data,
                    filtered_rows=data.get('filtered_rows'),
                    no_guide_positions=data.get('no_guide_positions') or [],
                    cadd_status='error',
                )
            cache.set(key, {'status': 'error', 'error': str(exc)}, JOB_TIMEOUT)

    thread = threading.Thread(target=_work, daemon=True, name=f'amber-cadd-{session_key[:8]}')
    thread.start()


def start_analysis(session_key: str, form: dict[str, Any]) -> bool:
    """Start analysis in a daemon thread. Returns True if a new job was started."""
    key = _job_cache_key(session_key)
    with _lock:
        state = cache.get(key)
        if state and state.get('status') == 'running':
            return False
        cache.set(key, {'status': 'running', 'error': None}, JOB_TIMEOUT)
        cache.delete(_cadd_cache_key(session_key))

    def _work() -> None:
        try:
            from .views import _normalize_analysis_result, parse_editing_window

            form_local = dict(form)
            form_local['uniprot_id'] = normalize_input_id(form_local.get('uniprot_id', ''))
            window_min, window_max = parse_editing_window(_form_adapter(form_local))
            analysis_result = _normalize_analysis_result(run_analysis(
                uniprot_id=form_local.get('uniprot_id', ''),
                editor=form_local.get('editor', 'BOTH'),
                alpha_threshold=form_local.get('alpha_threshold', '0.5'),
                top_sgrnas=form_local.get('top_sgrnas', '5'),
                window_min=window_min,
                window_max=window_max,
                pam_type=form_local.get('pam_type', 'NGG'),
            ))
            full_rows = analysis_result['guide_rows']
            no_guide_full = analysis_result.get('no_guide_positions') or []
            form_data = _build_form_data(form_local, analysis_result)
            filtered = _display_rows(full_rows, form_data)
            cadd_status = analysis_result.get('cadd_status', 'pending')
            form_data['cadd_status'] = cadd_status
            result_store.save_analysis_results_by_key(
                session_key,
                full_rows=full_rows,
                form_data=form_data,
                filtered_rows=filtered,
                no_guide_positions=no_guide_full,
                cadd_status=cadd_status,
            )
            cache.set(key, {'status': 'done', 'error': None}, JOB_TIMEOUT)
            if cadd_status == 'pending':
                _start_cadd_enrichment(session_key)
        except UserInputError as exc:
            cache.set(key, {'status': 'error', 'error': str(exc)}, JOB_TIMEOUT)
        except Exception as exc:  # noqa: BLE001
            cache.set(key, {'status': 'error', 'error': f'Analysis failed: {exc}'}, JOB_TIMEOUT)

    thread = threading.Thread(target=_work, daemon=True, name=f'amber-{session_key[:8]}')
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


def get_cadd_job_status(session_key: str) -> dict[str, Any]:
    """Return deferred CADD job status."""
    data = result_store.load_analysis_results_by_key(session_key)
    stored = data.get('cadd_status') or data.get('form_data', {}).get('cadd_status')
    if stored in ('done', 'error', 'skipped'):
        warning = None
        if stored == 'error':
            warning = data.get('form_data', {}).get('cadd_warning')
        return {'status': stored, 'error': warning}
    state = cache.get(_cadd_cache_key(session_key))
    if state:
        return dict(state)
    if stored == 'pending':
        return {'status': 'pending', 'error': None}
    return {'status': 'idle', 'error': None}


def clear_job_state(session_key: str) -> None:
    """Test helper."""
    cache.delete(_job_cache_key(session_key))
    cache.delete(_cadd_cache_key(session_key))
