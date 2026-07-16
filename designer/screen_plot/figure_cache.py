"""In-memory cache for full-screen Plotly figures (avoids nginx timeout on build)."""

from __future__ import annotations

import json
import threading
from typing import Any

from .plotly_builder import build_full_screen_figure
from .data_loader import ScreenDataStore

_cache_lock = threading.Lock()
_figure_cache: dict[tuple, dict[str, Any]] = {}
_building: set[tuple] = set()
_build_errors: dict[tuple, str] = {}


def cache_key(gene: str, guide_positions: list[int]) -> tuple:
    return gene, tuple(guide_positions)


def get_cached_figure(key: tuple) -> dict[str, Any] | None:
    with _cache_lock:
        return _figure_cache.get(key)


def is_building(key: tuple) -> bool:
    with _cache_lock:
        return key in _building


def get_build_error(key: tuple) -> str | None:
    with _cache_lock:
        return _build_errors.get(key)


def start_full_figure_build(
    gene: str,
    guide_positions: list[int],
    store: ScreenDataStore | None = None,
) -> bool:
    """Start background build. Returns True if a new build was started."""
    key = cache_key(gene, guide_positions)
    with _cache_lock:
        if key in _figure_cache or key in _building:
            return False
        _building.add(key)
        _build_errors.pop(key, None)

    def _work():
        try:
            store_local = store or ScreenDataStore.get()
            fig = build_full_screen_figure(gene, guide_positions=guide_positions, store=store_local)
            payload = json.loads(fig.to_json())
            with _cache_lock:
                _figure_cache[key] = payload
        except Exception as exc:
            with _cache_lock:
                _build_errors[key] = str(exc)
        finally:
            with _cache_lock:
                _building.discard(key)

    threading.Thread(target=_work, daemon=True).start()
    return True


def reset_for_tests():
    with _cache_lock:
        _figure_cache.clear()
        _building.clear()
        _build_errors.clear()
