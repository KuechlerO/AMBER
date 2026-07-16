"""Guide coverage statistics and stacked bar figure (editor × score)."""

from __future__ import annotations

from collections import defaultdict
from typing import Any

# Segment key -> (legend label, color) — muted palette that fits the AMBER UI
SEGMENT_STYLE = {
    ('ABE', 'high'): ('ABE ≥ cutoff', '#2f9e6b'),
    ('ABE', 'low'): ('ABE < cutoff', '#a8d5c0'),
    ('CBE', 'high'): ('CBE ≥ cutoff', '#3b7ea8'),
    ('CBE', 'low'): ('CBE < cutoff', '#a9c7db'),
    ('both', 'high'): ('ABE & CBE ≥ cutoff', '#7c5cbf'),
    ('both', 'low'): ('ABE & CBE < cutoff', '#c4b5e0'),
    ('uncovered', None): ('Not covered', '#d7d2c8'),
}

SEGMENT_ORDER = [
    ('ABE', 'high'),
    ('ABE', 'low'),
    ('CBE', 'high'),
    ('CBE', 'low'),
    ('both', 'high'),
    ('both', 'low'),
    ('uncovered', None),
]


def _editor_class(editors: set[str]) -> str:
    has_abe = 'ABE' in editors
    has_cbe = 'CBE' in editors
    if has_abe and has_cbe:
        return 'both'
    if has_abe:
        return 'ABE'
    if has_cbe:
        return 'CBE'
    return 'ABE'


def classify_positions(
    rows: list[dict[str, Any]],
    cutoff: float,
) -> dict[int, tuple[str, str]]:
    """
    Per AA position: (editor_class, score_bin) where score_bin is 'high' if max
    alpha_score at that position is >= cutoff, else 'low'.
    """
    by_pos: dict[int, dict[str, Any]] = defaultdict(lambda: {'editors': set(), 'scores': []})
    for row in rows:
        pos = row.get('position')
        if pos is None:
            continue
        try:
            pos_i = int(pos)
        except (TypeError, ValueError):
            continue
        editor = (row.get('editor_used') or '').strip().upper()
        if editor in ('ABE', 'CBE'):
            by_pos[pos_i]['editors'].add(editor)
        score = row.get('alpha_score')
        try:
            if score is not None:
                by_pos[pos_i]['scores'].append(float(score))
        except (TypeError, ValueError):
            pass

    classified: dict[int, tuple[str, str]] = {}
    for pos, info in by_pos.items():
        editors = info['editors'] or {'ABE'}
        best = max(info['scores']) if info['scores'] else 0.0
        score_bin = 'high' if best >= cutoff else 'low'
        classified[pos] = (_editor_class(editors), score_bin)
    return classified


def compute_coverage_stats(
    rows: list[dict[str, Any]],
    protein_length: int,
    cutoff: float,
    *,
    guide_count: int | None = None,
    no_guide_count: int = 0,
) -> dict[str, Any]:
    """Return summary numbers and segment counts/percentages for the coverage bar."""
    try:
        pro_len = int(protein_length)
    except (TypeError, ValueError):
        pro_len = 0
    pro_len = max(pro_len, 0)

    classified = classify_positions(rows, float(cutoff))
    editable = len(classified)
    uncovered = max(pro_len - editable, 0)

    segment_counts: dict[tuple, int] = {key: 0 for key in SEGMENT_ORDER}
    for editor_cls, score_bin in classified.values():
        key = (editor_cls, score_bin)
        if key in segment_counts:
            segment_counts[key] += 1
    segment_counts[('uncovered', None)] = uncovered

    segment_pct: dict[tuple, float] = {}
    for key, count in segment_counts.items():
        segment_pct[key] = (100.0 * count / pro_len) if pro_len else 0.0

    below_cutoff = sum(
        1 for _pos, (_ed, score_bin) in classified.items() if score_bin == 'low'
    )
    above_cutoff = sum(
        1 for _pos, (_ed, score_bin) in classified.items() if score_bin == 'high'
    )

    coverage_pct = (100.0 * editable / pro_len) if pro_len else 0.0
    return {
        'protein_length': pro_len,
        'editable_positions': editable,
        'uncovered_positions': uncovered,
        'coverage_pct': round(coverage_pct, 1),
        'above_cutoff_positions': above_cutoff,
        'below_cutoff_positions': below_cutoff,
        'guide_count': guide_count if guide_count is not None else len(rows),
        'no_guide_count': no_guide_count,
        'cutoff': float(cutoff),
        'segment_counts': segment_counts,
        'segment_pct': segment_pct,
        'classified': classified,
    }


def build_coverage_figure(stats: dict[str, Any]):
    """Horizontal stacked bar as % of protein length (0–100)."""
    import plotly.graph_objects as go

    pro_len = max(int(stats.get('protein_length') or 0), 1)
    counts = stats.get('segment_counts') or {}
    pcts = stats.get('segment_pct') or {}

    fig = go.Figure()
    for key in SEGMENT_ORDER:
        count = int(counts.get(key, 0))
        if count <= 0:
            continue
        label, color = SEGMENT_STYLE[key]
        pct = float(pcts.get(key, 100.0 * count / pro_len))
        fig.add_trace(
            go.Bar(
                y=[''],
                x=[pct],
                name=label,
                orientation='h',
                marker=dict(color=color, line=dict(color='rgba(255,255,255,0.85)', width=1)),
                hovertemplate=(
                    f'{label}<br>{count} AA · {pct:.1f}% of protein<extra></extra>'
                ),
            )
        )

    fig.update_layout(
        barmode='stack',
        height=160,
        margin=dict(l=8, r=8, t=8, b=56),
        showlegend=True,
        legend=dict(
            orientation='h',
            yanchor='top',
            y=-0.45,
            x=0,
            font=dict(size=11, color='#4b5563'),
            bgcolor='rgba(0,0,0,0)',
        ),
        xaxis=dict(
            title=dict(text='Fragment of covered positions (%)', font=dict(size=12, color='#4b5563')),
            range=[0, 100],
            ticksuffix='%',
            fixedrange=True,
            gridcolor='rgba(0,0,0,0.06)',
            zeroline=False,
        ),
        yaxis=dict(visible=False, fixedrange=True),
        plot_bgcolor='rgba(0,0,0,0)',
        paper_bgcolor='rgba(0,0,0,0)',
        font=dict(family='inherit', color='#374151'),
    )
    return fig
