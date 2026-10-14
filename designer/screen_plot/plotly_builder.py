"""Plotly screen enrichment figures (overview + multi-marker LFC panels)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from .data_loader import ScreenDataStore
from .gene_lookup import SCREEN_MARKERS
from .helpers import ext_loc_splice, is_KO

BE_TYPES = ['ABE', 'CBE']

MARKER_COLORS = {
    'TNFa': '#22c55e',
    'IFNG': '#f97316',
    'PD1': '#3b82f6',
    'CD25': '#a855f7',
}
MULTI_MARKER_COLOR = '#eab308'
_DEFAULT_GUIDE_COLOR = 'rgb(30,60,120)'
_DOMAIN_NEUTRAL = 'rgb(190,190,190)'
_CLINVAR_MAX_LANES = 8
_CLINVAR_COLORS = {
    'Pathogenic': '#b91c1c',
    'Likely pathogenic': '#ea580c',
    'Uncertain significance': '#2563eb',
    'Likely benign': '#65a30d',
    'Benign': '#16a34a',
    'Conflicting classifications': '#7c3aed',
    'Unknown': '#9ca3af',
}


def _diverging_norm(lfc_series: pd.Series):
    max_frame = lfc_series.replace([np.inf, -np.inf], np.nan).dropna()
    if max_frame.empty:
        y_min, y_max = -0.75, 0.75
    else:
        y_max = float(max_frame.nlargest(2).mean())
        y_min = float(max_frame.nsmallest(2).mean())
        if y_min >= -0.75:
            y_min = -0.75
        if y_max <= 0.75:
            y_max = 0.75
    return y_min * 0.66, 0.0, y_max * 0.66


def _lfc_to_color(value: float, vmin: float, vcenter: float, vmax: float) -> str:
    if value <= vcenter:
        t = 0.5 * (value - vmin) / (vcenter - vmin) if vcenter != vmin else 0.5
    else:
        t = 0.5 + 0.5 * (value - vcenter) / (vmax - vcenter) if vmax != vcenter else 0.5
    t = max(0.0, min(1.0, t))
    if t <= 0.5:
        t2 = t * 2
        r = int(69 + (245 - 69) * t2)
        g = int(117 + (245 - 117) * t2)
        b = int(180 + (245 - 180) * t2)
    else:
        t2 = (t - 0.5) * 2
        r = int(245 + (215 - 245) * t2)
        g = int(245 + (48 - 245) * t2)
        b = int(245 + (39 - 245) * t2)
    return f'rgb({r},{g},{b})'


def _collapse_duplicate_index(df: pd.DataFrame) -> pd.DataFrame:
    if not df.index.duplicated().any():
        return df
    agg = {c: 'mean' for c in ['LFC', 'control_mean', 'treat_mean', 'p.twosided'] if c in df.columns}
    for col in df.columns:
        if col not in agg:
            agg[col] = 'first'
    return df.groupby(df.index).agg(agg)


def _gene_marker_slice(store: ScreenDataStore, gene: str, marker: str, be_type: str) -> pd.DataFrame:
    pro_len = store.protein_length[gene]
    df = store.data_grouped.loc[gene, marker, be_type]
    df = _collapse_duplicate_index(df)
    df = df.reindex(np.arange(0, pro_len))
    df[['control_mean', 'LFC']] = df[['control_mean', 'LFC']].fillna(0)
    return df


def _ko_reference_lfc(store: ScreenDataStore, gene: str, marker: str, be_type: str) -> float | None:
    pro_len = store.protein_length.get(gene, 1)
    all_data = store.data_cutoff
    if all_data is None:
        return None
    ko = all_data[
        (all_data['marker'] == marker)
        & (all_data['be_type'] == be_type)
        & (all_data['gene'] == gene)
    ].copy()
    ko = ko[ko['sites_mapped'].apply(lambda x: any(is_KO(a) for a in x))]
    if ko.empty:
        return None
    ko['loc'] = ko['sites_mapped'].apply(
        lambda x: min(ext_loc_splice(a) for a in x if a != 'None')
    )
    ko = ko[ko['loc'] / pro_len <= 0.5]
    if ko.empty:
        return None
    ko_lfc = ko.groupby('sgrna')['LFC'].mean().sort_values(ascending=False, key=abs).iloc[:3]
    if ko_lfc.empty:
        return None
    return float(ko_lfc.mean())


def _guide_colors_and_hover(
    store: ScreenDataStore | None,
    gene: str,
    guide_positions: list[int],
) -> tuple[list[str], list[str]]:
    colors = []
    hovers = []
    overlaps = store.screen_residue_locs_by_marker(gene) if store and store.is_ready() else {}

    for pos in guide_positions:
        hits = [m for m in SCREEN_MARKERS if pos in overlaps.get(m, set())]
        if not hits:
            colors.append(_DEFAULT_GUIDE_COLOR)
            hovers.append(f'AMBER guide<br>Position {pos}<br>No screen overlap')
        elif len(hits) == 1:
            marker = hits[0]
            colors.append(MARKER_COLORS[marker])
            hovers.append(f'AMBER guide<br>Position {pos}<br>Overlaps {marker} screen')
        else:
            colors.append(MULTI_MARKER_COLOR)
            hovers.append(
                f'AMBER guide<br>Position {pos}<br>Overlaps: {", ".join(hits)}'
            )
    return colors, hovers


def _set_track_yaxis(
    fig: go.Figure,
    row: int,
    *,
    title: str,
    showticklabels: bool = False,
    y_range: tuple[float, float] = (0, 1),
    tickvals: list[float] | None = None,
) -> None:
    """Shared y-axis chrome for overview protein-map tracks."""
    kwargs = dict(
        title_text=title,
        title_font=dict(size=11),
        range=list(y_range),
        showticklabels=showticklabels,
        visible=True,
        zeroline=False,
        showgrid=False,
        row=row,
        col=1,
    )
    if showticklabels:
        kwargs['tickfont'] = dict(size=9)
        if tickvals is not None:
            kwargs['tickvals'] = tickvals
            kwargs['ticktext'] = [str(int(v) if float(v).is_integer() else v) for v in tickvals]
        elif y_range == (0, 1):
            kwargs['tickvals'] = [0, 0.5, 1]
    fig.update_yaxes(**kwargs)


def _add_domain_track(
    fig: go.Figure,
    domain_layout: list[dict],
    domain_row: int,
    plt_range: tuple[int, int],
    pro_len: int,
    lfc_vmin: float = -0.5,
    lfc_vmax: float = 0.5,
) -> None:
    _set_track_yaxis(fig, domain_row, title='Domains')
    fig.update_xaxes(range=[plt_range[0] - 1, plt_range[1] + 1], row=domain_row, col=1)

    if not domain_layout:
        fig.add_annotation(
            text='No UniProt domain annotations for this gene',
            x=(plt_range[0] + plt_range[1]) / 2,
            y=0.5,
            showarrow=False,
            font=dict(size=11, color='#666'),
            row=domain_row,
            col=1,
        )
        return

    centers, widths, colors, labels, hover = [], [], [], [], []
    for d in domain_layout:
        start, end = d['start'], d['end']
        length = end - start
        if length <= 0:
            continue
        centers.append(start + length / 2)
        widths.append(length)
        avg = d.get('avg_lfc')
        if avg is not None:
            colors.append(_lfc_to_color(float(avg), lfc_vmin, 0.0, lfc_vmax))
        else:
            colors.append(_DOMAIN_NEUTRAL)
        labels.append(d['name'])
        hover.append(f"{d['name']}<br>{start}–{end}")

    # Vertically centered band: base 0.25 + height 0.5 → spans 0.25–0.75.
    fig.add_trace(
        go.Bar(
            x=centers,
            y=[0.5] * len(centers),
            base=[0.25] * len(centers),
            width=widths,
            marker=dict(color=colors, line=dict(color='black', width=0.5)),
            text=labels,
            textposition='inside',
            insidetextanchor='middle',
            textfont=dict(size=8),
            hovertext=hover,
            hoverinfo='text',
            showlegend=False,
        ),
        row=domain_row,
        col=1,
    )


def _clinvar_color(classification: str | None) -> str:
    key = (classification or '').strip()
    if key in _CLINVAR_COLORS:
        return _CLINVAR_COLORS[key]
    lower = key.lower()
    if 'conflict' in lower:
        return _CLINVAR_COLORS['Conflicting classifications']
    if lower == 'pathogenic':
        return _CLINVAR_COLORS['Pathogenic']
    if lower == 'likely pathogenic':
        return _CLINVAR_COLORS['Likely pathogenic']
    if 'uncertain' in lower:
        return _CLINVAR_COLORS['Uncertain significance']
    if lower == 'likely benign':
        return _CLINVAR_COLORS['Likely benign']
    if lower == 'benign':
        return _CLINVAR_COLORS['Benign']
    return _CLINVAR_COLORS['Unknown']


def _assign_clinvar_lanes(
    variants: list[dict],
    *,
    max_lanes: int = _CLINVAR_MAX_LANES,
) -> list[tuple[dict, int, int]]:
    """
    Assign each variant a lane (0..max_lanes-1) for stacked display.

    Returns list of (variant, lane, overflow_count_at_residue).
    overflow_count is total variants at that residue minus those shown.
    """
    by_residue: dict[int, list[dict]] = {}
    for v in variants:
        try:
            residue = int(v.get('residue'))
        except (TypeError, ValueError):
            continue
        by_residue.setdefault(residue, []).append(v)

    assigned: list[tuple[dict, int, int]] = []
    for residue in sorted(by_residue):
        group = by_residue[residue]
        overflow = max(0, len(group) - max_lanes)
        for idx, variant in enumerate(group[:max_lanes]):
            assigned.append((variant, idx, overflow if idx == 0 else 0))
    return assigned


def _add_clinvar_track(
    fig: go.Figure,
    clinvar_variants: list[dict] | None,
    clinvar_row: int,
    plt_range: tuple[int, int],
    pro_len: int,
    *,
    empty_message: str | None = None,
) -> None:
    variants = list(clinvar_variants or [])
    assigned = _assign_clinvar_lanes(variants) if variants else []
    max_count = max((lane + 1 for _, lane, _ in assigned), default=1)
    y_top = max(max_count, 1) + 0.5
    # Label every 2nd count (2, 4, …); keep a single tick when max is 1.
    tickvals = list(range(2, max_count + 1, 2)) or ([1] if max_count >= 1 else [])
    _set_track_yaxis(
        fig,
        clinvar_row,
        title='# entries',
        showticklabels=True,
        y_range=(0.5, y_top),
        tickvals=tickvals,
    )
    fig.update_xaxes(range=[plt_range[0] - 1, plt_range[1] + 1], row=clinvar_row, col=1)

    if not assigned:
        fig.add_annotation(
            text=empty_message or 'No ClinVar protein variants mapped',
            x=(plt_range[0] + plt_range[1]) / 2,
            y=(0.5 + y_top) / 2,
            showarrow=False,
            font=dict(size=11, color='#666'),
            row=clinvar_row,
            col=1,
        )
        return

    xs, ys, colors, hovers, customdata = [], [], [], [], []
    for variant, lane, overflow in assigned:
        residue = int(variant['residue'])
        # 1-indexed entry count for stacked ClinVar ticks at this residue.
        y = lane + 1
        xs.append(residue)
        ys.append(y)
        colors.append(_clinvar_color(variant.get('classification')))
        title = variant.get('title') or 'ClinVar variant'
        protein_change = variant.get('protein_change') or '—'
        classification = variant.get('classification') or 'Unknown'
        review = variant.get('review_status') or '—'
        accession = variant.get('accession') or variant.get('uid') or '—'
        url = variant.get('url') or ''
        hover = (
            f'{title}<br>'
            f'Protein: {protein_change} (residue {residue})<br>'
            f'Classification: {classification}<br>'
            f'Review: {review}<br>'
            f'{accession}'
        )
        if url:
            hover += '<br>Click to open in ClinVar'
        if overflow:
            hover += f'<br>+{overflow} more at this site'
        hovers.append(hover)
        customdata.append(url)

    fig.add_trace(
        go.Scatter(
            x=xs,
            y=ys,
            mode='markers',
            marker=dict(
                symbol='line-ns',
                size=16,
                color=colors,
                line=dict(width=2.5, color=colors),
            ),
            customdata=customdata,
            hovertext=hovers,
            hoverinfo='text',
            showlegend=False,
        ),
        row=clinvar_row,
        col=1,
    )


def _guide_bar_heights(
    guide_positions: list[int],
    guide_scores: dict[int, float] | None,
) -> list[float]:
    """Bar height = max AM at position when scores provided; else fixed 0.5."""
    if guide_scores is None:
        return [0.5] * len(guide_positions)
    heights = []
    for pos in guide_positions:
        raw = guide_scores.get(pos)
        if raw is None:
            heights.append(0.02)
            continue
        try:
            heights.append(max(0.02, min(1.0, float(raw))))
        except (TypeError, ValueError):
            heights.append(0.02)
    return heights


def _add_guide_track(
    fig: go.Figure,
    guide_positions: list[int],
    guide_row: int,
    plt_range: tuple[int, int],
    store: ScreenDataStore | None,
    gene: str,
    *,
    guide_scores: dict[int, float] | None = None,
    y_axis_title: str = 'AlphaMissense',
) -> None:
    _set_track_yaxis(fig, guide_row, title=y_axis_title, showticklabels=True)
    fig.update_xaxes(range=[plt_range[0] - 1, plt_range[1] + 1], row=guide_row, col=1)

    if not guide_positions:
        fig.add_annotation(
            text='No guide RNAs in current results',
            x=(plt_range[0] + plt_range[1]) / 2,
            y=0.5,
            showarrow=False,
            font=dict(size=11, color='#666'),
            row=guide_row,
            col=1,
        )
        return

    colors, hovers = _guide_colors_and_hover(store, gene, guide_positions)
    heights = _guide_bar_heights(guide_positions, guide_scores)
    if guide_scores is not None:
        enriched = []
        for pos, hover, height in zip(guide_positions, hovers, heights):
            enriched.append(
                f'{hover}<br>Max AM {height:.3f}<br>Click for guide details'
            )
        hovers = enriched
    customdata = [{'type': 'guide', 'position': int(pos)} for pos in guide_positions]
    fig.add_trace(
        go.Bar(
            x=guide_positions,
            y=heights,
            width=[0.8] * len(guide_positions),
            marker=dict(color=colors),
            customdata=customdata,
            hovertext=hovers,
            hoverinfo='text',
            showlegend=False,
        ),
        row=guide_row,
        col=1,
    )


def _add_lfc_panel(
    fig: go.Figure,
    store: ScreenDataStore,
    gene: str,
    marker: str,
    be_type: str,
    row: int,
    plt_range: tuple[int, int],
    pro_len: int,
    *,
    mark_zero_coverage: bool = True,
) -> tuple[float, float]:
    df = _gene_marker_slice(store, gene, marker, be_type)
    df_plot = df.loc[plt_range[0]:plt_range[1]]
    vmin, vcenter, vmax = _diverging_norm(df_plot['LFC'])

    x_vals = df_plot.index.tolist()
    bar_colors = [_lfc_to_color(v, vmin, vcenter, vmax) for v in df_plot['LFC']]
    y_max = max(float(df_plot['LFC'].max()), 1.25)
    y_min = min(float(df_plot['LFC'].min()), -1.25)

    fig.add_trace(
        go.Bar(
            x=x_vals,
            y=df_plot['LFC'],
            marker_color=bar_colors,
            width=1.0,
            showlegend=False,
            hovertemplate=f'{marker} {be_type}<br>Position %{{x}}<br>LFC %{{y:.3f}}<extra></extra>',
        ),
        row=row,
        col=1,
    )

    ko_lfc = _ko_reference_lfc(store, gene, marker, be_type)
    if ko_lfc is not None:
        ko_color = _lfc_to_color(ko_lfc, vmin, vcenter, vmax)
        fig.add_shape(
            type='line',
            x0=0,
            x1=pro_len,
            y0=ko_lfc,
            y1=ko_lfc,
            line=dict(color=ko_color, width=1.5, dash='dash'),
            row=row,
            col=1,
        )

    if mark_zero_coverage:
        zero_idx = df_plot.index[df_plot['control_mean'] == 0.0]
        for idx in zero_idx:
            fig.add_shape(
                type='rect',
                x0=int(idx) - 0.5,
                x1=int(idx) + 0.5,
                y0=y_min - 2,
                y1=y_max + 2,
                fillcolor='rgba(180,180,180,0.5)',
                line_width=0,
                layer='below',
                row=row,
                col=1,
            )

    fig.add_hline(y=0, line_width=0.5, line_color='black', row=row, col=1)
    fig.update_yaxes(
        title_text=f'{marker} {be_type}',
        range=[y_min - 0.2, y_max + 0.2],
        row=row,
        col=1,
    )
    return vmin, vmax


def build_overview_figure(
    gene: str,
    guide_positions: list[int] | None = None,
    store: ScreenDataStore | None = None,
    *,
    protein_length: int | None = None,
    domain_layout: list[dict] | None = None,
    clinvar_variants: list[dict] | None = None,
    clinvar_empty_message: str | None = None,
    guide_scores: dict[int, float] | None = None,
) -> go.Figure:
    """Domains + ClinVar + AMBER guide tracks (shown below the results table)."""
    store = store or ScreenDataStore.get()

    if protein_length is not None:
        pro_len = int(protein_length)
    else:
        store.ensure_domains_only()
        pro_len = store.protein_length.get(gene)
        if not pro_len:
            raise KeyError(f'Unknown screen gene: {gene}')

    if domain_layout is None:
        try:
            domain_layout = store.get_structural_domain_layout(gene)
        except Exception:
            domain_layout = []

    plt_range = (0, pro_len)
    guide_positions = sorted({int(p) for p in (guide_positions or []) if p is not None})
    if guide_scores:
        guide_positions = sorted(set(guide_positions) | {int(p) for p in guide_scores})

    color_store = store if store.is_ready() else None
    guide_title = f'AMBER guides ({len(guide_positions)} positions)'
    clinvar_n = len(clinvar_variants or [])
    clinvar_title = f'ClinVar variants ({clinvar_n})'

    fig = make_subplots(
        rows=3,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.14,
        row_heights=[1, 1, 1],
        subplot_titles=['UniProt domains', clinvar_title, guide_title],
    )

    _add_domain_track(fig, domain_layout or [], 1, plt_range, pro_len)
    _add_clinvar_track(
        fig,
        clinvar_variants,
        2,
        plt_range,
        pro_len,
        empty_message=clinvar_empty_message,
    )
    _add_guide_track(
        fig,
        guide_positions,
        3,
        plt_range,
        color_store,
        gene,
        guide_scores=guide_scores,
        y_axis_title='AM score',
    )

    fig.update_xaxes(range=[plt_range[0] - 1, plt_range[1] + 1], row=1, col=1)
    fig.update_layout(
        title=f'{gene} — protein map',
        height=560,
        autosize=True,
        showlegend=False,
        margin=dict(l=70, r=20, t=55, b=40),
    )
    return fig


def default_lfc_panels() -> list[tuple[str, str]]:
    """All marker × editor LFC panels in display order."""
    return [(marker, be_type) for marker in SCREEN_MARKERS for be_type in BE_TYPES]


def parse_lfc_panel_key(key: str) -> tuple[str, str] | None:
    """Parse 'TNFa_ABE' → ('TNFa', 'ABE')."""
    raw = (key or '').strip()
    if '_' not in raw:
        return None
    marker, be_type = raw.rsplit('_', 1)
    if marker not in SCREEN_MARKERS or be_type not in BE_TYPES:
        return None
    return marker, be_type


def build_full_screen_figure(
    gene: str,
    guide_positions: list[int] | None = None,
    store: ScreenDataStore | None = None,
    *,
    lfc_panels: list[tuple[str, str]] | None = None,
    include_domains: bool = True,
    include_guides: bool = True,
) -> go.Figure:
    """Selected marker × ABE/CBE LFC panels, optional domains and AMBER guides."""
    store = store or ScreenDataStore.get()
    store.ensure_loaded()

    pro_len = store.protein_length[gene]
    plt_range = (0, pro_len)
    guide_positions = sorted({int(p) for p in (guide_positions or []) if p is not None})

    panels = list(lfc_panels) if lfc_panels is not None else default_lfc_panels()
    valid = []
    seen = set()
    for marker, be_type in panels:
        if marker not in SCREEN_MARKERS or be_type not in BE_TYPES:
            continue
        key = (marker, be_type)
        if key in seen:
            continue
        seen.add(key)
        valid.append(key)
    panels = valid

    n_rows = len(panels) + int(include_domains) + int(include_guides)
    if n_rows < 1:
        raise ValueError('Select at least one enrichment panel to display.')

    guide_title = f'AMBER guides ({len(guide_positions)} positions)'
    subplot_titles = [f'{marker} {be_type} LFC' for marker, be_type in panels]
    if include_domains:
        subplot_titles.append('UniProt domains')
    if include_guides:
        subplot_titles.append(guide_title)

    # Keep gaps modest; scale height with visible row count.
    spacing = min(0.045, 0.35 / max(n_rows - 1, 1)) if n_rows > 1 else 0.0
    fig = make_subplots(
        rows=n_rows,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=spacing,
        row_heights=[1] * n_rows,
        subplot_titles=subplot_titles,
    )

    last_vmin, last_vmax = -0.5, 0.5
    row = 1
    for marker, be_type in panels:
        vmin, vmax = _add_lfc_panel(
            fig,
            store,
            gene,
            marker,
            be_type,
            row,
            plt_range,
            pro_len,
            mark_zero_coverage=False,
        )
        last_vmin, last_vmax = vmin, vmax
        row += 1

    if include_domains:
        domain_layout = store.get_structural_domain_layout(gene)
        _add_domain_track(fig, domain_layout, row, plt_range, pro_len, last_vmin, last_vmax)
        row += 1

    if include_guides:
        _add_guide_track(fig, guide_positions, row, plt_range, store, gene)

    fig.update_xaxes(range=[plt_range[0] - 1, plt_range[1] + 1], row=1, col=1)
    fig.update_annotations(font_size=10)
    height = max(360, 140 * n_rows + 100)
    fig.update_layout(
        title=f'{gene} — screen enrichment (LFC)',
        height=height,
        autosize=False,
        showlegend=False,
        margin=dict(l=60, r=20, t=60, b=40),
        barmode='overlay',
    )
    return fig


def build_screen_enrichment_figure(
    gene: str,
    marker: str,
    store: ScreenDataStore | None = None,
    guide_positions: list[int] | None = None,
) -> go.Figure:
    """Backward-compatible alias: single-marker request still returns the full multi-marker figure."""
    return build_full_screen_figure(gene, guide_positions=guide_positions, store=store)
