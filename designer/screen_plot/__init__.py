"""Screen enrichment plots (NGG library, Supplementary Tables 1–6)."""

from .gene_lookup import (
    SCREEN_MARKERS,
    lookup_gene_by_uniprot,
    normalize_uniprot_accession,
    resolve_lookup_uniprot_id,
    screen_plot_eligibility,
)
from .data_loader import get_screen_data_store
from .plotly_builder import (
    MARKER_COLORS,
    build_full_screen_figure,
    build_overview_figure,
    build_screen_enrichment_figure,
)

__all__ = [
    'SCREEN_MARKERS',
    'MARKER_COLORS',
    'lookup_gene_by_uniprot',
    'normalize_uniprot_accession',
    'resolve_lookup_uniprot_id',
    'screen_plot_eligibility',
    'get_screen_data_store',
    'build_overview_figure',
    'build_full_screen_figure',
    'build_screen_enrichment_figure',
]
