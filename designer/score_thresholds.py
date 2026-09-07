"""Author-recommended pathogenicity / significance cutoffs for AMBER scores.

Display / marking only — does not change guide generation or the user AlphaMissense
filter slider on the home page.
"""

from __future__ import annotations

from typing import Any

# Cheng et al., Science 2023 — AlphaMissense class thresholds (90% ClinVar precision)
AM_LIKELY_PATHOGENIC = 0.564
AM_LIKELY_BENIGN = 0.34

# Brandes et al., Nat Genet 2023 — ESM1b LLR pathogenic/benign cutoff
ESM1B_DAMAGING = -7.5

# CADD PHRED: top 1% of SNVs (common interpretive mark; CADD authors advise
# against a single universal clinical cutoff)
CADD_SIGNIFICANT = 20.0


def _as_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def am_likely_pathogenic(score: Any) -> bool:
    """True if AlphaMissense score is in the likely_pathogenic class (≥ 0.564)."""
    val = _as_float(score)
    return val is not None and val >= AM_LIKELY_PATHOGENIC


def esm1b_damaging(llr: Any) -> bool:
    """True if ESM1b LLR is at or below the author damaging cutoff (≤ −7.5)."""
    val = _as_float(llr)
    return val is not None and val <= ESM1B_DAMAGING


def cadd_significant(phred: Any) -> bool:
    """True if CADD PHRED is at or above the top-1% mark (≥ 20)."""
    val = _as_float(phred)
    return val is not None and val >= CADD_SIGNIFICANT


def threshold_legend_context() -> dict[str, Any]:
    """Values for results / home templates."""
    return {
        'am_likely_pathogenic_threshold': AM_LIKELY_PATHOGENIC,
        'esm1b_damaging_threshold': ESM1B_DAMAGING,
        'cadd_significant_threshold': CADD_SIGNIFICANT,
    }
