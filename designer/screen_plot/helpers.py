"""Notebook helper functions for screen enrichment plots."""

from __future__ import annotations

from scipy.stats import combine_pvalues


def ext_loc_splice(s):
    if not s or len(s) < 3:
        return 'None'
    if 'ex' in s:
        return int(s.split('_')[-1][1:-1])
    if s == 'None':
        return 'None'
    return int(s[1:-1])


def is_KO(s):
    if s is None:
        return False
    if 'ex' in s or s.endswith('/'):
        return True
    if s.endswith('*') and not s.startswith('*'):
        return True
    if s.startswith('M0') and not s.endswith('M'):
        return True
    return False


def is_syn(s):
    if s is None:
        return False
    return s[0] == s[-1]


def combine_pvals(series):
    _, combined_p_value = combine_pvalues(list(series))
    return combined_p_value
