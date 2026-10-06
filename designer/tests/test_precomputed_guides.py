"""Tests for precomputed Parquet guide lookup."""

import json
from pathlib import Path
from unittest.mock import patch

import pandas as pd
from django.test import SimpleTestCase, override_settings

from designer.precomputed_guides import (
    _dataframe_to_pipeline_rows,
    _guide_intersects_window,
    try_load_precomputed_guides,
)


class PrecomputedFilterTests(SimpleTestCase):
    def test_window_intersection(self):
        self.assertTrue(_guide_intersects_window([4, 5], 4, 8))
        self.assertFalse(_guide_intersects_window([2], 4, 8))
        self.assertTrue(_guide_intersects_window([10], 1, 10))

    def test_dataframe_filters_editor_pam_threshold(self):
        df = pd.DataFrame([
            {
                'position': 1,
                'wt_aa': 'A',
                'wt_codon': 'GCA',
                'mut_aa': 'V',
                'alpha_score': 0.9,
                'avg_alpha_score': 0.9,
                'sgrna_seq': 'A' * 23,
                'protospacer': 'A' * 20,
                'pam': 'NGG',
                'strand': '+',
                'target_positions': '[5]',
                'num_editable_window': 1,
                'editor_used': 'ABE',
                'pam_type': 'NGG',
                'outcome_details': '[{"codon":"GTA","aa":"V","score":0.9,"is_pathogenic":true}]',
            },
            {
                'position': 2,
                'wt_aa': 'C',
                'wt_codon': 'TGC',
                'mut_aa': 'Y',
                'alpha_score': 0.1,
                'avg_alpha_score': 0.1,
                'sgrna_seq': 'C' * 23,
                'protospacer': 'C' * 20,
                'pam': 'NGG',
                'strand': '+',
                'target_positions': '[6]',
                'num_editable_window': 1,
                'editor_used': 'CBE',
                'pam_type': 'NGG',
                'outcome_details': '[{"codon":"TAC","aa":"Y","score":0.1,"is_pathogenic":true}]',
            },
        ])
        rows = _dataframe_to_pipeline_rows(
            df,
            editor='ABE',
            pam_type='NGG',
            window_min=4,
            window_max=8,
            alpha_threshold=0.5,
            top_sgrnas=5,
            uniprot_id='P00001',
        )
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]['editor_used'], 'ABE')


@override_settings(PRECOMPUTED_GUIDES_DISABLE=False)
class PrecomputedLookupIntegrationTests(SimpleTestCase):
    def test_loads_real_shard_when_present(self):
        root = Path(__file__).resolve().parents[2]
        guides_dir = root / 'data' / 'precomputed_guides' / 'guides'
        if not guides_dir.is_dir():
            self.skipTest('No precomputed guides on disk')
        shard = next(guides_dir.glob('P04439_*.parquet'), None)
        if shard is None:
            self.skipTest('P04439 shard missing')

        with override_settings(PRECOMPUTED_GUIDES_DIR=root / 'data' / 'precomputed_guides'):
            with patch('designer.precomputed_guides.lookup_esm1b', return_value=None):
                result = try_load_precomputed_guides(
                    'P04439',
                    'BOTH',
                    '0.0',
                    5,
                    4,
                    8,
                    'NGG',
                )
        self.assertIsNotNone(result)
        self.assertTrue(result.get('precomputed_source'))
        self.assertGreater(len(result['guide_rows']), 0)
        self.assertIn('outcomes', result['guide_rows'][0])

    def test_returns_none_when_window_outside_precomputed_range(self):
        root = Path(__file__).resolve().parents[2]
        guides_dir = root / 'data' / 'precomputed_guides' / 'guides'
        if not guides_dir.is_dir():
            self.skipTest('No precomputed guides on disk')
        with override_settings(PRECOMPUTED_GUIDES_DIR=root / 'data' / 'precomputed_guides'):
            result = try_load_precomputed_guides(
                'P04439',
                'ABE',
                '0.5',
                5,
                11,
                15,
                'NGG',
            )
        self.assertIsNone(result)

    def test_returns_none_when_disabled(self):
        with override_settings(PRECOMPUTED_GUIDES_DISABLE=True):
            result = try_load_precomputed_guides(
                'P04439',
                'ABE',
                '0.5',
                5,
                4,
                8,
                'NGG',
            )
        self.assertIsNone(result)
