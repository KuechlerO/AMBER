import json
from unittest.mock import MagicMock

import numpy as np
import pandas as pd
from django.test import TestCase

from designer.screen_plot.plotly_builder import (
    build_full_screen_figure,
    build_overview_figure,
)


class PlotlyBuilderSmokeTest(TestCase):
    def _mock_store(self, domain_layout=None):
        store = MagicMock()
        store.protein_length = {'TESTGENE': 5}
        store.is_ready = MagicMock(return_value=True)
        store.data_cutoff = pd.DataFrame({
            'marker': [], 'be_type': [], 'gene': [], 'sites_mapped': [], 'sgrna': [], 'LFC': [],
        })
        idx = pd.MultiIndex.from_product(
            [['TESTGENE'], ['TNFa', 'IFNG', 'PD1', 'CD25'], ['ABE', 'CBE'], [0, 1, 2, 3, 4]],
            names=['gene', 'marker', 'be_type', 'loc'],
        )
        df = pd.DataFrame({
            'LFC': [0.1, -0.2, 0.3, -0.1, 0.0] * 8,
            'control_mean': [10.0] * 40,
            'p.twosided': [0.05] * 40,
            'effect_size_base': [None] * 40,
            'SE_base': [None] * 40,
            'p-value_base': [None] * 40,
            'codon_pos': [None] * 40,
        }, index=idx)
        store.data_grouped = df
        store.get_structural_domain_layout = MagicMock(
            return_value=domain_layout
            if domain_layout is not None
            else [{'name': 'Test domain', 'start': 0, 'end': 3, 'avg_lfc': None, 'domain_type': 'Domain'}]
        )
        store.screen_residue_locs_by_marker = MagicMock(return_value={
            'TNFa': {1},
            'IFNG': set(),
            'PD1': {3},
            'CD25': set(),
        })
        store.ensure_loaded = MagicMock()
        store.ensure_domains_only = MagicMock()
        return store

    def test_overview_two_rows(self):
        fig = build_overview_figure('TESTGENE', guide_positions=[1, 3], store=self._mock_store())
        payload = json.loads(fig.to_json())
        titles = [a.get('text', '') for a in payload['layout'].get('annotations', [])]
        self.assertIn('UniProt domains', titles)
        self.assertTrue(any('ClinVar variants' in (t or '') for t in titles))
        self.assertFalse(any('AMBER guides' in (t or '') for t in titles))
        xaxis3 = payload['layout'].get('xaxis3') or {}
        x_title = xaxis3.get('title')
        if isinstance(x_title, dict):
            x_title = x_title.get('text', '')
        self.assertIn('AMBER guides', x_title or '')

    def test_overview_clinvar_markers(self):
        variants = [
            {
                'uid': '1',
                'accession': 'VCV1',
                'title': 'p.Leu2Val',
                'protein_change': 'L2V',
                'residue': 2,
                'classification': 'Pathogenic',
                'review_status': 'criteria provided',
                'url': 'https://www.ncbi.nlm.nih.gov/clinvar/variation/1/',
            },
            {
                'uid': '2',
                'accession': 'VCV2',
                'title': 'p.Arg2Cys',
                'protein_change': 'R2C',
                'residue': 2,
                'classification': 'Uncertain significance',
                'review_status': 'criteria provided',
                'url': 'https://www.ncbi.nlm.nih.gov/clinvar/variation/2/',
            },
        ]
        fig = build_overview_figure(
            'TESTGENE',
            guide_positions=[1],
            store=self._mock_store(),
            clinvar_variants=variants,
        )
        payload = json.loads(fig.to_json())
        scatter = [t for t in payload['data'] if t.get('type') == 'scatter']
        self.assertEqual(len(scatter), 1)
        self.assertEqual(sorted(scatter[0]['x']), [2, 2])
        colors = scatter[0]['marker']['color']
        self.assertIn('#b91c1c', colors)
        self.assertIn('#2563eb', colors)
        custom = scatter[0].get('customdata') or []
        self.assertEqual(len(custom), 2)
        self.assertTrue(all(isinstance(u, str) and u.startswith('http') for u in custom))
        hover = scatter[0].get('hovertext') or []
        self.assertTrue(any('Click to open in ClinVar' in (h or '') for h in hover))

    def test_overview_empty_clinvar_still_builds(self):
        fig = build_overview_figure(
            'TESTGENE',
            guide_positions=[1],
            store=self._mock_store(),
            clinvar_variants=[],
            clinvar_empty_message='No ClinVar protein variants mapped',
        )
        payload = json.loads(fig.to_json())
        self.assertTrue(payload['data'] or payload['layout'])
        titles = [a.get('text', '') for a in payload['layout'].get('annotations', [])]
        self.assertTrue(any('ClinVar variants' in (t or '') for t in titles))
        self.assertTrue(any('No ClinVar protein variants mapped' in (t or '') for t in titles))

    def test_full_screen_lfc_only_ten_rows(self):
        fig = build_full_screen_figure('TESTGENE', guide_positions=[1], store=self._mock_store())
        payload = json.loads(fig.to_json())
        titles = [a.get('text', '') for a in payload['layout'].get('annotations', [])]
        self.assertEqual(len(titles), 9)
        self.assertTrue(any('TNFa ABE LFC' in t for t in titles))
        self.assertIn('UniProt domains', titles)
        self.assertFalse(any('AMBER guides' in (t or '') for t in titles))
        self.assertFalse(any('log(p)' in t for t in titles))
        self.assertFalse(any('effect size' in t for t in titles))
        xaxis10 = payload['layout'].get('xaxis10') or {}
        x_title = xaxis10.get('title')
        if isinstance(x_title, dict):
            x_title = x_title.get('text', '')
        self.assertIn('AMBER guides', x_title or '')

    def test_guide_overlap_colors(self):
        store = self._mock_store()
        fig = build_overview_figure('TESTGENE', guide_positions=[1, 3], store=store)
        payload = json.loads(fig.to_json())
        guide_traces = [
            t for t in payload['data']
            if t.get('type') == 'bar' and list(t.get('x', [])) == [1, 3]
        ]
        self.assertEqual(len(guide_traces), 1)
        colors = guide_traces[0]['marker']['color']
        self.assertEqual(colors[0], '#22c55e')
        self.assertEqual(colors[1], '#3b82f6')
