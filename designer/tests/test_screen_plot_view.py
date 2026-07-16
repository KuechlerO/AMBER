import json
from pathlib import Path
from unittest.mock import MagicMock, patch

from django.test import Client, SimpleTestCase, TestCase, override_settings

SCREEN_PLOT_PATH = '/results/screen-plot/'
SCREEN_OVERVIEW_PATH = '/results/screen-plot/overview/'


class ScreenPlotViewTest(TestCase):
    def setUp(self):
        self.client = Client()

    @patch('designer.views.get_form_data')
    @patch('designer.views.resolve_domain_layout')
    @patch('designer.views.screen_plot_eligibility')
    def test_overview_returns_figure_without_library(
        self, mock_eligibility, mock_domains, mock_form_data
    ):
        mock_eligibility.return_value = {
            'available': False,
            'meta': None,
            'reason': 'not_in_library',
        }
        mock_domains.return_value = ([], 'none')
        mock_form_data.return_value = {
            'uniprot_id': 'P99999',
            'uniprot_accession': 'P99999',
            'protein_length': 120,
            'gene_symbol': 'TESTGENE',
        }
        mock_fig = MagicMock()
        mock_fig.to_json.return_value = json.dumps({
            'data': [],
            'layout': {'title': {'text': 'TESTGENE'}},
        })
        with patch('designer.views.build_overview_figure', return_value=mock_fig):
            with patch('designer.views.get_filtered_rows', return_value=[{'position': 1, 'sgrna_seq': 'A'}]):
                with patch('designer.views.get_screen_data_store') as mock_store:
                    store = MagicMock()
                    store.is_ready.return_value = False
                    mock_store.return_value = store
                    resp = self.client.post(
                        SCREEN_OVERVIEW_PATH,
                        {'uniprot_id': 'P99999'},
                    )
        self.assertEqual(resp.status_code, 200)
        data = resp.json()
        self.assertTrue(data['available'])
        self.assertIn('figure', data)

    @patch('designer.views.screen_plot_eligibility')
    def test_full_plot_returns_202_while_loading(self, mock_eligibility):
        mock_eligibility.return_value = {
            'available': True,
            'meta': {'gene': 'TESTGENE'},
            'reason': None,
        }
        mock_store = MagicMock()
        mock_store.is_ready.return_value = False
        mock_store.is_loading.return_value = True
        mock_store.start_background_load.return_value = False
        with patch('designer.views.get_screen_data_store', return_value=mock_store):
            resp = self.client.post(SCREEN_PLOT_PATH, {'uniprot_id': 'Q15116'})
        self.assertEqual(resp.status_code, 202)
        self.assertTrue(resp.json()['loading'])


class StructurePdbViewTest(SimpleTestCase):
    def setUp(self):
        self.client = Client()

    def test_rejects_path_traversal(self):
        resp = self.client.get('/results/structure/../secret.pdb')
        self.assertIn(resp.status_code, (404, 301, 302))

    @override_settings(SCREEN_DATA_DIR=Path('/tmp/amber-missing-screen-data'))
    @patch('designer.protein_map.resolve_alphafold_pdb_url')
    @patch('requests.get')
    def test_proxies_alphafold_when_local_missing(self, mock_get, mock_resolve):
        mock_resolve.return_value = 'https://alphafold.ebi.ac.uk/files/AF-P04439-F1-model_v6.pdb'
        mock_resp = MagicMock()
        mock_resp.status_code = 200
        mock_resp.content = b'ATOM      1  N   ALA A   1\n'
        mock_resp.raise_for_status = MagicMock()
        mock_get.return_value = mock_resp
        resp = self.client.get('/results/structure/P04439.pdb')
        self.assertEqual(resp.status_code, 200)
        self.assertIn(b'ATOM', resp.content)
        mock_get.assert_called()
