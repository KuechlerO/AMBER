"""Tests for gene symbol → UniProt search helper and JSON endpoint."""

from unittest.mock import MagicMock, patch

from django.test import Client, SimpleTestCase, TestCase

from designer.uniprot_gene_search import (
    GeneSearchError,
    annotate_am_availability,
    looks_like_gene_symbol,
    search_uniprot_by_gene_symbol,
)


def _fake_uniprot_payload(entries):
    return {'results': entries}


class LooksLikeGeneSymbolTest(SimpleTestCase):
    def test_accepts_symbols(self):
        self.assertTrue(looks_like_gene_symbol('TRBC1'))
        self.assertTrue(looks_like_gene_symbol('hla-a'))
        self.assertTrue(looks_like_gene_symbol('CTSH'))

    def test_rejects_accessions_and_ensembl(self):
        self.assertFalse(looks_like_gene_symbol('P04439'))
        self.assertFalse(looks_like_gene_symbol('A0A5B9'))
        self.assertFalse(looks_like_gene_symbol('ENST00000376809'))
        self.assertFalse(looks_like_gene_symbol(''))


class AnnotateAmAvailabilityTest(SimpleTestCase):
    @patch('designer.uniprot_gene_search.alphamissense_accessions_present')
    def test_flags_and_sorts_am_first_within_reviewed(self, mock_present):
        mock_present.return_value = {'P04439'}
        rows = annotate_am_availability([
            {
                'accession': 'A0A0G2JMB4',
                'gene': 'TRBC2',
                'protein_name': 'x',
                'length': 179,
                'reviewed': False,
            },
            {
                'accession': 'P04439',
                'gene': 'HLA-A',
                'protein_name': 'y',
                'length': 365,
                'reviewed': True,
            },
            {
                'accession': 'P01850',
                'gene': 'TRBC1',
                'protein_name': 'z',
                'length': 176,
                'reviewed': True,
            },
        ])
        self.assertTrue(rows[0]['am_available'])
        self.assertEqual(rows[0]['accession'], 'P04439')
        self.assertFalse(rows[1]['am_available'])
        self.assertEqual(rows[1]['accession'], 'P01850')
        self.assertFalse(rows[2]['am_available'])


class SearchUniprotByGeneSymbolTest(SimpleTestCase):
    def test_rejects_uniprot_shaped_input(self):
        with self.assertRaises(GeneSearchError):
            search_uniprot_by_gene_symbol('P04439')

    def test_empty_raises(self):
        with self.assertRaises(GeneSearchError):
            search_uniprot_by_gene_symbol('  ')

    @patch('designer.uniprot_gene_search.alphamissense_accessions_present', return_value=set())
    @patch('designer.uniprot_gene_search._http_get')
    def test_parses_and_sorts_reviewed_first(self, mock_get, _mock_am):
        response = MagicMock()
        response.raise_for_status = MagicMock()
        response.json.return_value = _fake_uniprot_payload([
            {
                'entryType': 'UniProtKB unreviewed (TrEMBL)',
                'primaryAccession': 'A0A5H1ZRT1',
                'proteinDescription': {
                    'submissionNames': [{'fullName': {'value': 'Fragment'}}],
                },
                'genes': [{'geneName': {'value': 'TRBC1'}}],
                'sequence': {'length': 177},
            },
            {
                'entryType': 'UniProtKB reviewed (Swiss-Prot)',
                'primaryAccession': 'P01850',
                'proteinDescription': {
                    'recommendedName': {'fullName': {'value': 'T cell receptor beta constant 1'}},
                },
                'genes': [{'geneName': {'value': 'TRBC1'}}],
                'sequence': {'length': 176},
            },
        ])
        mock_get.return_value = response

        rows = search_uniprot_by_gene_symbol('TRBC1')
        self.assertEqual([r['accession'] for r in rows], ['P01850', 'A0A5H1ZRT1'])
        self.assertTrue(rows[0]['reviewed'])
        self.assertFalse(rows[1]['reviewed'])
        self.assertFalse(rows[0]['am_available'])
        self.assertEqual(rows[0]['gene'], 'TRBC1')
        self.assertEqual(rows[0]['length'], 176)

    @patch('designer.uniprot_gene_search.alphamissense_accessions_present', return_value=set())
    @patch('designer.uniprot_gene_search._http_get')
    def test_falls_back_when_exact_empty(self, mock_get, _mock_am):
        empty = MagicMock()
        empty.raise_for_status = MagicMock()
        empty.json.return_value = {'results': []}

        filled = MagicMock()
        filled.raise_for_status = MagicMock()
        filled.json.return_value = _fake_uniprot_payload([
            {
                'entryType': 'UniProtKB reviewed (Swiss-Prot)',
                'primaryAccession': 'P09668',
                'proteinDescription': {
                    'recommendedName': {'fullName': {'value': 'Cathepsin H'}},
                },
                'genes': [{'geneName': {'value': 'CTSH'}}],
                'sequence': {'length': 335},
            },
        ])
        mock_get.side_effect = [empty, filled]

        rows = search_uniprot_by_gene_symbol('CTSH')
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]['accession'], 'P09668')
        self.assertEqual(mock_get.call_count, 2)

    @patch('designer.uniprot_gene_search.alphamissense_accessions_present', return_value=set())
    @patch('designer.uniprot_gene_search._http_get')
    def test_empty_results(self, mock_get, _mock_am):
        response = MagicMock()
        response.raise_for_status = MagicMock()
        response.json.return_value = {'results': []}
        mock_get.return_value = response
        self.assertEqual(search_uniprot_by_gene_symbol('NOGENEXYZ'), [])


class GeneUniprotSearchViewTest(TestCase):
    @patch('designer.views.search_uniprot_by_gene_symbol')
    def test_amber_endpoint_ok(self, mock_search):
        mock_search.return_value = [{
            'accession': 'P01850',
            'gene': 'TRBC1',
            'protein_name': 'T cell receptor beta constant 1',
            'length': 176,
            'reviewed': True,
            'am_available': False,
        }]
        c = Client()
        r = c.get('/api/gene-uniprot/', {'q': 'TRBC1'})
        self.assertEqual(r.status_code, 200)
        data = r.json()
        self.assertEqual(data['results'][0]['accession'], 'P01850')
        self.assertFalse(data['results'][0]['am_available'])

    @patch('designer.views.search_uniprot_by_gene_symbol')
    def test_amber_endpoint_bad_query(self, mock_search):
        mock_search.side_effect = GeneSearchError('Invalid gene symbol.')
        c = Client()
        r = c.get('/api/gene-uniprot/', {'q': '!!!'})
        self.assertEqual(r.status_code, 400)
        self.assertIn('error', r.json())

    @patch('saffron.views.search_uniprot_by_gene_symbol')
    def test_saffron_endpoint_ok(self, mock_search):
        mock_search.return_value = [{
            'accession': 'P09668',
            'gene': 'CTSH',
            'protein_name': 'Cathepsin H',
            'length': 335,
            'reviewed': True,
            'am_available': True,
        }]
        c = Client()
        r = c.get('/saffron/api/gene-uniprot/', {'q': 'CTSH'})
        self.assertEqual(r.status_code, 200)
        self.assertEqual(r.json()['results'][0]['accession'], 'P09668')


class HomePickerMarkupTest(SimpleTestCase):
    def test_amber_home_has_picker(self):
        c = Client()
        r = c.get('/')
        self.assertEqual(r.status_code, 200)
        self.assertContains(r, 'gene-picker-modal')
        self.assertContains(r, 'gene-lookup-btn')
        self.assertContains(r, 'gene_picker.js')
        self.assertContains(r, 'btn-form-aux')

    def test_saffron_home_has_picker(self):
        c = Client()
        r = c.get('/saffron/')
        self.assertEqual(r.status_code, 200)
        self.assertContains(r, 'gene-picker-modal')
        self.assertContains(r, 'Look up gene')


class EmptyAlphamissensePipelineTest(SimpleTestCase):
    @patch('designer.pipeline.get_exon_boundaries')
    @patch('designer.pipeline.get_alphamissense_from_db')
    @patch('designer.pipeline.get_cds_and_transcript_from_uniprot')
    def test_raises_when_am_missing(self, mock_cds, mock_am, mock_exons):
        import pandas as pd
        from designer.pipeline import run_pipeline

        mock_cds.return_value = ('ATG' * 10, 'ENST00000000001', 'TRBC2')
        mock_am.return_value = pd.DataFrame(
            columns=['protein_variant', 'pathogenicity score', 'a.a.1', 'position', 'a.a.2']
        )
        with self.assertRaises(ValueError) as ctx:
            run_pipeline(
                uniprot_id='A0A0G2JMB4',
                editor='ABE',
                alpha_threshold=0.0,
                top_sgrnas=5,
                window_min=4,
                window_max=8,
            )
        self.assertIn('No AlphaMissense scores found', str(ctx.exception))
        mock_exons.assert_not_called()
