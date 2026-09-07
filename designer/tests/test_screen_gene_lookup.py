from django.test import TestCase

from designer.screen_plot.gene_lookup import (
    list_screen_library_genes,
    lookup_gene_by_uniprot,
    normalize_uniprot_accession,
    screen_plot_eligibility,
)


class ScreenGeneLookupTest(TestCase):
    def test_normalize_isoform_suffix(self):
        self.assertEqual(normalize_uniprot_accession('P60709-1'), 'P60709')

    def test_lookup_pdcd1(self):
        meta = lookup_gene_by_uniprot('Q15116')
        if meta is None:
            self.skipTest('Supplementary Table 1 not available in test environment')
        self.assertEqual(meta['gene'], 'PDCD1')

    def test_list_screen_library_genes(self):
        genes = list_screen_library_genes()
        if not genes:
            self.skipTest('Supplementary Table 1 not available in test environment')
        self.assertGreaterEqual(len(genes), 100)
        self.assertIn('gene', genes[0])
        self.assertIn('uniprot_accession', genes[0])
        symbols = [g['gene'] for g in genes]
        self.assertEqual(symbols, sorted(symbols, key=str.upper))
        self.assertIn('PDCD1', symbols)

    def test_eligibility_not_in_library(self):
        status = screen_plot_eligibility('P99999')
        if status['reason'] == 'data_unavailable':
            self.skipTest('screen data unreadable')
        self.assertFalse(status['available'])
        self.assertEqual(status['reason'], 'not_in_library')

    def test_screen_library_page(self):
        response = self.client.get('/screen-library/')
        self.assertEqual(response.status_code, 200)
        if response.context and response.context.get('data_available'):
            self.assertContains(response, 'NGG screen library genes')
            self.assertContains(response, 'PDCD1')
