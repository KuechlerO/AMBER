from pathlib import Path

from django.test import SimpleTestCase, Client, override_settings

from saffron.pipeline import parse_mutation, apply_mutation, normalize_aa_sequence, SaffronUserError
from saffron.patho_spv import (
    annotate_guide_rows,
    clear_patho_spv_index,
    lookup_patho_spv,
)
from saffron.signalp_client import (
    parse_prediction_results,
    parse_region_gff3,
    compute_deltas,
    run_signalp,
    SignalPPrediction,
)

FIXTURES = Path(__file__).parent / 'fixtures'


class MutationParsingTests(SimpleTestCase):
    def test_parse_ok(self):
        self.assertEqual(parse_mutation('a23v'), ('A', 23, 'V'))

    def test_parse_rejects_bad(self):
        with self.assertRaises(SaffronUserError):
            parse_mutation('23V')

    def test_apply_mutation(self):
        seq = 'MABCDEFGHIJKLMNOP'
        out = apply_mutation(seq, 'A', 2, 'V')
        self.assertEqual(out[1], 'V')
        with self.assertRaises(SaffronUserError):
            apply_mutation(seq, 'X', 2, 'V')

    def test_normalize_sequence(self):
        self.assertEqual(normalize_aa_sequence('acdefghiklm'), 'ACDEFGHIKLM')
        with self.assertRaises(SaffronUserError):
            normalize_aa_sequence('SHORT')


class SignalPParserTests(SimpleTestCase):
    def test_parse_prediction_results(self):
        preds = parse_prediction_results(FIXTURES / 'prediction_results.txt')
        self.assertIn('WT', preds)
        self.assertEqual(preds['WT'].prediction, 'SP')
        self.assertAlmostEqual(preds['WT'].probabilities['SP'], 0.9710)
        self.assertEqual(preds['WT'].cs_before, 22)
        self.assertEqual(preds['m5_A5D'].prediction, 'OTHER')

    def test_parse_regions(self):
        regions = parse_region_gff3(FIXTURES / 'region_output.gff3')
        self.assertTrue(regions['WT'])
        starts = [r['start'] for r in regions['WT'] if r['feature'] == 'signal_peptide']
        self.assertEqual(starts, [1])

    def test_compute_deltas(self):
        wt = SignalPPrediction(
            seq_id='WT',
            prediction='SP',
            probabilities={'SP': 0.9, 'OTHER': 0.1},
            cs_before=22,
        )
        mut = SignalPPrediction(
            seq_id='M',
            prediction='OTHER',
            probabilities={'SP': 0.2, 'OTHER': 0.8},
            cs_before=None,
        )
        d = compute_deltas(wt, mut)
        self.assertTrue(d['sp_lost'])
        self.assertAlmostEqual(d['delta_wt_class_prob'], -0.7)

    @override_settings()
    def test_mock_signalp_runs(self):
        import os
        os.environ['SIGNALP_MOCK'] = '1'
        try:
            # Classical-ish N-term hydrophobic stretch
            seq = 'MKLLVVVILILALALALALAVSSSDDDDEEEEE' + 'G' * 40
            preds = run_signalp({'WT': seq}, organism='eukarya', want_plots=False)
            self.assertIn('WT', preds)
            self.assertEqual(preds['WT'].prediction, 'SP')
            self.assertIsNotNone(preds['WT'].sp_end)
            preds2 = run_signalp(
                {'WT': seq}, organism='eukarya', job_id='test_mock_plots', want_plots=True
            )
            self.assertTrue(preds2['WT'].plot_path)
            self.assertTrue(Path(preds2['WT'].plot_path).is_file())
        finally:
            os.environ.pop('SIGNALP_MOCK', None)


class PathoSpvCatalogueTests(SimpleTestCase):
    def setUp(self):
        clear_patho_spv_index()

    def tearDown(self):
        clear_patho_spv_index()

    @override_settings(PATHO_SPV_CSV=str(FIXTURES / 'patho_SPVs_mini.csv'))
    def test_lookup_and_annotate(self):
        hit = lookup_patho_spv('P05067', 8, 'L', 'P')
        self.assertIsNotNone(hit)
        self.assertTrue(hit['paper_pathogenic'])
        self.assertEqual(hit['clin_sig'], 'pathogenic')

        rows = [
            {'position': 8, 'wt_aa': 'L', 'mut_aa': 'P'},
            {'position': 17, 'wt_aa': 'A', 'mut_aa': 'T'},
            {'position': 99, 'wt_aa': 'G', 'mut_aa': 'A'},
        ]
        annotate_guide_rows(rows, accession='P05067')
        self.assertTrue(rows[0]['paper_pathogenic'])
        self.assertTrue(rows[0]['in_patho_catalogue'])
        self.assertFalse(rows[1]['paper_pathogenic'])
        self.assertTrue(rows[1]['in_patho_catalogue'])
        self.assertEqual(rows[1]['clin_sig'], 'uncertain_significance')
        self.assertFalse(rows[2]['in_patho_catalogue'])


class SaffronViewSmokeTests(SimpleTestCase):
    def test_home_ok(self):
        c = Client()
        r = c.get('/saffron/')
        self.assertEqual(r.status_code, 200)
        self.assertContains(r, 'SAFFRON')

    def test_about_ok(self):
        c = Client()
        r = c.get('/saffron/about/')
        self.assertEqual(r.status_code, 200)
