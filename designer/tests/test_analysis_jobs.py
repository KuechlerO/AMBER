"""Tests for AMBER async analysis jobs, AM threshold filter, and deferred CADD."""

import time
from unittest.mock import patch

import pandas as pd
from django.core.cache import cache
from django.test import Client, SimpleTestCase, TestCase

from designer.analysis_jobs import (
    clear_job_state,
    get_status,
    start_analysis,
)
from designer.pipeline import enrich_guide_rows_with_cadd, get_alphamissense_from_db, run_pipeline
from designer.result_store import analysis_cache_key, load_analysis_results_by_key
from designer.services import UserInputError


class AlphamissenseThresholdFilterTests(SimpleTestCase):
    @patch('designer.pipeline.Alpha_missense')
    def test_filters_scores_below_threshold(self, mock_am):
        mock_am.objects.using.return_value.filter.return_value.values_list.return_value = [
            ('A10B', 0.2),
            ('A11C', 0.8),
            ('A12D', 0.5),
        ]
        df = get_alphamissense_from_db('P12345', 0.5)
        self.assertEqual(len(df), 2)
        self.assertTrue((df['pathogenicity score'] >= 0.5).all())

    @patch('designer.pipeline.get_exon_boundaries')
    @patch('designer.pipeline.get_alphamissense_from_db')
    @patch('designer.pipeline.get_cds_and_transcript_from_uniprot')
    def test_raises_when_nothing_meets_threshold(self, mock_cds, mock_am, mock_exons):
        mock_cds.return_value = ('ATG' * 10, 'ENST00000000001', 'TEST')
        empty = pd.DataFrame(
            columns=['protein_variant', 'pathogenicity score', 'a.a.1', 'position', 'a.a.2']
        )
        raw = pd.DataFrame({
            'protein_variant': ['A1B'],
            'pathogenicity score': [0.1],
            'a.a.1': ['A'],
            'position': [1],
            'a.a.2': ['B'],
        })
        mock_am.side_effect = [empty, raw]
        with self.assertRaises(ValueError) as ctx:
            run_pipeline(
                uniprot_id='P12345',
                editor='ABE',
                alpha_threshold=0.5,
                top_sgrnas=5,
                window_min=4,
                window_max=8,
            )
        self.assertIn('meet the score threshold', str(ctx.exception))
        mock_exons.assert_not_called()


class DeferredCaddPipelineTests(SimpleTestCase):
    @patch('designer.pipeline.annotate_sgrna_binding', side_effect=lambda rows: rows)
    @patch('designer.pipeline.generate_candidates')
    @patch('designer.pipeline.get_exon_boundaries')
    @patch('designer.pipeline.get_alphamissense_from_db')
    @patch('designer.pipeline.get_cds_and_transcript_from_uniprot')
    @patch('designer.pipeline.fetch_cadd_scores')
    def test_default_skips_cadd_fetch(
        self, mock_cadd, mock_cds, mock_am, mock_exons, mock_gen, mock_bind,
    ):
        mock_cds.return_value = ('ATG' * 10, 'ENST00000000001', 'TEST')
        mock_am.return_value = pd.DataFrame({
            'protein_variant': ['M1K'],
            'pathogenicity score': [0.9],
            'a.a.1': ['M'],
            'position': [1],
            'a.a.2': ['K'],
        })
        mock_exons.return_value = {'exon_regions': [], 'exon_boundaries': []}
        mock_gen.return_value = pd.DataFrame([{
            'position': 1,
            'WT_codon': 'ATG',
            'WT_AA': 'M',
            'ABE_mutated_AA': 'K',
            'ABE_mutated_AA_score': 0.9,
            'sgRNA_1_seq': 'A' * 20 + 'AGG',
            'sgRNA_1_protospacer': 'A' * 20,
            'sgRNA_1_pam': 'AGG',
            'sgRNA_1_strand': '+',
            'sgRNA_1_TargetApos': [5],
            'sgRNA_1_numAsWindow': 1,
            'sgRNA_1_outcome_details': [{
                'codon': 'AAG', 'aa': 'K', 'score': 0.9, 'is_pathogenic': True,
                'esm1b': -1.0, 'snvs': [('1', 10, 'A', 'G')],
            }],
            'sgRNA_1_local_mut_aa': 'K',
            'sgRNA_1_local_am_score': 0.9,
            'sgRNA_1_avg_alpha_score': 0.9,
        }])
        result = run_pipeline(
            uniprot_id='P12345',
            editor='ABE',
            alpha_threshold=0.5,
            top_sgrnas=5,
            window_min=4,
            window_max=8,
        )
        mock_cadd.assert_not_called()
        self.assertEqual(result['cadd_status'], 'pending')
        self.assertTrue(result['guide_rows'])
        self.assertIn('outcome_details', result['guide_rows'][0])

    def test_enrich_guide_rows_with_cadd(self):
        rows = [{
            'mut_aa': 'K',
            'alpha_score': 0.9,
            'outcome_details': [{
                'codon': 'AAG', 'aa': 'K', 'score': 0.9, 'is_pathogenic': True,
                'esm1b': -1.0, 'snvs': [('1', 10, 'A', 'G')],
            }],
        }]
        with patch('designer.pipeline.fetch_cadd_scores') as mock_cadd:
            mock_cadd.return_value = (
                {('1', 10, 'A', 'G'): {'phred': 22.0, 'raw': 1.0}},
                None,
            )
            enriched, warning = enrich_guide_rows_with_cadd(rows)
        self.assertIsNone(warning)
        self.assertAlmostEqual(enriched[0]['cadd_phred'], 22.0)
        self.assertIn('outcome_details', enriched[0])


class AnalysisJobsTests(TestCase):
    session_key = 'test-amber-session'

    def setUp(self):
        clear_job_state(self.session_key)
        cache.delete(analysis_cache_key(self.session_key))

    def tearDown(self):
        clear_job_state(self.session_key)
        cache.delete(analysis_cache_key(self.session_key))

    def test_idle_when_no_job(self):
        self.assertEqual(get_status(self.session_key)['status'], 'idle')

    @patch('designer.analysis_jobs._start_cadd_enrichment')
    @patch('designer.analysis_jobs.run_analysis')
    def test_job_completes_and_saves_results(self, mock_run, mock_cadd_start):
        mock_run.return_value = {
            'guide_rows': [{
                'position': 1,
                'alpha_score': 0.9,
                'sgrna_seq': 'A' * 23,
                'outcome_details': [],
            }],
            'no_guide_positions': [],
            'gene_symbol': 'TEST',
            'uniprot_accession': 'P12345',
            'protein_length': 10,
            'transcript_id': 'ENST1',
            'cadd_warning': 'CADD pending',
            'cadd_status': 'pending',
        }
        self.assertTrue(start_analysis(self.session_key, {
            'uniprot_id': 'P12345',
            'editor': 'ABE',
            'alpha_threshold': '0.5',
            'top_sgrnas': '5',
            'window_min': '4',
            'window_max': '8',
            'pam_type': 'NGG',
            'selected_columns': ['position', 'alpha_score'],
        }))
        self.assertEqual(get_status(self.session_key)['status'], 'running')

        for _ in range(100):
            if get_status(self.session_key)['status'] == 'done':
                break
            time.sleep(0.05)

        self.assertEqual(get_status(self.session_key)['status'], 'done')
        stored = load_analysis_results_by_key(self.session_key)
        self.assertTrue(stored.get('full_rows'))
        self.assertEqual(stored.get('cadd_status'), 'pending')
        mock_cadd_start.assert_called_once_with(self.session_key)

    @patch('designer.analysis_jobs.run_analysis')
    def test_job_user_error(self, mock_run):
        mock_run.side_effect = UserInputError('Invalid UniProt ID/ Ensembl provided.')
        start_analysis(self.session_key, {'uniprot_id': 'BAD', 'editor': 'ABE'})

        for _ in range(100):
            status = get_status(self.session_key)
            if status['status'] == 'error':
                break
            time.sleep(0.05)

        status = get_status(self.session_key)
        self.assertEqual(status['status'], 'error')
        self.assertIn('Invalid', status['error'])

    @patch('designer.analysis_jobs.run_analysis')
    def test_second_start_while_running_is_noop(self, mock_run):
        def slow_run(**_kwargs):
            time.sleep(0.3)
            return {
                'guide_rows': [],
                'no_guide_positions': [],
                'cadd_status': 'pending',
            }

        mock_run.side_effect = slow_run
        self.assertTrue(start_analysis(self.session_key, {
            'uniprot_id': 'P12345',
            'editor': 'ABE',
            'selected_columns': ['position'],
        }))
        self.assertFalse(start_analysis(self.session_key, {
            'uniprot_id': 'P12345',
            'editor': 'ABE',
            'selected_columns': ['position'],
        }))


class AmberAsyncViewTests(TestCase):
    def test_loading_post_returns_polling_page(self):
        c = Client()
        with patch('designer.views.start_analysis', return_value=True) as mock_start:
            r = c.post(
                '/loading/',
                {
                    'uniprot_id': 'P12345',
                    'editor': 'ABE',
                    'alpha_threshold': '0.5',
                    'top_sgrnas': '5',
                    'window_min': '4',
                    'window_max': '8',
                    'pam_type': 'NGG',
                },
            )
        self.assertEqual(r.status_code, 200)
        self.assertContains(r, 'Running Analysis')
        self.assertContains(r, 'analysis-status')
        mock_start.assert_called_once()

    def test_analysis_status_returns_json(self):
        c = Client()
        r = c.get('/analysis-status/')
        self.assertEqual(r.status_code, 200)
        self.assertEqual(r.json()['status'], 'idle')

    def test_cadd_status_returns_json(self):
        c = Client()
        r = c.get('/results/cadd-status/')
        self.assertEqual(r.status_code, 200)
        self.assertIn(r.json()['status'], ('idle', 'pending', 'done', 'error', 'running', 'skipped'))

    def test_results_get_without_data_redirects_home(self):
        c = Client()
        r = c.get('/results/')
        self.assertEqual(r.status_code, 302)
        self.assertEqual(r.url, '/')
