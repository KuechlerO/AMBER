"""Tests for SAFFRON background analysis jobs."""

import time
from unittest.mock import patch

from django.core.cache import cache
from django.test import Client, TestCase

from designer.services import UserInputError

from saffron.analysis_jobs import clear_job_state, get_status, start_analysis
from saffron.result_store import analysis_cache_key, load_analysis_results_by_key


class AnalysisJobsTests(TestCase):
    session_key = 'test-saffron-session'

    def setUp(self):
        clear_job_state(self.session_key)
        cache.delete(analysis_cache_key(self.session_key))

    def tearDown(self):
        clear_job_state(self.session_key)
        cache.delete(analysis_cache_key(self.session_key))

    def test_idle_when_no_job(self):
        self.assertEqual(get_status(self.session_key)['status'], 'idle')

    @patch('saffron.analysis_jobs.run_saffron_analysis')
    def test_job_completes_and_saves_results(self, mock_run):
        mock_run.return_value = {'guide_rows': [], 'mode': 'sequence', 'wt_signalp': {}}
        self.assertTrue(start_analysis(self.session_key, {'input_mode': 'sequence', 'aa_sequence': 'M' * 30}))
        self.assertEqual(get_status(self.session_key)['status'], 'running')

        for _ in range(100):
            if get_status(self.session_key)['status'] == 'done':
                break
            time.sleep(0.05)

        self.assertEqual(get_status(self.session_key)['status'], 'done')
        self.assertTrue(load_analysis_results_by_key(self.session_key))

    @patch('saffron.analysis_jobs.run_saffron_analysis')
    def test_job_user_error(self, mock_run):
        mock_run.side_effect = UserInputError('Invalid input')
        start_analysis(self.session_key, {'input_mode': 'sequence'})

        for _ in range(100):
            status = get_status(self.session_key)
            if status['status'] == 'error':
                break
            time.sleep(0.05)

        status = get_status(self.session_key)
        self.assertEqual(status['status'], 'error')
        self.assertIn('Invalid input', status['error'])

    @patch('saffron.analysis_jobs.run_saffron_analysis')
    def test_second_start_while_running_is_noop(self, mock_run):
        def slow_run(_form):
            time.sleep(0.3)
            return {'guide_rows': [], 'mode': 'sequence'}

        mock_run.side_effect = slow_run
        self.assertTrue(start_analysis(self.session_key, {'input_mode': 'sequence'}))
        self.assertFalse(start_analysis(self.session_key, {'input_mode': 'sequence'}))


class SaffronAsyncViewTests(TestCase):
    def test_loading_post_returns_polling_page(self):
        c = Client()
        with patch('saffron.views.start_analysis', return_value=True) as mock_start:
            r = c.post(
                '/saffron/loading/',
                {
                    'input_mode': 'sequence',
                    'aa_sequence': 'ACDEFGHIKLM' * 3,
                },
            )
        self.assertEqual(r.status_code, 200)
        self.assertContains(r, 'Running SAFFRON')
        self.assertContains(r, 'analysis-status')
        mock_start.assert_called_once()

    def test_analysis_status_returns_json(self):
        c = Client()
        r = c.get('/saffron/analysis-status/')
        self.assertEqual(r.status_code, 200)
        self.assertEqual(r.json()['status'], 'idle')

    def test_results_get_without_data_redirects_home(self):
        c = Client()
        r = c.get('/saffron/results/')
        self.assertEqual(r.status_code, 302)
        self.assertIn('/saffron/', r.url)
