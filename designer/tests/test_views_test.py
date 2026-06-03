from django.test import TestCase
from django.urls import reverse

from unittest.mock import patch
from designer.services import UserInputError
from designer.views import filter_rows_by_score, resolve_score_filter

import math

class SimpleViewsTest(TestCase):

    # tests Template
    def test_home_view(self):
        response = self.client.get(reverse('home'))
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, 'designer/home.html')

    def test_loading_view(self):
        response = self.client.get(reverse('loading'))
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, 'designer/loading.html')

    # tests results
    # test get -> redirect
    def test_results_get_redirects(self):
        response = self.client.get(reverse('results'))
        self.assertEqual(response.status_code, 302)


    # simulate run_analysis
    @patch('designer.views.run_analysis')
    def test_results_post_success(self, mock_run_analysis):
        mock_run_analysis.return_value = [
            {'position': 1, 'alpha_score': 0.8}
        ]

        response = self.client.post(reverse('results'), {
            'uniprot_id': 'P04439',
            'editor': 'ABE',
            'alpha_threshold': '0.5',
            'top_sgrnas': '5',
            'window_min': '4',
            'window_max': '8',
        })

        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, 'designer/results.html')

        from designer.result_store import load_analysis_results
        cached = load_analysis_results(response.wsgi_request)
        self.assertIn('filtered_rows', cached)
        self.assertIn('full_rows', cached)

    # tests if only_over_threshold works
    @patch('designer.views.run_analysis')
    def test_threshold_filtering(self, mock_run_analysis):
        mock_run_analysis.return_value = [
            {'position': 1, 'alpha_score': 0.9},
            {'position': 2, 'alpha_score': 0.2},
        ]

        response = self.client.post(reverse('results'), {
            'uniprot_id': 'P04439',
            'editor': 'ABE',
            'alpha_threshold': '0.5',
            'output': 'only_over_threshold',
        })

        results = response.context['results']
        self.assertEqual(len(results), 1)

    # User Input error test
    @patch('designer.views.run_analysis')
    def test_user_input_error(self, mock_run_analysis):
        mock_run_analysis.side_effect = UserInputError('Invalid input')

        response = self.client.post(reverse('results'), {
            'uniprot_id': 'BAD',
        })

        self.assertTemplateUsed(response, 'designer/error.html')
        self.assertEqual(response.context['error_type'], 'user')

    # invlaid editing window
    @patch('designer.views.run_analysis')
    def test_invalid_window_values_become_none(self, mock_run_analysis):
        mock_run_analysis.return_value = []

        self.client.post(reverse('results'), {
            'uniprot_id': 'P04439',
            'editor': 'ABE',
            'window_min': 'abc',
            'window_max': 'xyz',
        })

        mock_run_analysis.assert_called_once()
        kwargs = mock_run_analysis.call_args.kwargs

        self.assertIsNone(kwargs['window_min'])
        self.assertIsNone(kwargs['window_max'])



    # regular exception test
    @patch('designer.views.run_analysis')
    def test_server_error(self, mock_run_analysis):
        mock_run_analysis.side_effect = Exception('Crash')

        response = self.client.post(reverse('results'), {
            'uniprot_id': 'P04439',
        })

        self.assertTemplateUsed(response, 'designer/error.html')
        self.assertEqual(response.context['error_type'], 'server')

    # download excel test
    @patch('designer.views.run_analysis')
    def test_download_excel(self, mock_run_analysis):
        mock_run_analysis.return_value = [{'position': 1, 'alpha_score': 0.8}]
        self.client.post(reverse('results'), {
            'uniprot_id': 'P04439',
            'editor': 'ABE',
        })
        response = self.client.get(reverse('download_excel'))

        self.assertEqual(response.status_code, 200)
        self.assertEqual(
            response['Content-Type'],
            'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
        )

    @patch('designer.views.run_analysis')
    def test_download_excel_has_content(self, mock_run_analysis):
        mock_run_analysis.return_value = [{'position': 1, 'alpha_score': 0.8}]
        self.client.post(reverse('results'), {
            'uniprot_id': 'P04439',
            'editor': 'ABE',
        })
        response = self.client.get(reverse('download_excel'))

        self.assertGreater(len(response.content), 0)

    # download csv test
    @patch('designer.views.run_analysis')
    def test_download_csv(self, mock_run_analysis):
        mock_run_analysis.return_value = [{'position': 1, 'alpha_score': 0.8}]
        self.client.post(reverse('results'), {
            'uniprot_id': 'P04439',
            'editor': 'ABE',
        })
        response = self.client.get(reverse('download_csv'))

        self.assertEqual(response.status_code, 200)
        self.assertEqual(response['Content-Type'], 'text/csv')

    @patch('designer.views.run_analysis')
    def test_download_csv_content(self, mock_run_analysis):
        mock_run_analysis.return_value = [{'position': 1, 'alpha_score': 0.8}]
        self.client.post(reverse('results'), {
            'uniprot_id': 'P04439',
            'editor': 'ABE',
            'columns': ['position'],
        })
        response = self.client.get(reverse('download_csv'))

        content = response.content.decode()
        self.assertIn('Position', content)
        self.assertIn('1', content)

# tests results
class ResultsCacheTest(TestCase):

    @patch('designer.views.run_analysis')
    def test_cache_is_saved(self, mock_run_analysis):
        mock_run_analysis.return_value = [
            {'position': 1, 'alpha_score': 0.8}
        ]

        response = self.client.post(reverse('results'), {
            'uniprot_id': 'P04439',
            'editor': 'ABE',
        })

        from designer.result_store import load_analysis_results
        cached = load_analysis_results(response.wsgi_request)
        self.assertIn('filtered_rows', cached)
        self.assertIn('full_rows', cached)
        self.assertIn('form_data', cached)
        self.assertEqual(cached['form_data']['uniprot_id'], 'P04439')

    # tests if filtering with math works
    @patch('designer.views.run_analysis')
    def test_nan_and_none_filtering(self, mock_run_analysis):
        mock_run_analysis.return_value = [
            {'position': 1, 'alpha_score': 0.9},
            {'position': 2, 'alpha_score': None},
            {'position': 3, 'alpha_score': float('nan')},
            {'position': 4, 'alpha_score': 0.3},
        ]

        response = self.client.post(reverse('results'), {
            'uniprot_id': 'P04439',
            'editor': 'ABE',
            'alpha_threshold': '0.5',
            'output': 'only_over_threshold',
        })

        results = response.context['results']

        # only values >= 0.5
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0]['position'], 1)

    # test position_count, guide_count
    @patch('designer.views.run_analysis')
    def test_counts_in_context(self, mock_run_analysis):
        mock_run_analysis.return_value = [
            {'position': 1, 'alpha_score': 0.9},
            {'position': 1, 'alpha_score': 0.8},
            {'position': 2, 'alpha_score': 0.7},
        ]

        response = self.client.post(reverse('results'), {
            'uniprot_id': 'P04439',
            'editor': 'ABE',
        })

        self.assertEqual(response.context['position_count'], 2)
        self.assertEqual(response.context['guide_count'], 3)

    # test editor display
    @patch('designer.views.run_analysis')
    def test_editor_display_both(self, mock_run_analysis):
        mock_run_analysis.return_value = []

        response = self.client.post(reverse('results'), {
            'uniprot_id': 'P04439',
            'editor': 'BOTH',
        })

        self.assertEqual(
            response.context['form_data']['editor_display'],
            'ABE & CBE'
        )

    # tests selected_columns
    @patch('designer.views.run_analysis')
    def test_default_columns_used_when_empty(self, mock_run_analysis):
        mock_run_analysis.return_value = [
            {'position': 1, 'alpha_score': 0.8}
        ]

        response = self.client.post(reverse('results'), {
            'uniprot_id': 'P04439',
            'editor': 'ABE',
            'columns': []
        })

        selected_columns = response.context['form_data']['selected_columns']

        self.assertIn('position', selected_columns)
        self.assertIn('alpha_score', selected_columns)

    # tests if editor_used is selected correct
    @patch('designer.views.run_analysis')
    def test_editor_both_adds_editor_used_column(self, mock_run_analysis):
        mock_run_analysis.return_value = [
            {'position': 1, 'alpha_score': 0.8}
        ]

        response = self.client.post(reverse('results'), {
            'uniprot_id': 'P04439',
            'editor': 'BOTH',
            'columns': ['position', 'alpha_score']
        })

        selected_columns = response.context['form_data']['selected_columns']

        self.assertIn('editor_used', selected_columns)

        self.assertEqual(
            selected_columns,
            ['position', 'alpha_score', 'editor_used']
        )


class ResultsRefreshTest(TestCase):

    @patch('designer.views.run_analysis')
    def test_refresh_only_skips_run_analysis(self, mock_run_analysis):
        mock_run_analysis.return_value = [
            {'position': 1, 'alpha_score': 0.9},
            {'position': 2, 'alpha_score': 0.2},
        ]
        self.client.post(reverse('results'), {
            'uniprot_id': 'P04439',
            'editor': 'ABE',
            'alpha_threshold': '0.5',
        })
        mock_run_analysis.reset_mock()

        response = self.client.post(reverse('results'), {
            'refresh_only': '1',
            'uniprot_id': 'P04439',
            'editor': 'ABE',
            'alpha_threshold': '0.5',
            'filter_threshold': 'true',
            'columns': ['position', 'alpha_score'],
        })

        mock_run_analysis.assert_not_called()
        self.assertEqual(response.status_code, 200)
        self.assertEqual(len(response.context['results']), 1)
        self.assertEqual(
            response.context['form_data']['selected_columns'],
            ['position', 'alpha_score'],
        )

    def test_refresh_without_cache_redirects_home(self):
        response = self.client.post(reverse('results'), {
            'refresh_only': '1',
            'uniprot_id': 'P04439',
            'editor': 'ABE',
        })
        self.assertEqual(response.status_code, 302)

    @patch('designer.views.run_analysis')
    def test_refresh_below_cutoff_filters_rows(self, mock_run_analysis):
        mock_run_analysis.return_value = [
            {'position': 1, 'alpha_score': 0.9},
            {'position': 2, 'alpha_score': 0.2},
            {'position': 3, 'alpha_score': 0.4},
        ]
        self.client.post(reverse('results'), {
            'uniprot_id': 'P04439',
            'editor': 'ABE',
            'alpha_threshold': '0.5',
            'score_filter_mode': 'all',
        })

        response = self.client.post(reverse('results'), {
            'refresh_only': '1',
            'uniprot_id': 'P04439',
            'editor': 'ABE',
            'alpha_threshold': '0.5',
            'score_filter_mode': 'below',
            'display_score_cutoff': '0.5',
            'columns': ['position', 'alpha_score'],
        })

        self.assertEqual(response.status_code, 200)
        scores = [r['alpha_score'] for r in response.context['results']]
        self.assertEqual(scores, [0.2, 0.4])
        self.assertEqual(response.context['form_data']['score_filter_mode'], 'below')


class ScoreFilterTest(TestCase):

    def test_filter_rows_below_cutoff(self):
        rows = [
            {'alpha_score': 0.9},
            {'alpha_score': 0.2},
            {'alpha_score': None},
        ]
        filtered = filter_rows_by_score(rows, 'below', 0.5)
        self.assertEqual([r['alpha_score'] for r in filtered], [0.2])

    def test_resolve_score_filter_legacy_and_new_fields(self):
        mode, cutoff = resolve_score_filter(
            {'filter_threshold': 'true', 'alpha_threshold': '0.6'},
            default_cutoff='0.5',
        )
        self.assertEqual((mode, cutoff), ('above', 0.6))

        mode, cutoff = resolve_score_filter(
            {
                'score_filter_mode': 'below',
                'display_score_cutoff': '0.35',
            },
            default_cutoff='0.5',
        )
        self.assertEqual((mode, cutoff), ('below', 0.35))
