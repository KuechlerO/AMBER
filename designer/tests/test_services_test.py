from django.test import TestCase
from unittest.mock import patch
from designer.services import run_analysis
from designer.services import UserInputError
import requests

class RunAnalysisTests(TestCase):

    # tests services.py
    @patch('designer.services.run_pipeline')
    def test_run_analysis_success(self, mock_pipeline):
        mock_pipeline.return_value = [{'position': 1, 'alpha_score': 0.8}]

        result = run_analysis('P04439', 'ABE', '0.5', '5', '1', '10')
        self.assertIsInstance(result, list)
        self.assertEqual(result[0]['position'], 1)

    # tests invalid uniprot_id
    def test_empty_uniprot_id(self):
        with self.assertRaises(UserInputError):
            run_analysis('', 'ABE', '0.5', '5', '1', '10')

    @patch('designer.services.is_valid_uniprot_id')
    @patch('designer.services.run_pipeline')
    def test_run_analysis_valid_uniprot(self, mock_pipeline, mock_valid_id):
        mock_valid_id.return_value = True
        mock_pipeline.return_value = [{'position': 1, 'alpha_score': 0.8}]

        result = run_analysis('P12345', 'ABE', '0.5', '5', '1', '10')
        self.assertIsInstance(result, list)
        self.assertEqual(result[0]['position'], 1)

    @patch('designer.services.is_valid_uniprot_id')
    def test_run_analysis_invalid_uniprot(self, mock_valid_id):
        mock_valid_id.return_value = False

        with self.assertRaises(UserInputError) as e:
            run_analysis('INVALID', 'ABE', '0.5', '5', '1', '10')
        self.assertIn('Invalid UniProt ID/ Ensembl provided.', str(e.exception))

    # test ensembl
    @patch('designer.services.is_valid_uniprot_id')
    @patch('designer.services.is_valid_ensembl_id')
    @patch('designer.services.run_pipeline')
    def test_valid_ensembl_id(self, mock_pipeline, mock_ensembl, mock_uniprot):
        mock_uniprot.return_value = False
        mock_ensembl.return_value = True
        mock_pipeline.return_value = []

        result = run_analysis('ENST000003', 'ABE', '0.5', '5', '1', '10')
        self.assertEqual(result, [])

    # tests invalid editor input
    def test_invalid_editor(self):
        with self.assertRaises(UserInputError):
            run_analysis('P04439', 'XYZ', '0.5', '5', '1', '10')

    # tests invalid alpha_threshold input
    def test_invalid_alpha_threshold(self):
        with self.assertRaises(UserInputError):
            run_analysis('P04439', 'ABE', 'abc', '5', '1', '10')

    # tests invalid top_sgRNA input
    def test_invalid_top_sgrnas(self):
        with self.assertRaises(UserInputError):
            run_analysis('P04439', 'ABE', '0.5', 'xyz', '1', '10')

    # tests invalid editing window
    @patch('designer.services.run_pipeline')
    @patch('designer.services.is_valid_uniprot_id')
    def test_valid_window_range(self, mock_valid_id, mock_pipeline):
        mock_valid_id.return_value = True
        mock_pipeline.return_value = []

        result = run_analysis('P12345', 'ABE', '0.5', '5', '3', '10')
        self.assertEqual(result, [])

    def test_invalid_window_values_type(self):
        with self.assertRaises(UserInputError):
            run_analysis('P04439', 'ABE', '0.5', '5', 'a', '10')

    def test_window_out_of_range(self):
        with self.assertRaises(UserInputError):
            run_analysis('P04439', 'ABE', '0.5', '5', '0', '25')

    def test_window_min_greater_than_max(self):
        with self.assertRaises(UserInputError):
            run_analysis('P04439', 'ABE', '0.5', '5', '10', '5')



    # tests exception in pipeline
    @patch('designer.services.run_pipeline')
    def test_pipeline_file_not_found(self, mock_pipeline):
        mock_pipeline.side_effect = FileNotFoundError('file missing')

        with self.assertRaises(UserInputError) as e:
            run_analysis('P04439', 'ABE', '0.5', '5', '1', '10')
        self.assertIn('file missing', str(e.exception))

    @patch('designer.services.is_valid_uniprot_id')
    @patch('designer.services.run_pipeline')
    def test_requests_exception(self, mock_pipeline, mock_valid_id):
        mock_valid_id.return_value = True
        mock_pipeline.side_effect = requests.RequestException('API down')

        with self.assertRaises(UserInputError) as e:
            run_analysis('P12345', 'ABE', '0.5', '5', '1', '10')

        self.assertIn('Could not retrieve sequence data', str(e.exception))

    # test pipeline
    @patch('designer.services.run_pipeline')
    @patch('designer.services.is_valid_uniprot_id')
    def test_pipeline_called_with_correct_args(self, mock_valid_id, mock_pipeline):
        mock_valid_id.return_value = True

        run_analysis('P12345', 'ABE', '0.5', '5', '2', '8', 'hide')

        mock_pipeline.assert_called_once_with(
            uniprot_id='P12345',
            editor='ABE',
            alpha_threshold='0.5',
            top_sgrnas='5',
            window_min=2,
            window_max=8,
            duplicate_mode='hide',
        )
