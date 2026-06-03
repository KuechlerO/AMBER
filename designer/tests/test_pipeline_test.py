from django.test import TestCase
import unittest
from unittest.mock import patch, MagicMock
from django.test import TestCase
import pandas as pd


class PipelineTests(TestCase):

    # tests simple help functions
    def test_complementary_base(self):
        from designer.pipeline import complementary_base

        self.assertEqual(complementary_base('A'), 'G')
        self.assertEqual(complementary_base('T'), 'C')
        self.assertEqual(complementary_base('G'), 'A')
        self.assertEqual(complementary_base('C'), 'T')

    def test_translate_codon(self):
        from designer.pipeline import translate_codon

        self.assertEqual(translate_codon('ATG'), 'M')  # Startcodon
        self.assertEqual(translate_codon('XXX'), 'X')  # Error


    # tests flatten_candidate_dataframe
    def test_flatten_candidate_dataframe(self):
        from designer.pipeline import flatten_candidate_dataframe

        df = pd.DataFrame([{
            'position': 1,
            'WT_codon': 'ATG',
            'WT_AA': 'M',
            'ABE_mutated_AA': 'V',
            'ABE_mutated_AA_score': 0.9,
            'sgRNA_1_seq': 'AAA',
            'sgRNA_1_protospacer': 'AAA',
            'sgRNA_1_pam': 'GGG',
            'sgRNA_1_strand': '+',
            'sgRNA_1_TargetApos': 5,
            'sgRNA_1_numAsWindow': 2,
            'sgRNA_1_outcomes': ['AAA(V)']
        }])

        result = flatten_candidate_dataframe(df, 'ABE', 5)

        self.assertEqual(len(result), 1)
        self.assertEqual(result[0]['position'], 1)
        self.assertEqual(result[0]['sgrna_seq'], 'AAA')


    # tests API
    @patch('designer.pipeline.requests.get')
    def test_get_cds_and_transcript_from_uniprot(self, mock_get):
        from designer.pipeline import get_cds_and_transcript_from_uniprot

        # Mock für UniProt JSON
        mock_uniprot_response = MagicMock()
        mock_uniprot_response.json.return_value = {
            'uniProtKBCrossReferences': [
                {'database': 'Ensembl', 'id': 'ENST00000335137.4'}
            ]
        }
        mock_uniprot_response.raise_for_status = MagicMock()

        # Mock for Ensembl CDS
        mock_ensembl_response = MagicMock()
        mock_ensembl_response.status_code = 200
        mock_ensembl_response.text = 'ATGGCC'


        mock_get.side_effect = [mock_uniprot_response, mock_ensembl_response]

        cds, transcript = get_cds_and_transcript_from_uniprot('P12345')

        self.assertEqual(transcript, 'ENST00000335137')
        self.assertEqual(cds, 'ATGGCC')

    # tests DB
    @patch('designer.pipeline.Alpha_missense')
    def test_get_alphamissense_from_db(self, mock_model):
        from designer.pipeline import get_alphamissense_from_db

        mock_model.objects.using.return_value.filter.return_value.values_list.return_value = [
            ('A123B', 0.8),
            ('C456D', 0.2)
        ]

        df = get_alphamissense_from_db('P12345', 0.5)

        self.assertIn('position', df.columns)
        self.assertEqual(len(df), 2)

    # Pipeline Test
    @patch('designer.pipeline.get_cds_and_transcript_from_uniprot')
    @patch('designer.pipeline.get_alphamissense_from_db')
    @patch('designer.pipeline.generate_candidates')
    @patch('designer.pipeline.flatten_candidate_dataframe')
    def test_run_pipeline(
        self,
        mock_flatten,
        mock_generate,
        mock_db,
        mock_cds
    ):
        from designer.pipeline import run_pipeline

        # --- Mock CDS fetch ---
        mock_cds.return_value = ('ATGGCC', 'ENST0001')

        # --- Mock DB DataFrame (must satisfy kill_doubles) ---
        mock_db.return_value = pd.DataFrame({
            'position': [1],
            'WT_codon': ['ATG'],
            'WT_AA': ['M'],
            'avg_alpha_score': [0.9],

            'ABE_mutated_AA': ['I'],
            'ABE_mutated_AA_score': [0.9],

            # Critical for kill_doubles
            'sgRNA_1_seq': ['AAA'],
        })

        # --- Mock candidate generation (must also be DataFrame-like) ---
        mock_generate.return_value = pd.DataFrame({
            'position': [1],
            'sgRNA_1_seq': ['AAA'],
        })

        # --- Mock flattening (final output we assert on) ---
        mock_flatten.return_value = [
            {
                'position': 1,
                'alpha_score': 0.9,
                'sgrna_seq': 'AAA',
                'pam': 'GGG',
                'strand': '+'
            }
        ]

        # --- Run pipeline ---
        result = run_pipeline(
            'P12345',
            'ABE',
            '0.5',
            '5',
            1,
            10,
            pam_type='NGG',
        )

        # --- Assertions ---
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0]['position'], 1)
        self.assertEqual(result[0]['sgrna_seq'], 'AAA')

        # --- Optional: verify orchestration ---
        mock_cds.assert_called_once_with('P12345')
        mock_db.assert_called_once()
        mock_generate.assert_called_once()
        mock_flatten.assert_called_once()
