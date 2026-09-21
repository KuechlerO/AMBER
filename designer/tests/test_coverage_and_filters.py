"""Tests for ID normalization, score filters, and coverage stats."""

from unittest.mock import patch

from django.test import SimpleTestCase

from designer.coverage_stats import classify_positions, compute_coverage_stats
from designer.pipeline import extract_ensembl_transcript
from designer.services import looks_like_uniprot_accession, normalize_input_id
from designer.views import filter_rows_by_score


class NormalizeInputIdTest(SimpleTestCase):
    def test_strips_and_removes_internal_spaces(self):
        self.assertEqual(normalize_input_id('  P04439  '), 'P04439')
        self.assertEqual(normalize_input_id('P0 44 39'), 'P04439')
        self.assertEqual(normalize_input_id(' ENST000001 '), 'ENST000001')

    def test_uppercases_ids(self):
        self.assertEqual(normalize_input_id('p04439'), 'P04439')
        self.assertEqual(normalize_input_id('a0a5b9'), 'A0A5B9')


class UniprotAccessionSyntaxTest(SimpleTestCase):
    def test_accepts_swissprot_and_trembl_forms(self):
        self.assertTrue(looks_like_uniprot_accession('P04439'))
        self.assertTrue(looks_like_uniprot_accession('P01850'))
        self.assertTrue(looks_like_uniprot_accession('A0A5B9'))
        self.assertTrue(looks_like_uniprot_accession('A0A075B6S1'))
        self.assertTrue(looks_like_uniprot_accession('P60709-1'))
        self.assertTrue(looks_like_uniprot_accession('q9y5w3'))

    def test_rejects_non_accessions(self):
        self.assertFalse(looks_like_uniprot_accession('TRBC1'))
        self.assertFalse(looks_like_uniprot_accession('ENST00000376809'))
        self.assertFalse(looks_like_uniprot_accession('ABC'))
        self.assertFalse(looks_like_uniprot_accession(''))


class EnsemblTranscriptResolutionTest(SimpleTestCase):
    def test_uses_direct_ensembl_xref(self):
        data = {
            'uniProtKBCrossReferences': [
                {'database': 'Ensembl', 'id': 'ENST00000376809.9'},
            ],
            'genes': [{'geneName': {'value': 'HLA-A'}}],
        }
        self.assertEqual(extract_ensembl_transcript(data), 'ENST00000376809')

    def test_falls_back_to_gene_symbol_lookup(self):
        data = {
            'uniProtKBCrossReferences': [],
            'genes': [{'geneName': {'value': 'TRBC1'}}],
        }
        with patch(
            'designer.pipeline.lookup_ensembl_canonical_transcript',
            return_value='ENST00000633705',
        ) as mock_lookup:
            tid = extract_ensembl_transcript(data)
        self.assertEqual(tid, 'ENST00000633705')
        mock_lookup.assert_called_once_with(gene_symbol='TRBC1')

    def test_falls_back_to_ensg_lookup(self):
        data = {
            'uniProtKBCrossReferences': [
                {'database': 'HPA', 'id': 'ENSG00000211751'},
            ],
            'genes': [],
        }
        with patch(
            'designer.pipeline.lookup_ensembl_canonical_transcript',
            return_value='ENST00000633705',
        ) as mock_lookup:
            tid = extract_ensembl_transcript(data)
        self.assertEqual(tid, 'ENST00000633705')
        mock_lookup.assert_called_once_with(gene_id='ENSG00000211751')

    def test_raises_when_unresolvable(self):
        with self.assertRaises(ValueError):
            extract_ensembl_transcript({'uniProtKBCrossReferences': [], 'genes': []})


class ScoreFilterTest(SimpleTestCase):
    def setUp(self):
        self.rows = [
            {'alpha_score': 0.5, 'position': 1},
            {'alpha_score': 0.6, 'position': 2},
            {'alpha_score': 0.7, 'position': 3},
        ]

    def test_above_is_inclusive(self):
        out = filter_rows_by_score(self.rows, 'above', 0.6)
        self.assertEqual([r['position'] for r in out], [2, 3])

    def test_below_is_inclusive(self):
        out = filter_rows_by_score(self.rows, 'below', 0.6)
        self.assertEqual([r['position'] for r in out], [1, 2])

    def test_all_keeps_everything(self):
        out = filter_rows_by_score(self.rows, 'all', 0.6)
        self.assertEqual(len(out), 3)


class CoverageStatsTest(SimpleTestCase):
    def test_editor_and_score_segmentation(self):
        rows = [
            {'position': 10, 'editor_used': 'ABE', 'alpha_score': 0.9, 'sgrna_seq': 'A'},
            {'position': 10, 'editor_used': 'CBE', 'alpha_score': 0.4, 'sgrna_seq': 'B'},
            {'position': 20, 'editor_used': 'ABE', 'alpha_score': 0.3, 'sgrna_seq': 'C'},
            {'position': 30, 'editor_used': 'CBE', 'alpha_score': 0.8, 'sgrna_seq': 'D'},
        ]
        classified = classify_positions(rows, 0.6)
        self.assertEqual(classified[10], ('both', 'high'))  # max score 0.9
        self.assertEqual(classified[20], ('ABE', 'low'))
        self.assertEqual(classified[30], ('CBE', 'high'))

        stats = compute_coverage_stats(rows, protein_length=100, cutoff=0.6)
        self.assertEqual(stats['editable_positions'], 3)
        self.assertEqual(stats['uncovered_positions'], 97)
        self.assertEqual(stats['coverage_pct'], 3.0)
        self.assertEqual(stats['segment_counts'][('both', 'high')], 1)
        self.assertEqual(stats['segment_counts'][('ABE', 'low')], 1)
        self.assertEqual(stats['segment_counts'][('CBE', 'high')], 1)
        self.assertEqual(stats['segment_counts'][('uncovered', None)], 97)
        self.assertAlmostEqual(stats['segment_pct'][('both', 'high')], 1.0)
        self.assertEqual(stats['below_cutoff_positions'], 1)
        self.assertEqual(stats['above_cutoff_positions'], 2)
