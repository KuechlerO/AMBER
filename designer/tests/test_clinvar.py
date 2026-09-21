"""Tests for ClinVar E-utilities client and residue parsing."""

from __future__ import annotations

import tempfile
from unittest.mock import patch

from django.test import SimpleTestCase, override_settings

from designer.clinvar import (
    attach_clinvar_to_rows,
    fetch_clinvar_protein_variants,
    format_clinvar_export,
    index_clinvar_by_residue,
    parse_esummary_record,
    parse_residue_from_protein_change,
    parse_residue_from_title,
)


class ClinvarParseTest(SimpleTestCase):
    def test_one_letter_missense(self):
        self.assertEqual(parse_residue_from_protein_change('L260V'), 260)

    def test_nonsense_star(self):
        self.assertEqual(parse_residue_from_protein_change('R514*'), 514)

    def test_multi_isoform_respects_length(self):
        # Prefer largest residue within protein_length.
        self.assertEqual(
            parse_residue_from_protein_change('H139fs, H19fs, H178fs, H46fs', protein_length=200),
            178,
        )
        self.assertEqual(
            parse_residue_from_protein_change('H139fs, H19fs, H178fs, H46fs', protein_length=100),
            46,
        )

    def test_title_fallback(self):
        title = 'NM_000206.3(IL2RG):c.778T>G (p.Leu260Val)'
        self.assertEqual(parse_residue_from_title(title), 260)
        self.assertIsNone(parse_residue_from_title(title, protein_length=100))

    def test_esummary_record_maps_fields(self):
        record = {
            'uid': '4839707',
            'accession': 'VCV004839707',
            'title': 'NM_000206.3(IL2RG):c.778T>G (p.Leu260Val)',
            'protein_change': 'L260V',
            'germline_classification': {
                'description': 'Uncertain significance',
                'review_status': 'criteria provided, single submitter',
            },
            'molecular_consequence_list': ['missense variant'],
        }
        parsed = parse_esummary_record(record, protein_length=369)
        self.assertIsNotNone(parsed)
        self.assertEqual(parsed['residue'], 260)
        self.assertEqual(parsed['classification'], 'Uncertain significance')
        self.assertEqual(parsed['uid'], '4839707')
        self.assertIn('clinvar/variation/4839707', parsed['url'])

    def test_esummary_drops_unmapped(self):
        record = {
            'uid': '1',
            'title': 'NM_000206.3(IL2RG):c.1-10A>G',
            'protein_change': '',
            'germline_classification': {'description': 'Pathogenic'},
        }
        self.assertIsNone(parse_esummary_record(record, protein_length=369))


class ClinvarFetchTest(SimpleTestCase):
    def test_empty_gene_soft_fails(self):
        variants, err = fetch_clinvar_protein_variants('', protein_length=100)
        self.assertEqual(variants, [])
        self.assertIn('gene symbol', err or '')

    def test_fetch_parses_and_caches(self):
        cache_dir = tempfile.mkdtemp(prefix='clinvar_test_')
        esearch = {
            'esearchresult': {'count': '1', 'idlist': ['4839707']},
        }
        esummary = {
            'result': {
                'uids': ['4839707'],
                '4839707': {
                    'uid': '4839707',
                    'accession': 'VCV004839707',
                    'title': 'NM_000206.3(IL2RG):c.778T>G (p.Leu260Val)',
                    'protein_change': 'L260V',
                    'germline_classification': {
                        'description': 'Uncertain significance',
                        'review_status': 'criteria provided, single submitter',
                    },
                    'molecular_consequence_list': ['missense variant'],
                },
            }
        }

        with override_settings(CLINVAR_CACHE_DIR=cache_dir):
            with patch('designer.clinvar._http_get', side_effect=[esearch, esummary]) as mock_get:
                variants, err = fetch_clinvar_protein_variants('IL2RG', protein_length=369)
                self.assertIsNone(err)
                self.assertEqual(len(variants), 1)
                self.assertEqual(variants[0]['residue'], 260)
                self.assertEqual(mock_get.call_count, 2)

            # Second call should hit cache (no HTTP).
            with patch('designer.clinvar._http_get') as mock_get2:
                variants2, err2 = fetch_clinvar_protein_variants('IL2RG', protein_length=369)
                self.assertIsNone(err2)
                self.assertEqual(len(variants2), 1)
                mock_get2.assert_not_called()

    def test_network_error_soft_fails(self):
        with patch('designer.clinvar._http_get', side_effect=RuntimeError('down')):
            with patch('designer.clinvar._read_cache', return_value=None):
                variants, err = fetch_clinvar_protein_variants('IL2RG', protein_length=369, use_cache=False)
        self.assertEqual(variants, [])
        self.assertIn('ClinVar unavailable', err or '')


class ClinvarAttachTest(SimpleTestCase):
    def test_index_and_attach_by_position(self):
        variants = [
            {
                'uid': '1',
                'accession': 'VCV1',
                'protein_change': 'L10V',
                'residue': 10,
                'classification': 'Pathogenic',
                'url': 'https://www.ncbi.nlm.nih.gov/clinvar/variation/1/',
            },
            {
                'uid': '2',
                'accession': 'VCV2',
                'protein_change': 'R10C',
                'residue': 10,
                'classification': 'Uncertain significance',
                'url': 'https://www.ncbi.nlm.nih.gov/clinvar/variation/2/',
            },
            {
                'uid': '3',
                'accession': 'VCV3',
                'protein_change': 'G20A',
                'residue': 20,
                'classification': 'Benign',
                'url': 'https://www.ncbi.nlm.nih.gov/clinvar/variation/3/',
            },
        ]
        by_res = index_clinvar_by_residue(variants)
        self.assertEqual(len(by_res[10]), 2)
        self.assertEqual(by_res[10][0]['css_class'], 'clinvar-pathogenic')

        rows = attach_clinvar_to_rows(
            [{'position': 10, 'sgrna_seq': 'A'}, {'position': 99, 'sgrna_seq': 'B'}],
            by_res,
        )
        self.assertEqual(len(rows[0]['clinvar_at_position']), 2)
        self.assertIn('L10V|Pathogenic|VCV1', rows[0]['clinvar'])
        self.assertEqual(rows[1]['clinvar_at_position'], [])
        self.assertEqual(rows[1]['clinvar'], '')

    def test_format_export(self):
        text = format_clinvar_export([
            {'protein_change': 'L10V', 'classification': 'Pathogenic', 'accession': 'VCV1'},
        ])
        self.assertEqual(text, 'L10V|Pathogenic|VCV1')
