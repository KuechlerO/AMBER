"""ESM1b lookup, CDS→genome mapping, and mocked CADD client."""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import Mock

import pandas as pd
from django.test import SimpleTestCase, override_settings

from designer.cadd import (
    fetch_cadd_scores,
    parse_cadd_payload,
)
from designer.esm1b import clear_esm1b_cache, lookup_esm1b
from designer.pipeline import (
    cds_index_to_genomic,
    coding_edit_to_genomic_snv,
    finalize_guide_scores,
    flatten_candidate_dataframe,
    generate_candidates,
    parse_cds_mappings,
    wc_complement,
)

FIXTURE_DIR = Path(__file__).resolve().parent / 'fixtures' / 'esm1b'


class Esm1bLookupTest(SimpleTestCase):
    def setUp(self):
        clear_esm1b_cache()

    def tearDown(self):
        clear_esm1b_cache()

    @override_settings(ESM1B_DATA_DIR=str(FIXTURE_DIR))
    def test_lookup_wt_match(self):
        self.assertAlmostEqual(lookup_esm1b('P99999', 2, 'A', 'K'), -11.552)
        self.assertAlmostEqual(lookup_esm1b('P99999', 2, 'K', 'K'), 0.0)

    @override_settings(ESM1B_DATA_DIR=str(FIXTURE_DIR))
    def test_isoform_fallback_to_canonical(self):
        self.assertAlmostEqual(lookup_esm1b('P99999-2', 1, 'A', 'M'), -7.012)

    @override_settings(ESM1B_DATA_DIR=str(FIXTURE_DIR))
    def test_wt_mismatch_returns_none(self):
        self.assertIsNone(lookup_esm1b('P99999', 1, 'A', 'A'))

    @override_settings(ESM1B_DATA_DIR=str(FIXTURE_DIR))
    def test_missing_file_returns_none(self):
        self.assertIsNone(lookup_esm1b('P00000', 1, 'A', 'M'))


class GenomicMappingTest(SimpleTestCase):
    def test_plus_strand_block(self):
        data = parse_cds_mappings([
            {
                'seq_region_name': 'chr6',
                'start': 1000,
                'end': 1009,
                'strand': 1,
            }
        ])
        self.assertEqual(data['exon_regions'], [(0, 10)])
        self.assertEqual(data['genomic_blocks'][0]['chrom'], '6')
        self.assertEqual(cds_index_to_genomic(0, data['genomic_blocks']), ('6', 1000, 1))
        self.assertEqual(cds_index_to_genomic(3, data['genomic_blocks']), ('6', 1003, 1))
        snv = coding_edit_to_genomic_snv('ATGCCCCCCC', 0, 'G', data['genomic_blocks'])
        self.assertEqual(snv, ('6', 1000, 'A', 'G'))

    def test_minus_strand_block(self):
        data = parse_cds_mappings([
            {
                'seq_region_name': 'X',
                'start': 500,
                'end': 509,
                'strand': -1,
            }
        ])
        self.assertEqual(cds_index_to_genomic(0, data['genomic_blocks']), ('X', 509, -1))
        self.assertEqual(cds_index_to_genomic(2, data['genomic_blocks']), ('X', 507, -1))
        snv = coding_edit_to_genomic_snv('ATGCCCCCCC', 0, 'G', data['genomic_blocks'])
        self.assertEqual(snv, ('X', 509, 'T', 'C'))

    def test_wc_complement(self):
        self.assertEqual(wc_complement('A'), 'T')
        self.assertEqual(wc_complement('G'), 'C')


class CaddClientTest(SimpleTestCase):
    def test_parse_list_of_dicts(self):
        records = parse_cadd_payload([
            {'Chrom': '1', 'Pos': '100', 'Ref': 'A', 'Alt': 'G', 'RawScore': '1.2', 'PHRED': '23.4'},
            {'Chrom': '1', 'Pos': '100', 'Ref': 'A', 'Alt': 'C', 'RawScore': '0.1', 'PHRED': '2.2'},
        ])
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0]['alt'], 'G')
        self.assertAlmostEqual(records[0]['phred'], 23.4)

    def test_parse_header_rows(self):
        records = parse_cadd_payload([
            ['Chrom', 'Pos', 'Ref', 'Alt', 'RawScore', 'PHRED'],
            ['6', '200', 'C', 'T', '0.5', '12.0'],
        ])
        self.assertEqual(records[0]['alt'], 'T')
        self.assertAlmostEqual(records[0]['phred'], 12.0)

    def test_cache_hit_skips_http(self):
        tmp = Path(self._tmp())
        with override_settings(CADD_CACHE_DIR=str(tmp)):
            path = tmp / 'GRCh38-v1.7' / '1_55.json'
            path.parent.mkdir(parents=True)
            path.write_text(json.dumps([
                {'chrom': '1', 'pos': 55, 'ref': 'A', 'alt': 'G', 'phred': 18.0, 'raw': 0.9},
            ]), encoding='utf-8')
            http_get = Mock(side_effect=AssertionError('should not hit network'))
            scores, warning = fetch_cadd_scores([('1', 55, 'A', 'G')], http_get=http_get)
            self.assertIsNone(warning)
            self.assertAlmostEqual(scores[('1', 55, 'A', 'G')]['phred'], 18.0)
            http_get.assert_not_called()

    def test_api_failure_does_not_abort(self):
        with override_settings(CADD_CACHE_DIR=self._tmp()):
            http_get = Mock(side_effect=RuntimeError('down'))
            scores, warning = fetch_cadd_scores([('2', 10, 'C', 'T')], http_get=http_get)
            self.assertEqual(scores, {})
            self.assertIn('CADD scores unavailable', warning or '')

    def test_alt_selection_from_api(self):
        payload = [
            {'Chrom': '7', 'Pos': 99, 'Ref': 'A', 'Alt': 'C', 'PHRED': '4.0', 'RawScore': '0.1'},
            {'Chrom': '7', 'Pos': 99, 'Ref': 'A', 'Alt': 'G', 'PHRED': '22.0', 'RawScore': '1.1'},
        ]
        response = Mock()
        response.status_code = 200
        response.json.return_value = payload
        response.raise_for_status = Mock()
        with override_settings(CADD_CACHE_DIR=self._tmp(), CADD_MAX_WORKERS=1):
            scores, warning = fetch_cadd_scores(
                [('7', 99, 'A', 'G')],
                http_get=Mock(return_value=response),
            )
        self.assertIsNone(warning)
        self.assertAlmostEqual(scores[('7', 99, 'A', 'G')]['phred'], 22.0)

    def test_429_retries_with_backoff(self):
        limited = Mock()
        limited.status_code = 429
        ok = Mock()
        ok.status_code = 200
        ok.json.return_value = [
            {'Chrom': '3', 'Pos': 8, 'Ref': 'C', 'Alt': 'T', 'PHRED': '9.5', 'RawScore': '0.2'},
        ]
        ok.raise_for_status = Mock()
        sleeps = []
        with override_settings(CADD_CACHE_DIR=self._tmp(), CADD_MAX_WORKERS=1):
            scores, warning = fetch_cadd_scores(
                [('3', 8, 'C', 'T')],
                http_get=Mock(side_effect=[limited, ok]),
                sleeper=sleeps.append,
            )
        self.assertIsNone(warning)
        self.assertAlmostEqual(scores[('3', 8, 'C', 'T')]['phred'], 9.5)
        self.assertEqual(sleeps, [1.5])

    def _tmp(self) -> str:
        import tempfile
        return tempfile.mkdtemp(prefix='cadd-test-')


class FlattenAndJoinTest(SimpleTestCase):
    def test_flatten_keeps_outcome_details(self):
        details = [{
            'codon': 'GAA', 'aa': 'E', 'score': 0.81, 'is_pathogenic': True,
            'esm1b': -5.5, 'snvs': [('6', 100, 'A', 'G')],
        }]
        df = pd.DataFrame([{
            'position': 8,
            'WT_codon': 'AAA',
            'WT_AA': 'K',
            'ABE_mutated_AA': 'G',
            'ABE_mutated_AA_score': 0.4,
            'sgRNA_1_seq': 'A' * 20 + 'AGG',
            'sgRNA_1_protospacer': 'A' * 20,
            'sgRNA_1_pam': 'AGG',
            'sgRNA_1_strand': '+',
            'sgRNA_1_TargetApos': [5],
            'sgRNA_1_numAsWindow': 2,
            'sgRNA_1_outcome_details': details,
            'sgRNA_1_local_mut_aa': 'E',
            'sgRNA_1_local_am_score': 0.81,
            'sgRNA_1_avg_alpha_score': 0.81,
        }])
        rows = flatten_candidate_dataframe(df, 'ABE', 5)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]['mut_aa'], 'E')
        self.assertAlmostEqual(rows[0]['alpha_score'], 0.81)
        self.assertEqual(rows[0]['_outcome_details'][0]['aa'], 'E')

    def test_finalize_builds_lists_and_per_guide_max(self):
        rows = [{
            'mut_aa': 'G',
            'alpha_score': 0.4,
            '_outcome_details': [
                {
                    'codon': 'GGA', 'aa': 'G', 'score': 0.4, 'is_pathogenic': True,
                    'esm1b': -1.0, 'snvs': [('1', 11, 'C', 'T')],
                },
                {
                    'codon': 'GAA', 'aa': 'E', 'score': 0.9, 'is_pathogenic': True,
                    'esm1b': -8.0, 'snvs': [('1', 10, 'A', 'G')],
                },
            ],
        }]
        score_map = {
            ('1', 10, 'A', 'G'): {'phred': 20.0, 'raw': 1.0},
            ('1', 11, 'C', 'T'): {'phred': 10.0, 'raw': 0.4},
        }
        finalize_guide_scores(rows, score_map)
        self.assertEqual(rows[0]['mut_aa'], 'E')
        self.assertAlmostEqual(rows[0]['alpha_score'], 0.9)
        self.assertAlmostEqual(rows[0]['esm1b_score'], -8.0)
        self.assertAlmostEqual(rows[0]['cadd_phred'], 20.0)
        self.assertTrue(rows[0]['alpha_score_significant'])
        self.assertTrue(rows[0]['esm1b_score_significant'])
        self.assertTrue(rows[0]['cadd_phred_significant'])
        self.assertTrue(any('GAA(E)[0.900]*' in line for line in rows[0]['outcomes']))
        self.assertFalse(any('GGA(G)[0.400]*' in line for line in rows[0]['outcomes']))
        self.assertTrue(any('GAA(E)[-8.000]*' in line for line in rows[0]['esm1b_outcomes']))
        self.assertFalse(any('GGA(G)[-1.000]*' in line for line in rows[0]['esm1b_outcomes']))
        self.assertTrue(any('GAA(E)[20.00]*' in line for line in rows[0]['cadd_outcomes']))
        self.assertFalse(any('GGA(G)[10.00]*' in line for line in rows[0]['cadd_outcomes']))
        self.assertNotIn('_outcome_details', rows[0])

    def test_author_cutoff_stars_and_scalar_flags(self):
        """Outcomes just above/below each author cutoff get or omit *."""
        rows = [{
            'mut_aa': 'A',
            'alpha_score': 0.563,
            '_outcome_details': [
                {
                    'codon': 'GCA', 'aa': 'A', 'score': 0.563, 'is_pathogenic': True,
                    'esm1b': -7.49, 'snvs': [('1', 1, 'A', 'G')],
                },
                {
                    'codon': 'GCG', 'aa': 'A', 'score': 0.564, 'is_pathogenic': True,
                    'esm1b': -7.5, 'snvs': [('1', 2, 'A', 'G')],
                },
            ],
        }]
        finalize_guide_scores(rows, {
            ('1', 1, 'A', 'G'): {'phred': 19.99, 'raw': 0.9},
            ('1', 2, 'A', 'G'): {'phred': 20.0, 'raw': 1.0},
        })
        am_lines = rows[0]['outcomes']
        esm_lines = rows[0]['esm1b_outcomes']
        cadd_lines = rows[0]['cadd_outcomes']
        self.assertTrue(any(line.endswith('*') and '0.564' in line for line in am_lines))
        self.assertTrue(any(not line.endswith('*') and '0.563' in line for line in am_lines))
        self.assertTrue(any(line.endswith('*') and '-7.500' in line for line in esm_lines))
        self.assertTrue(any(not line.endswith('*') and '-7.490' in line for line in esm_lines))
        self.assertTrue(any(line.endswith('*') and '20.00' in line for line in cadd_lines))
        self.assertTrue(any(not line.endswith('*') and '19.99' in line for line in cadd_lines))
        self.assertTrue(rows[0]['alpha_score_significant'])
        self.assertTrue(rows[0]['esm1b_score_significant'])
        self.assertTrue(rows[0]['cadd_phred_significant'])

        below = [{
            'mut_aa': 'G',
            'alpha_score': 0.5,
            '_outcome_details': [{
                'codon': 'GGA', 'aa': 'G', 'score': 0.5, 'is_pathogenic': True,
                'esm1b': -7.0, 'snvs': [('3', 1, 'C', 'T')],
            }],
        }]
        finalize_guide_scores(below, {
            ('3', 1, 'C', 'T'): {'phred': 15.0, 'raw': 0.2},
        })
        self.assertFalse(below[0]['alpha_score_significant'])
        self.assertFalse(below[0]['esm1b_score_significant'])
        self.assertFalse(below[0]['cadd_phred_significant'])
        self.assertEqual(below[0]['outcomes'], ['GGA(G)[0.500]'])
        self.assertEqual(below[0]['esm1b_outcomes'], ['GGA(G)[-7.000]'])
        self.assertEqual(below[0]['cadd_outcomes'], ['GGA(G)[15.00]'])

    def test_cadd_uses_guide_outcomes_not_residue_best(self):
        rows = [{
            'mut_aa': 'E',
            'alpha_score': 0.9,
            '_outcome_details': [{
                'codon': 'GGA', 'aa': 'G', 'score': 0.4, 'is_pathogenic': True,
                'esm1b': -2.0, 'snvs': [('2', 50, 'A', 'G')],
            }],
        }]
        finalize_guide_scores(rows, {
            ('2', 50, 'A', 'G'): {'phred': 12.5, 'raw': 0.3},
        })
        self.assertAlmostEqual(rows[0]['cadd_phred'], 12.5)
        self.assertEqual(rows[0]['cadd_outcomes'], ['GGA(G)[12.50]'])

    def test_generate_candidates_attaches_esm_and_snvs(self):
        cds = ('N' * 16) + ('A' * 20) + 'AGG' + 'NN'
        patho_df = pd.DataFrame([
            {'protein_variant': 'K8E', 'pathogenicity score': 0.9, 'a.a.1': 'K', 'position': '8', 'a.a.2': 'E'},
            {'protein_variant': 'K8G', 'pathogenicity score': 0.4, 'a.a.1': 'K', 'position': '8', 'a.a.2': 'G'},
        ])
        exon_data = parse_cds_mappings([{
            'seq_region_name': '1',
            'start': 1000,
            'end': 1000 + len(cds) - 1,
            'strand': 1,
        }])

        def stub_esm(_uid, _pos, mut, _wt=None):
            return -5.5 if mut == 'E' else -1.0

        df = generate_candidates(
            cds,
            patho_df,
            'ABE',
            0.5,
            'ENST00000000000',
            4,
            8,
            pam_type='NGG',
            exon_data=exon_data,
            uniprot_id='P99999',
            esm_lookup=stub_esm,
        )
        self.assertFalse(df.empty)
        rows = flatten_candidate_dataframe(df, 'ABE', 5)
        self.assertTrue(rows)
        self.assertTrue(any(r.get('_outcome_details') for r in rows))
        score_map = {}
        for row in rows:
            for detail in row.get('_outcome_details') or []:
                for snv in detail.get('snvs') or []:
                    score_map[snv] = {'phred': 15.0, 'raw': 0.5}
        finalize_guide_scores(rows, score_map)
        self.assertTrue(any(r.get('esm1b_score') == -5.5 for r in rows))
        self.assertTrue(any(r.get('cadd_phred') == 15.0 for r in rows))
        self.assertTrue(any(r.get('cadd_outcomes') for r in rows))
