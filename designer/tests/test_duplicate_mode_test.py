from django.test import TestCase

from designer.pipeline import (
    annotate_sgrna_binding,
    apply_duplicate_mode,
    normalize_duplicate_mode,
)


def _row(position, seq, score):
    return {
        'position': position,
        'sgrna_seq': seq,
        'alpha_score': score,
        'pam': 'NGG',
        'strand': '+',
    }


class DuplicateModeTest(TestCase):

    def test_normalize_legacy_hide(self):
        self.assertEqual(normalize_duplicate_mode('hide'), 'best')

    def test_unique_keeps_only_single_site_guides(self):
        rows = [
            _row(1, 'SEQ-A', 0.9),
            _row(2, 'SEQ-A', 0.8),
            _row(3, 'SEQ-B', 0.7),
        ]
        result = apply_duplicate_mode(rows, 'unique')
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0]['sgrna_seq'], 'SEQ-B')
        self.assertEqual(result[0]['binding_count'], 1)

    def test_best_keeps_highest_score_per_guide(self):
        rows = [
            _row(10, 'SEQ-A', 0.7),
            _row(20, 'SEQ-A', 0.95),
            _row(30, 'SEQ-B', 0.6),
        ]
        result = apply_duplicate_mode(rows, 'best')
        self.assertEqual(len(result), 2)
        by_seq = {r['sgrna_seq']: r for r in result}
        self.assertEqual(by_seq['SEQ-A']['position'], 20)
        self.assertEqual(by_seq['SEQ-A']['binding_count'], 2)
        self.assertEqual(by_seq['SEQ-B']['binding_count'], 1)

    def test_group_keeps_all_and_sorts_by_sequence(self):
        rows = [
            _row(10, 'SEQ-A', 0.7),
            _row(20, 'SEQ-A', 0.95),
            _row(30, 'SEQ-B', 0.6),
        ]
        result = apply_duplicate_mode(rows, 'group')
        self.assertEqual(len(result), 3)
        self.assertEqual([r['sgrna_seq'] for r in result], ['SEQ-A', 'SEQ-A', 'SEQ-B'])
        self.assertEqual(result[0]['binding_count'], 2)

    def test_annotate_binding_positions(self):
        rows = annotate_sgrna_binding([
            _row(5, 'SEQ-A', 0.5),
            _row(12, 'SEQ-A', 0.6),
        ])
        self.assertEqual(rows[0]['binding_positions'], [5, 12])
        self.assertEqual(rows[0]['binding_count'], 2)
