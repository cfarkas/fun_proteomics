import copy
import json
from pathlib import Path
import sqlite3
import sys
import tempfile
from types import SimpleNamespace
import unittest
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'pipeline_lib'))
sys.path.insert(0, str(ROOT / 'tests'))
import pipeline
import reporting
from make_fixture import make_fixture


class TableTests(unittest.TestCase):
    def setUp(self):
        self.samples = [{'sample_id': 'AnySample', 'condition': 'Drug', 'replicate': 1}]
        self.rows = [{'Protein': 'P1', 'Protein ID': 'P1', 'Organism': 'Homo sapiens', 'Gene': 'A B',
                      'Indistinguishable Proteins': 'P2,P3', 'AnySample Spectral Count': '2', 'AnySample Intensity': '90'},
                     {'Protein': 'P4', 'Organism': 'Homo sapiens', 'AnySample Spectral Count': '0', 'AnySample Intensity': '50'},
                     {'Protein': 'contam_P5', 'Organism': 'Homo sapiens', 'AnySample Spectral Count': '5', 'AnySample Intensity': '900'}]

    def summarize(self):
        return reporting.summarize(list(self.rows[0]), self.rows, self.samples, 'cell', 'DDA+', 'Homo sapiens')

    def test_single_sample_intensity_never_becomes_maxlfq(self):
        counts, cohorts, matrices = self.summarize()
        self.assertEqual(counts[0]['direct_protein_groups'], 1)
        self.assertEqual(counts[0]['intensity_protein_groups'], 2)
        self.assertEqual(counts[0]['maxlfq_protein_groups'], '')
        self.assertFalse(counts[0]['maxlfq_available'])
        self.assertNotIn('maxlfq_intensity', matrices)
        self.assertEqual(cohorts[0]['groups_in_combined_table'], 2)
        self.assertEqual(matrices['spectral_counts'][0]['Gene'], 'A B')

    def test_mbr_and_direct_counts_are_distinct(self):
        for row, value in zip(self.rows, ['0', '50', '900']):
            row['AnySample MaxLFQ Intensity'] = value
        counts, _, _ = self.summarize()
        self.assertEqual(counts[0]['maxlfq_without_direct_ms2'], 1)
        self.assertEqual(counts[0]['direct_without_maxlfq'], 1)

    def test_incomplete_cohort_quantification_is_rejected(self):
        self.samples.append({'sample_id': 'Second', 'condition': 'Drug', 'replicate': 2})
        for row in self.rows:
            row['Second Spectral Count'] = '1'
            row['Second Intensity'] = '10'
        with self.assertRaisesRegex(ValueError, 'multi-run'):
            self.summarize()

    def test_wrong_sample_column_is_rejected(self):
        self.samples[0]['sample_id'] = 'Wrong'
        with self.assertRaises(ValueError):
            self.summarize()

    def test_negative_or_fractional_counts_are_rejected(self):
        for value in ['-1', '2.5', 'inf']:
            self.rows[0]['AnySample Spectral Count'] = value
            with self.assertRaises(ValueError):
                self.summarize()


class InputTests(unittest.TestCase):
    def validate_fixture(self, config):
        args = SimpleNamespace(config=config, schema=ROOT/'assets/config.schema.json', base=config.parent)
        pipeline.validate(args)

    def test_duplicate_ids_and_dia_frames_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            config = make_fixture(tmp)
            cfg = json.loads(config.read_text())
            cfg['samples'][1]['sample_id'] = cfg['samples'][0]['sample_id']
            config.write_text(json.dumps(cfg))
            with self.assertRaisesRegex(ValueError, 'Duplicate sample'):
                self.validate_fixture(config)
            config = make_fixture(tmp)
            cfg = json.loads(config.read_text())
            dbfile = Path(cfg['samples'][0]['raw_path'])/'analysis.tdf'
            with sqlite3.connect(dbfile) as db:
                db.execute('UPDATE Frames SET MsMsType=9 WHERE MsMsType=8')
            with self.assertRaisesRegex(ValueError, 'DDA-PASEF'):
                self.validate_fixture(config)

    def test_archive_traversal_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            archive = Path(tmp)/'unsafe.zip'
            with zipfile.ZipFile(archive, 'w') as z:
                z.writestr('../escape.txt', 'bad')
            with self.assertRaises(ValueError):
                pipeline.safe_extract(archive, Path(tmp)/'tools')
            self.assertFalse((Path(tmp)/'escape.txt').exists())


if __name__ == '__main__':
    unittest.main()
