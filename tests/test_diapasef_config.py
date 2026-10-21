"""Regression checks for acquisition routing and DIA study input integrity."""
import csv
import importlib.util
import json
from pathlib import Path
import shutil
import sqlite3
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('dia_config', ROOT / 'scripts/diapasef_config.py')
dia = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dia)


def make_fixture(root):
    root = Path(root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    (root / 'raw').mkdir(exist_ok=True)
    (root / 'cache').mkdir(exist_ok=True)
    rows = []
    for group in ('serum', 'cell_pellet'):
        for index in range(2):
            raw = root / 'raw' / group / f"sample {index + 1}'s.d"
            raw.mkdir(parents=True, exist_ok=True)
            database = raw / 'analysis.tdf'
            if database.exists():
                database.unlink()
            with sqlite3.connect(database) as con:
                con.execute('CREATE TABLE Frames (MsMsType INTEGER)')
                con.executemany('INSERT INTO Frames VALUES (?)', [(0,), (9,)])
                con.execute('CREATE TABLE GlobalMetadata (Key TEXT, Value TEXT)')
                con.execute("INSERT INTO GlobalMetadata VALUES ('AcquisitionDateTime', ?)", (f'2026-01-0{index+1}',))
            (raw / 'analysis.tdf_bin').write_bytes(b'fixture-only')
            rows.append({'sample_id': f'00{index + 1}', 'sample_type': group,
                         'raw_path': str(raw.relative_to(root / 'raw')),
                         'library_reference': str(index == 0).lower()})
    with (root / 'samples.tsv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    (root / 'targets.fasta').write_text('>sp|TEST|TEST Test protein\nMPEPTIDEKPEPTIDER\n')
    header = ['Protein.Group', 'Protein.Names', 'Genes', 'Description', 'N.Sequences', 'N.Proteotypic.Sequences']
    header += [Path(r['raw_path']).name for r in rows if r['sample_type'] == 'serum']
    with (root / 'service.tsv').open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(header)
        writer.writerow(['TEST', 'TEST', 'TEST', 'fixture', 1, 1, 1, 1])
    cfg = {'samplesheet': 'samples.tsv', 'raw_root': 'raw', 'canonical_fasta': 'targets.fasta',
           'install_cache': 'cache', 'service_matrices': {'serum': 'service.tsv'},
           'expected_samples': 4, 'cpus': 2, 'memory_gib': 4, 'outdir': 'results',
           'python': shutil.which('python3.12') or '/usr/bin/python3'}
    path = root / 'study.local.json'
    path.write_text(json.dumps(cfg))
    return path


class TestDiaConfiguration(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="dia config's ")
        self.root = Path(self.temp.name)
        self.config = make_fixture(self.root)

    def tearDown(self):
        self.temp.cleanup()

    def validate(self):
        return dia.normalize_config(self.config, ROOT / 'modules/diapasef')

    def change_config(self, **updates):
        cfg = json.loads(self.config.read_text())
        cfg.update(updates)
        self.config.write_text(json.dumps(cfg))

    def change_rows(self, callback):
        sheet = self.root / 'samples.tsv'
        with sheet.open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        callback(rows)
        with sheet.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)

    def test_relative_paths_spaces_and_shared_subject_ids(self):
        cfg = self.validate()
        self.assertEqual(cfg['expected_samples'], 4)
        self.assertEqual(cfg['groups']['serum'], {'samples': 2, 'library_references': 1})
        self.assertEqual(cfg['groups']['cell_pellet'], {'samples': 2, 'library_references': 1})
        self.assertEqual(cfg['raw_root'], str(self.root / 'raw'))

    def test_mixed_acquisitions_rejected(self):
        path = next((self.root / 'raw').rglob('analysis.tdf'))
        with sqlite3.connect(path) as con:
            con.execute('INSERT INTO Frames VALUES (8)')
        with self.assertRaisesRegex(ValueError, 'expected DIA-PASEF'):
            self.validate()

    def test_dda_acquisitions_rejected(self):
        path = next((self.root / 'raw').rglob('analysis.tdf'))
        with sqlite3.connect(path) as con:
            con.execute('UPDATE Frames SET MsMsType=8 WHERE MsMsType=9')
        with self.assertRaisesRegex(ValueError, 'expected DIA-PASEF'):
            self.validate()

    def test_missing_binary_rejected(self):
        next((self.root / 'raw').rglob('analysis.tdf_bin')).unlink()
        with self.assertRaisesRegex(ValueError, 'Incomplete raw'):
            self.validate()

    def test_no_library_reference_rejected(self):
        self.change_rows(lambda rows: [row.update(library_reference='false') for row in rows])
        with self.assertRaisesRegex(ValueError, 'library reference'):
            self.validate()

    def test_duplicate_id_within_type_rejected(self):
        self.change_rows(lambda rows: rows[1].update(sample_id=rows[0]['sample_id']))
        with self.assertRaisesRegex(ValueError, 'Duplicate sample ID'):
            self.validate()

    def test_incomplete_service_matrix_rejected(self):
        self.change_config(service_matrices={'serum': 'targets.fasta'})
        with self.assertRaisesRegex(ValueError, 'service matrix'):
            self.validate()

    def test_output_must_not_overlap_raws(self):
        self.change_config(outdir='raw/results')
        with self.assertRaisesRegex(ValueError, 'separate directory trees'):
            self.validate()

    def test_count_and_unknown_option_checks(self):
        self.change_config(expected_samples=5)
        with self.assertRaisesRegex(ValueError, 'Expected 5 samples'):
            self.validate()
        self.change_config(expected_samples=4, misspelled_setting=True)
        with self.assertRaisesRegex(ValueError, 'Unknown DIA'):
            self.validate()

    def test_qvalue_cannot_be_relaxed(self):
        method = self.root / 'method.toml'
        method.write_text((ROOT / 'modules/diapasef/method.toml').read_text().replace('precursor_q = 0.01', 'precursor_q = 0.05'))
        self.change_config(method_config='method.toml')
        with self.assertRaisesRegex(ValueError, 'precursor_q'):
            self.validate()


if __name__ == '__main__':
    unittest.main()
