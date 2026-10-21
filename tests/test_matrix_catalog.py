"""Matrix exports preserve exact cell text and expose both file formats."""
import csv
from pathlib import Path
import sys
import tempfile
import unittest
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'pipeline_lib'))
from matrix_catalog import export_matrix_formats, build_matrix_catalog


def cells(path):
    with path.open(newline='', encoding='utf-8') as handle:
        return list(csv.reader(handle, delimiter='\t' if path.suffix == '.tsv' else ','))


class MatrixCatalogTests(unittest.TestCase):
    def test_lossless_text_and_missing_values(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/'matrix.csv'
            values = [['FeatureID','S1','S2','annotation'],['00012','0','','comma, tab\t newline\nµ'],['NA','1.234567890123456789','1e-300','quote"']]
            with path.open('w',newline='') as handle:
                csv.writer(handle).writerows(values)
            export_matrix_formats(path)
            self.assertEqual(cells(path.with_suffix('.tsv')), values)
            export_matrix_formats(path.with_suffix('.tsv'))
            self.assertEqual(cells(path), values)

    def test_invalid_table_does_not_replace_existing_companion(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'broken.csv'; path.write_text('ID,S1\nP1,2,3\n')
            path.with_suffix('.tsv').write_text('unchanged')
            with self.assertRaisesRegex(ValueError,'Ragged'):
                export_matrix_formats(path)
            self.assertEqual(path.with_suffix('.tsv').read_text(),'unchanged')

    def test_all_layouts_and_idempotent_catalog(self):
        with tempfile.TemporaryDirectory() as temp:
            output=Path(temp)
            locations=['main_results/counts/cells/dda/raw/protein_counts.csv',
                       'optional_analyses/statistics/counts/cells/dda/normalized/isoform_log2_cpm.csv',
                       'matrices/tissue/gene_intensities.csv',
                       'differential_expression/tissue/gene/log2_intensities.csv',
                       'benchmark/tissue/service_quantification.csv']
            for name in locations:
                path=output/name; path.parent.mkdir(parents=True,exist_ok=True)
                path.write_text('ID,S1\n0001,0\nP2,\n')
            first=build_matrix_catalog(output)
            second=build_matrix_catalog(output)
            self.assertEqual(first,second)
            self.assertEqual(len(first),5)
            for row in first:
                self.assertEqual(cells(output/row['csv_path']),cells(output/row['tsv_path']))
                self.assertIn(row['tsv_path'],(output/'MATRICES.html').read_text())
            self.assertTrue((output/'MATRICES.tsv').is_file())
            self.assertTrue((output/locations[0]).with_suffix('.tsv').is_file())
