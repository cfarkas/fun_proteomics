import csv,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'pipeline_lib'))
from final_results import peptide_counts,gene_counts
class PeptideExportTests(unittest.TestCase):
 def test_filters_proteins_and_preserves_sample_order(self):
  with tempfile.TemporaryDirectory() as d:
   root=Path(d);(root/'combined_peptide.tsv').write_text('Peptide Sequence\tProtein\tGene\tS2 Spectral Count\tS1 Spectral Count\nPEPTIDE\tP1\tG1\t7\t3\nCONTAM\tcontam_P2\tC\t5\t5\n')
   output=root/'peptides.csv';peptide_counts(root,[{'sample_id':'S1'},{'sample_id':'S2'}],{'P1'},output)
   rows=list(csv.DictReader(output.read_text().splitlines()));self.assertEqual(rows,[{'Peptide':'PEPTIDE','Protein':'P1','Gene':'G1','Mapped Genes':'','S1':'3','S2':'7'}])

 def test_gene_counts_exclude_ambiguous_assignments(self):
  rows=[{'Peptide':'A','Gene':'G1','Mapped Genes':'G1','S1':'3','S2':'0'},
        {'Peptide':'B','Gene':'G1','Mapped Genes':'','S1':'2','S2':'4'},
        {'Peptide':'C','Gene':'G1','Mapped Genes':'G2','S1':'100','S2':'200'},
        {'Peptide':'D','Gene':'','Mapped Genes':'','S1':'8','S2':'9'}]
  counts,excluded=gene_counts(rows,['S1','S2'])
  self.assertEqual(counts,[{'FeatureID':'G1','Gene':'G1','S1':5,'S2':4}])
  self.assertEqual(len(excluded),2)
 def test_legacy_peptides_require_mapping_evidence(self):
  with self.assertRaisesRegex(ValueError,'Mapped Genes'):gene_counts([{'Peptide':'A','Gene':'G1','S1':'2'}],['S1'])
