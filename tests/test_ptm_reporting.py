import csv
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'pipeline_lib'))
sys.path.insert(0,str(Path(__file__).resolve().parent))
import test_ptm_counts
from reporting import report

class ReportIntegration(unittest.TestCase):
 def test_real_report_exports_gene_tables(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);cfg,run=test_ptm_counts.PTMCounts().fixture(root)
   (run/'combined_protein.tsv').write_text('Protein\tGene\tS1 Spectral Count\tS2 Spectral Count\tS1 MaxLFQ Intensity\tS2 MaxLFQ Intensity\nP1\tGENE1\t4\t1\t100\t100\n')
   (run/'combined_peptide.tsv').write_text('Peptide Sequence\tProtein\tGene\tS1 Spectral Count\tS2 Spectral Count\nPEPTIDE\tP1\tGENE1\t4\t1\n')
   (run/'run_provenance.json').write_text('{}')
   plan=root/'plan.json';plan.write_text(json.dumps(cfg))
   out=root/'tables';report(SimpleNamespace(plan=plan,run=run,cohort='cell',mode='DDA+',out=out))
   with (out/'ptm_counts/gene_modified_peptide_counts.csv').open() as f:
    self.assertEqual(next(csv.DictReader(f)),{'Gene':'GENE1','S1':'2','S2':'0'})
   self.assertTrue((out/'count_tables.xlsx').is_file())
