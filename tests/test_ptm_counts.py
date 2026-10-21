import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'pipeline_lib'))
from ptm_counts import export

class PTMCounts(unittest.TestCase):
 def fixture(self,root):
  run=root/'run';(run/'S1').mkdir(parents=True);(run/'S2').mkdir()
  (run/'SEARCH_COMPLETE').touch()
  config={'organism':None,'samples':[{'sample_id':s,'cohort':'cell','condition':'Control','replicate':i} for i,s in enumerate(['S1','S2'],1)]}
  base={'Spectrum':'scan1','Peptide':'MAMK','Assigned Modifications':'1M(15.9949)','Protein':'P1','Is Decoy':'false','Is Contaminant':'false','Gene':'GENE1','Mapped Genes':'','Mapped Proteins':'','Charge':'2'}
  rows=[base,{**base,'Spectrum':'scan2'}, {**base,'Spectrum':'scan3','Charge':'3','Assigned Modifications':'3M(15.9949)'},
   {**base,'Spectrum':'scan4','Peptide':'MCAMK','Assigned Modifications':'1M(15.9949), 2C(57.02146), 4M(15.9949)'},
   {**base,'Spectrum':'scan5','Peptide':'MCK','Assigned Modifications':'2C(57.02146)'},
   {**base,'Spectrum':'scan6','Peptide':'MAMK','Mapped Genes':'GENE2'},
   {**base,'Spectrum':'scan7','Gene':''}, {**base,'Spectrum':'scan8','Is Decoy':'true'},
   {**base,'Spectrum':'scan9','Is Contaminant':'true'}, {**base,'Protein':'P2'}]
  for sample,rr in [('S1',rows),('S2',[{**base,'Assigned Modifications':''}])]:
   with (run/sample/'psm.tsv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(base),delimiter='\t');w.writeheader();w.writerows(rr)
  return config,run
 def read(self,path):
  with path.open() as f:return list(csv.DictReader(f))
 def test_unique_modified_sequences_and_ambiguity(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);cfg,run=self.fixture(root);out=root/'counts';export(cfg,run,'cell',out)
   row=self.read(out/'gene_modified_peptide_counts.csv')[0];self.assertEqual(row,{'Gene':'GENE1','S1':'2','S2':'0'})
   self.assertEqual(self.read(out/'gene_modified_peptidoform_counts.csv')[0]['S1'],'3')
   self.assertEqual(self.read(out/'gene_total_identified_peptide_counts.csv')[0]['S1'],'3')
   self.assertEqual(self.read(out/'gene_modification_peptide_counts.csv')[0]['S1'],'2')
   self.assertEqual(len(self.read(out/'ambiguous_modified_peptides.csv')),2)
   # Multi-M oxidation counts once per PSM, even though two residues modified.
   self.assertEqual(self.read(out/'ptm_psm_counts.csv')[0]['S1'],'6')
 def test_protein_counts_collapse_spectra_and_exclude_shared_assignments(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);cfg,run=self.fixture(root);out=root/'counts'
   path=run/'S1/psm.tsv'
   with path.open() as h:reader=csv.DictReader(h,delimiter='\t');fields=reader.fieldnames;rows=list(reader)
   rows.append({**rows[0],'Spectrum':'shared','Peptide':'SHARED','Mapped Proteins':'P2'})
   rows.append({**rows[0],'Spectrum':'phospho','Peptide':'MAMK','Assigned Modifications':'2S(79.9663)'})
   with path.open('w') as h:w=csv.DictWriter(h,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
   export(cfg,run,'cell',out)
   types=self.read(out/'protein_modification_peptide_counts.csv')
   self.assertEqual(next(r for r in types if r['Modification']=='Oxidation (M)')['S1'],'2')
   self.assertEqual(next(r for r in types if r['Modification']=='Phosphorylation (STY)')['S1'],'1')
   self.assertEqual(self.read(out/'ambiguous_protein_modified_peptides.csv')[0]['Peptide'],'SHARED')
 def test_incomplete_search_rejected(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);cfg,run=self.fixture(root);(run/'SEARCH_COMPLETE').unlink()
   with self.assertRaisesRegex(ValueError,'not complete'):export(cfg,run,'cell',root/'counts')
 def test_missing_psm_not_reported_as_zero(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);cfg,run=self.fixture(root);(run/'S2/psm.tsv').unlink()
   with self.assertRaises(FileNotFoundError):export(cfg,run,'cell',root/'counts')
 def test_organism_filter(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);cfg,run=self.fixture(root);cfg['organism']='Mouse'
   (run/'combined_protein.tsv').write_text('Protein\tOrganism\nP1\tHuman\n')
   export(cfg,run,'cell',root/'counts');self.assertEqual(self.read(root/'counts/gene_modified_peptide_counts.csv'),[])
if __name__=='__main__':unittest.main()
