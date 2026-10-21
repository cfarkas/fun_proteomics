import json
import sys
import tempfile
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
sys.path.insert(0,str(ROOT/'pipeline_lib'))
from config_server import absolute_config
from reporting import summarize
from validate_repo import source_files

class PortabilityTests(unittest.TestCase):
 def test_empty_editor_preserves_missing_paths(self):
  cfg=absolute_config(ROOT/'assets/empty_study.json')
  self.assertEqual(cfg['samples'],[])
  self.assertEqual(cfg['database']['fasta'],'')
  self.assertTrue(all(v=='' for v in cfg['tools'].values()))
  self.assertIsNone(cfg['organism'])
 def test_nonhuman_proteins_retained_without_species_filter(self):
  rows=[{'Protein':'P1','Gene':'GENE1','Organism':'Mus musculus','Sample Spectral Count':'3','Sample Intensity':'50'}]
  samples=[{'sample_id':'Sample','condition':'Control','replicate':1}]
  counts,groups,matrices=summarize(list(rows[0]),rows,samples,'tissue','DDA',None)
  self.assertEqual(counts[0]['direct_protein_groups'],1)
  self.assertEqual(len(matrices['spectral_counts']),1)
 def test_publication_set_excludes_local_artifacts(self):
  paths={str(f.relative_to(ROOT)) for f in source_files(ROOT)}
  self.assertIn('assets/empty_study.json',paths)
  self.assertNotIn('config.local.json',paths)
  self.assertTrue(all(not p.startswith(('results','work/','.runtime/','ui_configs/')) for p in paths))
