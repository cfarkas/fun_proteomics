"""Assignments must govern analysis rows, provenance, exclusions and provider matches."""
import copy
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
if not os.environ.get('FUN_PROTEOMICS_TEST_INSTALLED'):sys.path.insert(0,str(ROOT/'src'))
from fun_proteomics import assignments, study


class AssignmentTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(prefix='proteomics-assignment-');self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        for src,dst in [('quantification_matrix.csv','matrix.csv'),('quantification_samples.csv','samples.csv'),('service_quantification.csv','service.csv')]:
            shutil.copy2(ROOT/'examples'/src,self.root/dst)
        self.cfg=study.normalize(dict(project=str(self.root/'study'),matrix=str(self.root/'matrix.csv'),independent_samples=True))

    def test_matrix_discovery_starts_unassigned_and_preserves_source_when_renamed(self):
        board=assignments.discover(self.cfg,ROOT)
        self.assertEqual([r['sample_id'] for r in board['rows']],['C1','C2','A1','A2','B1','B2'])
        self.assertTrue(all(not r['condition'] for r in board['rows']))
        cfg=dict(self.cfg,assignment_board=board)
        with self.assertRaisesRegex(ValueError,'Assign at least one'):study.check(cfg,ROOT)
        for i,row in enumerate(board['rows']):row['condition']='Control' if i<2 else 'Treatment'
        board['rows'][0]['sample_id']='Renamed1'
        board['rows'][-1]['condition']=''
        report=study.check(cfg,ROOT)
        self.assertEqual((report['samples'],report['comparisons'],report['unassigned']),(5,1,1))
        self.assertEqual(study.samples_for(cfg,ROOT).iloc[0].run,'C1')
        path,_=study.save(cfg,ROOT)
        self.assertEqual(study.load(path)['assignment_board'],board)
        import pandas as pd
        for suffix,sep in [('csv',','),('tsv','\t')]:
            snapshot=pd.read_csv(path.parent/('sample_assignments.'+suffix),sep=sep,keep_default_na=False)
            self.assertEqual(len(snapshot),6);self.assertFalse(snapshot.iloc[-1].included)
        self.assertFalse((self.root/'study/results').exists())

    def test_imported_table_is_prefill_but_board_has_authority(self):
        cfg=dict(self.cfg,samples=str(self.root/'samples.csv'))
        board=assignments.discover(cfg,ROOT)
        self.assertEqual(board['conditions'],['Control','Treatment A','Treatment B'])
        for row in board['rows']:
            if row['condition']=='Treatment B':row['condition']=''
        cfg['assignment_board']=board
        report=study.check(cfg,ROOT)
        self.assertEqual(report['samples'],4);self.assertEqual(report['comparisons'],1)
        self.assertEqual(list(study.samples_for(cfg,ROOT).source_replicate),['1','2','1','2'])

    def test_renaming_to_another_column_name_keeps_the_original_quantities(self):
        cfg=dict(self.cfg,samples=str(self.root/'samples.csv'))
        board=assignments.discover(cfg,ROOT)
        board['rows'][0]['sample_id']='C2';board['rows'][1]['sample_id']='C1'
        frame=study.samples_for(dict(cfg,assignment_board=board),ROOT)
        from quantification import load_matrix
        values,_,mapping=load_matrix(cfg['matrix'],frame)
        self.assertEqual(list(mapping.source_column)[:2],['C1','C2'])
        self.assertEqual(list(values.columns)[:2],['C2','C1'])

    def test_changed_inputs_or_quantity_invalidates_discovery(self):
        board=assignments.discover(self.cfg,ROOT);cfg=dict(self.cfg,assignment_board=board)
        board['rows'][0]['condition']='Control'
        with self.assertRaisesRegex(ValueError,'changed after sample discovery'):
            study.samples_for(dict(cfg,quantity='value'),ROOT)
        p=self.root/'matrix.csv';p.write_text(p.read_text()+'\n')
        with self.assertRaisesRegex(ValueError,'changed after sample discovery'):study.samples_for(cfg,ROOT)

    def test_unknown_duplicate_or_empty_assignments_refused(self):
        board=assignments.discover(self.cfg,ROOT)
        board['rows'][0]['condition']='Control'
        for field,value in [('conditions',['Control','Control']),('conditions',['']),('rows',[])]:
            bad=copy.deepcopy(board);bad[field]=value
            with self.assertRaises(ValueError):assignments.rows_for(dict(self.cfg,assignment_board=bad))
        bad=copy.deepcopy(board);bad['rows'][1]['sample_id']=bad['rows'][0]['sample_id']
        with self.assertRaisesRegex(ValueError,'unique'):assignments.rows_for(dict(self.cfg,assignment_board=bad))
        bad=copy.deepcopy(board);bad['rows'][0]['condition']='Unknown'
        with self.assertRaisesRegex(ValueError,'unknown condition'):assignments.rows_for(dict(self.cfg,assignment_board=bad))

    def test_provider_mapping_follows_selected_and_renamed_samples(self):
        cfg=dict(self.cfg,samples=str(self.root/'samples.csv'),service_matrix=str(self.root/'service.csv'),service_samples=str(self.root/'samples.csv'))
        board=assignments.discover(cfg,ROOT)
        for row in board['rows']:
            if row['condition']=='Treatment B':row['condition']=''
        board['rows'][0]['sample_id']='Patient1';cfg['assignment_board']=board
        frame=study.samples_for(cfg,ROOT);provider=assignments.service_samples_for(cfg,frame,ROOT)
        self.assertEqual(list(provider.sample_id),['Patient1','C2','A1','A2'])
        self.assertEqual(provider.iloc[0].service_sample,'C1')
        report=study.check(cfg,ROOT)
        self.assertEqual(len(report['benchmark']),5)

    def test_fragpipe_columns_discovered_once_and_raw_tables_keep_absolute_paths(self):
        p=self.root/'fragpipe.tsv';p.write_text('Protein\tGenes\tC1 MaxLFQ Intensity\tC1 Intensity\tC2 MaxLFQ Intensity\nP1\tGENE\t1\t2\t3\n')
        cfg=dict(self.cfg,matrix=str(p));board=assignments.discover(cfg,ROOT)
        self.assertEqual([r['sample_id'] for r in board['rows']],['C1','C2'])
        raw=self.root/'raw';raw.mkdir();(raw/'S1.d').mkdir();(raw/'S2.d').mkdir();(raw/'BLANK.d').mkdir()
        cfg=study.normalize(dict(project=str(self.root/'raw-study'),mode='dia',raw_root=str(raw),differential=False))
        board=assignments.discover(cfg,ROOT)
        self.assertEqual(len(board['rows']),2)
        self.assertEqual(board['rows'][0]['raw_path'],str(raw/'S1.d'))
        board['rows'][0]['condition']='Control'
        frame=study.samples_for(dict(cfg,assignment_board=board),ROOT)
        self.assertEqual(len(frame),1);self.assertEqual(frame.iloc[0].cohort,'study')

    def test_assignment_export_can_prefill_a_new_project(self):
        cfg=dict(self.cfg,samples=str(self.root/'samples.csv'))
        cfg['assignment_board']=assignments.discover(cfg,ROOT)
        cfg['assignment_board']['rows'][0]['sample_id']='Renamed1'
        path,_=study.save(cfg,ROOT)
        new=dict(self.cfg,project=str(self.root/'second-study'),samples=str(path.parent/'sample_assignments.csv'))
        new['assignment_board']=assignments.discover(new,ROOT)
        self.assertEqual(new['assignment_board']['rows'][0]['detected_sample'],'C1')
        second,_=study.save(new,ROOT)
        self.assertTrue((second.parent/'sample_assignments.tsv').exists())

    def test_editing_a_saved_draft_keeps_previous_assignments(self):
        cfg=dict(self.cfg,samples=str(self.root/'samples.csv'))
        cfg['assignment_board']=assignments.discover(cfg,ROOT)
        path,_=study.save(cfg,ROOT);original=path.read_bytes()
        cfg['assignment_board']['rows'][0]['sample_id']='Renamed1'
        study.save(cfg,ROOT)
        backups=list((path.parent/'backups').glob('*/run.json'))
        self.assertEqual(len(backups),1);self.assertEqual(backups[0].read_bytes(),original)
        self.assertEqual(study.load(path)['assignment_board']['rows'][0]['sample_id'],'Renamed1')
        (path.parent.parent/'.fun_proteomics/launch.json').write_text('{}')
        cfg['assignment_board']['rows'][0]['sample_id']='Another'
        with self.assertRaisesRegex(ValueError,'already been launched'):study.save(cfg,ROOT)

    def test_matrix_path_basename_collision_requires_explicit_mapping(self):
        p=self.root/'ambiguous.tsv';p.write_text('Protein.Group\t/run1/S1.d\t/run2/S1.d\nP1\t1\t2\n')
        with self.assertRaisesRegex(ValueError,'unique'):assignments.discover(dict(self.cfg,matrix=str(p)),ROOT)


if __name__=='__main__':unittest.main()
