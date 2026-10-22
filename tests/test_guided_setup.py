"""User-facing setup must share execution, preserve inputs, and reject unsafe reuse."""
import contextlib
import copy
import io
import json
import os
from pathlib import Path
import shutil
import sqlite3
import sys
import tempfile
import threading
import unittest
from unittest.mock import patch
from urllib.error import HTTPError
from urllib.request import Request,urlopen

ROOT=Path(__file__).resolve().parents[1]
if not os.environ.get('FUN_PROTEOMICS_TEST_INSTALLED'):sys.path.insert(0,str(ROOT/'src'))
from fun_proteomics import cli,study
from fun_proteomics.setup_web import Application,make_server


class GuidedSetupTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(prefix='proteomics-guided-')
        self.addCleanup(self.tmp.cleanup);self.root=Path(self.tmp.name)
        for name in ('quantification_matrix.csv','quantification_samples.csv','service_quantification.csv'):
            shutil.copy2(ROOT/'examples'/name,self.root/name)
        self.cfg=study.normalize(dict(project=str(self.root/'project'),mode='matrix',independent_samples=True,
             matrix=str(self.root/'quantification_matrix.csv'),samples=str(self.root/'quantification_samples.csv'),
             service_matrix=str(self.root/'service_quantification.csv')))
        self.args=cli.parser().parse_args(['setup','--project',self.cfg['project']])

    def test_check_does_not_create_project_and_reports_every_pair(self):
        report=study.check(self.cfg,ROOT)
        self.assertEqual(report['comparisons'],3)
        self.assertEqual(report['samples'],6)
        self.assertEqual(report['genes'],10)
        self.assertFalse(Path(self.cfg['project']).exists())
        self.assertEqual(report['benchmark'][-1]['sample_id'],'__union__')
        json.dumps(report,allow_nan=False)

    def test_unconfirmed_independence_refused_but_matrix_only_allowed(self):
        cfg=dict(self.cfg,independent_samples=False)
        with self.assertRaisesRegex(ValueError,'independent'):study.check(cfg,ROOT)
        report=study.check(dict(cfg,differential=False),ROOT)
        self.assertEqual(report['comparisons'],0)

    def test_missing_or_ambiguous_sample_columns_fail_before_save(self):
        samples=Path(self.cfg['samples']);samples.write_text(samples.read_text().replace('C1,','Unknown,'))
        with self.assertRaisesRegex(ValueError,'matched 0 columns'):study.save(self.cfg,ROOT)
        self.assertFalse(Path(self.cfg['project']).exists())

    def test_save_and_reload_is_identical_and_launched_settings_are_locked(self):
        path,_=study.save(self.cfg,ROOT)
        original=path.read_bytes()
        study.save(study.load(path),ROOT)
        self.assertEqual(path.read_bytes(),original)
        (path.parent.parent/'.fun_proteomics/launch.json').write_text('{}')
        with self.assertRaisesRegex(ValueError,'new project'):study.save(dict(self.cfg,normalization='median'),ROOT)
        self.assertEqual(path.read_bytes(),original)

    def test_unrelated_existing_folder_is_not_adopted(self):
        p=Path(self.cfg['project']);p.mkdir();(p/'important.txt').write_text('keep')
        with self.assertRaisesRegex(ValueError,'empty'):study.save(self.cfg,ROOT)
        self.assertEqual((p/'important.txt').read_text(),'keep')

    def test_saved_configuration_cannot_be_silently_moved(self):
        path,_=study.save(self.cfg,ROOT)
        moved=self.root/'moved';shutil.copytree(path.parent.parent,moved)
        with self.assertRaisesRegex(ValueError,'recorded project location'):study.load(moved/'config/run.json')

    def test_an_active_project_cannot_be_resaved(self):
        path,_=study.save(self.cfg,ROOT)
        with study.locked(path.parent.parent):
            with self.assertRaisesRegex(ValueError,'already running'):study.save(self.cfg,ROOT)

    def test_independent_assignment_retains_original_replicates(self):
        frame=study.samples_for(self.cfg,ROOT)
        self.assertEqual(list(frame.source_replicate),['1','2','1','2','1','2'])
        self.assertTrue(frame.replicate.is_unique)

    def test_print_command_does_not_run_or_create_native_inputs(self):
        path,_=study.save(self.cfg,ROOT);args=cli.parser().parse_args(['--print-command','run','--config',str(path)])
        with patch.object(cli,'main') as execute,contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(study.run_study(path,args,ROOT),0)
        execute.assert_not_called()
        self.assertFalse((path.parent/'samples.csv').exists())
        self.assertFalse((path.parent.parent/'results').exists())
        self.assertFalse((path.parent.parent/'.fun_proteomics/launch.json').exists())

    def test_zero_service_detections_have_valid_json_preview(self):
        import pandas as pd
        p=Path(self.cfg['service_matrix']);frame=pd.read_csv(p)
        for name in ['C1','C2','A1','A2','B1','B2']:frame[name]=0
        frame.to_csv(p,index=False)
        result=study.check(self.cfg,ROOT)
        self.assertIsNone(result['benchmark'][-1]['percent_change'])
        json.dumps(result,allow_nan=False)

    def test_raw_scan_and_acquisition_validation(self):
        raw=self.root/'raw';raw.mkdir()
        rows=[]
        for i in range(4):
            d=raw/f'S{i}.d';d.mkdir();(d/'analysis.tdf_bin').write_bytes(b'synthetic')
            with sqlite3.connect(d/'analysis.tdf') as db:
                db.execute('CREATE TABLE Frames (MsMsType INTEGER)');db.executemany('INSERT INTO Frames VALUES (?)',[(0,),(9,)])
            rows.append(dict(sample_id=f'S{i}',condition='A' if i<2 else 'B',raw_path=str(d)))
        self.assertEqual(len(study.scan(raw)),4)
        fasta=self.root/'targets.fasta';fasta.write_text('>sp|P00001|DEMO GN=GENE_A\nMPEPTIDEK\n')
        cfg=study.normalize(dict(project=str(self.root/'raw-project'),mode='dia',independent_samples=True,
                  raw_root=str(raw),fasta=str(fasta),sample_rows=rows,cpus=1,memory_gib=4,python='/usr/bin/python3.12'))
        if not Path('/usr/bin/python3.12').is_file():self.skipTest('System Python 3.12 not available')
        report=study.check(cfg,ROOT)
        self.assertEqual(report['comparisons'],1)
        self.assertEqual(report['raw']['samples'],4)
        self.assertFalse(Path(cfg['project']).exists())
        with sqlite3.connect(raw/'S0.d/analysis.tdf') as db:db.execute('UPDATE Frames SET MsMsType=8 WHERE MsMsType=9')
        with self.assertRaisesRegex(ValueError,'DIA-PASEF'):study.check(cfg,ROOT)

    def test_local_r_installation_and_interactive_alias_are_reused(self):
        runtime=self.root/'workflow/.runtime/r-env/bin';runtime.mkdir(parents=True)
        r=runtime/'Rscript';r.write_text('synthetic executable')
        args=cli.parser().parse_args(['doctor'])
        with patch.dict(os.environ,{},clear=True):
            env=cli.environment(self.root/'workflow',args)
        self.assertEqual(env['PROTEOMICS_RSCRIPT'],str(r))
        with patch('fun_proteomics.study.dispatch',return_value=0) as dispatch:
            cli.main(['--rscript','/example/Rscript','--interactive','--project',str(self.root/'study')])
        self.assertEqual(dispatch.call_args.args[0].rscript,'/example/Rscript')
        self.assertEqual(dispatch.call_args.args[0].study_project,str(self.root/'study'))

    def test_raw_dda_generates_a_valid_count_workflow_design(self):
        import jsonschema
        raw=self.root/'dda';raw.mkdir();rows=[]
        for i in range(4):
            d=raw/f'S{i}.d';d.mkdir();(d/'analysis.tdf_bin').write_bytes(b'synthetic')
            with sqlite3.connect(d/'analysis.tdf') as db:
                db.execute('CREATE TABLE Frames (MsMsType INTEGER)');db.executemany('INSERT INTO Frames VALUES (?)',[(0,),(8,)])
            rows.append(dict(sample_id=f'S{i}',condition='A' if i<2 else 'B',raw_path=str(d)))
        fasta=self.root/'search.fasta';fasta.write_text('>P1 GN=A\nMPEPTIDEK\n>rev_P1 GN=A\nMEDITPEPK\n>contam_C\nMPEPTIDEK\n>rev_contam_C\nMEDITPEPK\n')
        archives={name:str(self.root/(name+'.zip')) for name in ('msfragger_zip','ionquant_zip','diatracer_zip')}
        for name in archives.values():Path(name).write_bytes(b'archive existence fixture')
        cfg=study.normalize(dict(project=str(self.root/'dda-project'),mode='dda',raw_root=str(raw),fasta=str(fasta),
                        sample_rows=rows,independent_samples=True,cpus=1,memory_gib=4,**archives))
        report=study.check(cfg,ROOT);self.assertEqual(report['comparisons'],1)
        native=study.raw_study(cfg,study.samples_for(cfg,ROOT),ROOT)
        jsonschema.validate(native,json.loads((ROOT/'assets/config.schema.json').read_text()))
        self.assertTrue(native['optional_analyses']['gene_counts'])
        self.assertEqual(native['differential']['comparisons'],'all_pairwise')
        self.assertIsInstance(native['samples'][0]['replicate'],int)

    def test_web_session_auth_and_save_use_the_same_configuration(self):
        app=Application(self.args,ROOT);server=make_server(app)
        thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
        self.addCleanup(server.server_close);self.addCleanup(server.shutdown)
        url=f'http://127.0.0.1:{server.server_port}'
        def request(path,payload=None,token=app.token,origin=None):
            headers={'X-Session-Token':token,'Content-Type':'application/json'}
            if origin:headers['Origin']=origin
            return urlopen(Request(url+path,json.dumps(payload).encode() if payload is not None else None,headers=headers),timeout=15)
        for token,origin in [('',None),(app.token,'https://example.invalid')]:
            with self.assertRaises(HTTPError) as caught:request('/api/initial',token=token,origin=origin)
            self.assertEqual(caught.exception.code,403)
        with request('/api/initial') as result:self.assertEqual(result.status,200)
        with request('/api/save',{'config':self.cfg}) as result:data=json.load(result)
        self.assertEqual(study.load(data['path']),self.cfg)
        with request('/api/check',{'config':self.cfg}) as result:self.assertEqual(json.load(result)['comparisons'],3)
        with self.assertRaises(HTTPError) as caught:request('/results/../../config/run.json')
        self.assertEqual(caught.exception.code,403)

    def test_dead_controller_is_interrupted_instead_of_blocking_resume(self):
        path,_=study.save(self.cfg,ROOT)
        state=path.parent.parent/'results/.fun_proteomics/run.json';state.parent.mkdir(parents=True)
        state.write_text(json.dumps({'status':'running','pid':999999999,'process_identity':'stale',
                                     'outdir':str(path.parent.parent/'results')}))
        self.assertEqual(study.status(path.parent.parent)['state'],'interrupted')
        app=Application(self.args,ROOT)
        self.assertEqual(app.status()['state'],'interrupted')

    def test_top_level_native_split_does_not_consume_setup_values(self):
        args=['setup','--terminal','--non-interactive','--project',str(self.root/'new'),
              '--matrix','counts','--samples','quantify','--independent-samples']
        with patch('fun_proteomics.study.dispatch',return_value=0) as dispatch:
            self.assertEqual(cli.main(args),0)
        self.assertEqual(dispatch.call_args.args[0].matrix,'counts')


if __name__=='__main__':unittest.main()
