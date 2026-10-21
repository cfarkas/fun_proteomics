import copy
import json
import csv
import os
import shutil
import subprocess
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'scripts'),str(ROOT/'pipeline_lib')]
from continue_run import select_session
from search_checkpoint import identity, candidate_matches, find_completed, write_checkpoint


class ContinueTests(unittest.TestCase):
    def setUp(self):
        temp=tempfile.TemporaryDirectory();self.addCleanup(temp.cleanup)
        self.root=Path(temp.name)

    def test_session_selection_includes_success_and_ignores_other_studies_and_previews(self):
        session='17ae0893-c3fe-4e41-9181-5ea3feefbc1b'
        (self.root/'.nextflow/cache'/session).mkdir(parents=True)
        config=self.root/'study with spaces.json';config.write_text(json.dumps({'outdir':'results'}))
        other=self.root/'other.json';other.write_text(json.dumps({'outdir':'other_results'}))
        command="nextflow run main.nf --config 'study with spaces.json' -w 'work with spaces'"
        line=f'date\t1m\tcompleted\tOK\trev\t{session}\t{command}\n'
        (self.root/'.nextflow/history').write_text(line+line.replace(command,command+' -preview')+
            line.replace(command,'nextflow run main.nf --config other.json'))
        found=select_session(self.root,self.root,config)
        self.assertEqual(found['session'],session)
        self.assertEqual(found['workdir'],str(self.root/'work with spaces'))
        self.assertIsNone(select_session(self.root,self.root,config,'new_output')['session'])

    def fixture(self):
        self.run=self.root/'work/ab/123/run';self.run.mkdir(parents=True)
        self.cfg={'outdir':str(self.root/'results'),'samples':[{'sample_id':'S1','raw_path':'/inputs/S1.d','cohort':'cells','condition':'Control'}],
            'raw_audit':[{'sample_id':'S1','analysis_tdf_sha256':'tdf','raw_files':{'analysis.tdf_bin':{'size':100,'mtime_ns':1}}}]}
        self.props={'database.db-path':'/database.fasta','workflow.threads':'1','msfragger.search_enzyme_name_1':'Trypsin'}
        self.versions={'MSFragger':{'version':'4.4.1','archive_sha256':'tools'},'runtime':{'python':'test'}}
        (self.run/'SEARCH_COMPLETE').write_text('complete\n')
        (self.run/'combined_protein.tsv').write_text('Protein\nP1\n')
        (self.run/'combined_peptide.tsv').write_text('Peptide\nABC\n')
        (self.run/'S1').mkdir();(self.run/'S1/psm.tsv').write_text('Peptide\nABC\n')
        (self.run/'run_provenance.json').write_text(json.dumps({'cohort':'cells','mode':'DDA','database_sha256':'db','stub':False,'dry_run':False}))
        (self.run.parent/'plan.json').write_text(json.dumps(self.cfg))
        (self.run/'fragpipe.workflow').write_text('\n'.join(k+'='+v for k,v in self.props.items()))
        (self.run/'software_versions.json').write_text(json.dumps(self.versions))

    def wanted(self):
        return copy.deepcopy(identity(self.cfg,'cells','DDA',self.props,'db',self.versions))

    def test_legacy_completed_search_reused_after_downstream_and_resource_changes(self):
        self.fixture()
        self.cfg.update(differential={'fdr':.01},optional_analyses={'gene_counts':False})
        self.cfg['samples'][0]['condition']='Renamed condition'
        self.props['workflow.threads']='8';self.props['database.db-path']='/elsewhere/db.fasta'
        self.props['workdir']='/another/workdir';self.props['fragpipe-config.tools-folder']='/another/tools'
        self.props['workflow.saved-with-ver']='24.0'
        self.assertTrue(candidate_matches(self.run,self.wanted()))
        self.assertEqual(find_completed(self.cfg,'cells','DDA',self.props,'db',self.versions,self.root/'work'),self.run)

    def test_changed_search_inputs_are_not_reused(self):
        self.fixture()
        for field,value in [('database_sha256','changed'),('mode','DDA+'),('analysis','immunopeptidomics')]:
            wanted=self.wanted();wanted[field]=value
            self.assertFalse(candidate_matches(self.run,wanted),field)
        wanted=self.wanted();wanted['samples'][0]['audit']['raw_files']['analysis.tdf_bin']['size']=101
        self.assertFalse(candidate_matches(self.run,wanted))
        wanted=self.wanted();wanted['properties']['msfragger.search_enzyme_name_1']='nonspecific'
        self.assertFalse(candidate_matches(self.run,wanted))

    def test_staged_symlink_uses_original_search_plan(self):
        self.fixture()
        reporting=self.root/'work/cd/456';reporting.mkdir(parents=True)
        (reporting/'run').symlink_to(self.run,target_is_directory=True)
        (reporting/'plan.json').write_text('{}')
        self.assertTrue(candidate_matches(reporting/'run',self.wanted()))

    def test_checkpoints_verify_output_content_and_reject_partial_or_synthetic_runs(self):
        self.fixture()
        write_checkpoint(self.run,self.cfg,'cells','DDA',self.props,'db',self.versions)
        self.assertTrue(candidate_matches(self.run,self.wanted()))
        (self.run/'S1/psm.tsv').write_text('Peptide\nXYZ\n')
        self.assertFalse(candidate_matches(self.run,self.wanted()))
        (self.run/'SEARCH_COMPLETE').unlink()
        with self.assertRaisesRegex(ValueError,'incomplete'):candidate_matches(self.run,self.wanted())

    def test_missing_raw_audit_cannot_authorize_reuse(self):
        self.fixture();self.cfg.pop('raw_audit')
        with self.assertRaisesRegex(ValueError,'acquisition audit'):self.wanted()


@unittest.skipUnless((ROOT/'.runtime/env/bin/nextflow').is_file(), 'Local Nextflow runtime required')
class NextflowContinueTests(unittest.TestCase):
    def test_continue_uses_native_cache_after_an_intervening_preview(self):
        with tempfile.TemporaryDirectory(prefix='nextflow_continue_') as temp:
            root=Path(temp);(root/'scripts').mkdir()
            shutil.copy2(ROOT/'scripts/continue_run.py',root/'scripts/continue_run.py')
            shutil.copy2(ROOT/'nextflow.config',root/'nextflow.config')
            (root/'study.json').write_text(json.dumps({'outdir':str(root/'results')}))
            (root/'main.nf').write_text('nextflow.enable.dsl = 2\nprocess CACHE_CHECK {\noutput:\npath "done.txt"\nscript:\n"echo synthetic > done.txt"\n}\nworkflow { CACHE_CHECK() }\n')
            env=dict(os.environ,JAVA_HOME=str(ROOT/'.runtime/env'),NXF_HOME=str(ROOT/'.runtime/nextflow'),
                     NXF_VER='26.04.2',NXF_DISABLE_CHECK_LATEST='true')
            base=[str(ROOT/'.runtime/env/bin/nextflow'),'run','main.nf','--runtime_dir',str(ROOT/'.runtime'),
                  '--config','study.json','-ansi-log','false']
            for extra in (['-with-trace','first.tsv'],['-preview'],['--continue','-with-trace','continued.tsv']):
                result=subprocess.run(base+extra,cwd=root,env=env,capture_output=True,text=True,timeout=60)
                self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            with (root/'continued.tsv').open() as handle:rows=list(csv.DictReader(handle,delimiter='\t'))
            self.assertEqual(len(rows),1)
            self.assertEqual(rows[0]['status'],'CACHED')
            history=(root/'.nextflow/history').read_text().splitlines()
            self.assertEqual(history[0].split('\t')[5],history[-1].split('\t')[5])


if __name__=='__main__':unittest.main()
