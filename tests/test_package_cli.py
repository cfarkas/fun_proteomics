"""Portable launchers must reject stale inputs and select the correct saved run."""
import argparse
import contextlib
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
if not os.environ.get('FUN_PROTEOMICS_TEST_INSTALLED'):
    sys.path.insert(0, str(ROOT/'src'))
from fun_proteomics import cli
from fun_proteomics.continuation import dia_settings, nextflow_session, verify_delivery, verify_recorded_inputs
from fun_proteomics.project import digest, verify_files


class PackageContinuationTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory(prefix='portable study ')
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)

    def delivery(self):
        final = self.root/'final_results'
        final.mkdir()
        (final/'matrix.tsv').write_text('Protein\tS1\nP1\t10\n')
        (final/'FINAL_RESULTS_COMPLETE').write_text('complete\n')
        (final/'OUTPUT_SHA256SUMS').write_text(''.join(digest(p)+'  '+p.name+'\n' for p in sorted(final.iterdir())))
        return final

    def test_delivery_checks_every_file_and_completion_marker(self):
        final = self.delivery()
        self.assertEqual(verify_delivery(self.root), [])
        (final/'matrix.tsv').write_text('corrupted')
        self.assertIn('matrix.tsv', verify_delivery(self.root))
        (final/'extra.csv').write_text('unexpected')
        self.assertIn('extra.csv', verify_delivery(self.root))
        (final/'FINAL_RESULTS_COMPLETE').unlink()
        self.assertTrue(verify_delivery(self.root))

    def test_manifest_cannot_escape_or_duplicate_files(self):
        final = self.delivery()
        manifest = final/'OUTPUT_SHA256SUMS'
        good = manifest.read_text()
        for invalid in ('a'*64+'  ../outside\n', good+good, 'z'*64+'  matrix.tsv\n'):
            manifest.write_text(invalid)
            with self.assertRaises(ValueError): verify_delivery(self.root)
        manifest.write_text('a'*64+'  symlink\n')
        (final/'symlink').symlink_to(self.root/'outside')
        with self.assertRaisesRegex(ValueError,'Invalid'): verify_delivery(self.root)

    def test_raw_reuse_rejects_changed_content_and_allows_identical_touch(self):
        audit = self.root/'audit'; audit.mkdir()
        (audit/'run.json').write_text('{}')
        raw = self.root/'raw.bin'; raw.write_bytes(b'original')
        item = dict(path=str(raw), bytes=raw.stat().st_size, mtime_ns=raw.stat().st_mtime_ns, sha256=digest(raw))
        (audit/'raw_integrity.json').write_text(json.dumps([{'files':[item]}]))
        verify_recorded_inputs(self.root, hash_raw=True)
        os.utime(raw, ns=(item['mtime_ns']+1000000,item['mtime_ns']+1000000))
        verify_recorded_inputs(self.root)
        raw.write_bytes(b'tampered')
        os.utime(raw, ns=(item['mtime_ns']+2000000,item['mtime_ns']+2000000))
        with self.assertRaisesRegex(ValueError,'Raw acquisition changed'): verify_recorded_inputs(self.root)

    def test_hash_raw_requires_a_recorded_baseline(self):
        audit=self.root/'audit'; audit.mkdir(); (audit/'run.json').write_text('{}')
        raw=self.root/'raw.bin'; raw.write_bytes(b'synthetic')
        item=dict(path=str(raw),bytes=raw.stat().st_size,mtime_ns=raw.stat().st_mtime_ns)
        (audit/'raw_integrity.json').write_text(json.dumps([{'files':[item]}]))
        verify_recorded_inputs(self.root)
        with self.assertRaisesRegex(ValueError,'No recorded raw SHA256'): verify_recorded_inputs(self.root,True)

    def test_count_sample_metadata_is_an_input(self):
        final=self.delivery()
        samples=self.root/'samples.csv';samples.write_text('sample_id,condition\nS1,A\n')
        provenance=self.root/'analysis_details/reporting/result_audit/main_statistics/input_provenance.json'
        provenance.parent.mkdir(parents=True)
        provenance.write_text(json.dumps([dict(path=str(samples),sha256=digest(samples))]))
        verify_recorded_inputs(self.root)
        samples.write_text('sample_id,condition\nS1,B\n')
        with self.assertRaisesRegex(ValueError,'count analysis input'): verify_recorded_inputs(self.root)

    def test_legacy_session_matches_config_and_output_ignores_previews(self):
        cache=self.root/'.nextflow/cache'; cache.mkdir(parents=True)
        commands=[('correct',"nextflow run main.nf --dia_config 'my study.json' --outdir results"),
                  ('other-study',"nextflow run main.nf --dia_config other.json --outdir results"),
                  ('other-output',"nextflow run main.nf --dia_config 'my study.json' --outdir elsewhere"),
                  ('preview',"nextflow run main.nf --dia_config 'my study.json' --outdir results -preview")]
        for session, _ in commands: (cache/session).mkdir()
        (self.root/'.nextflow/history').write_text(''.join('date\tduration\tname\tOK\trev\t'+session+'\t'+cmd+'\n' for session,cmd in commands))
        command=['nextflow','run',str(self.root/'main.nf'),'--dia_config',str(self.root/'my study.json'),'--outdir',str(self.root/'results')]
        self.assertEqual(nextflow_session(self.root,command),'correct')
        command[-1]=str(self.root/'missing')
        self.assertIsNone(nextflow_session(self.root,command))

    def test_partial_dia_cohort_cannot_skip_unfinished_searches(self):
        prepared=self.root/'prepared';prepared.mkdir()
        (prepared/'samples.tsv').write_text('sample_type\nserum\ntissue\n')
        (prepared/'search.fasta').write_text('>P1\nPEPTIDE\n')
        report=self.root/'search/serum/cohort/report.parquet';report.parent.mkdir(parents=True);report.touch()
        self.assertIsNone(dia_settings(self.root))
        report=self.root/'search/tissue/cohort/report.parquet';report.parent.mkdir(parents=True);report.touch()
        self.assertEqual(dia_settings(self.root)['source_run'],str(self.root))

    def test_owned_nextflow_is_not_redirected_into_dia_reporting(self):
        controller=self.root/'.fun_proteomics';controller.mkdir()
        launch=controller/'nextflow';(launch/'.nextflow/cache/session').mkdir(parents=True)
        state=dict(engine='nextflow',outdir=str(self.root),project=str(ROOT),status='failed',session='session',launch_dir=str(launch),command=['nextflow','run',str(ROOT/'main.nf')])
        (controller/'run.json').write_text(json.dumps(state))
        args=cli.parser().parse_args(['--continue',str(self.root),'--print-command'])
        with patch.object(cli,'dia_settings',side_effect=AssertionError('should use Nextflow')), patch.object(cli,'nextflow_binary',return_value='nextflow'), patch.object(cli,'execute',return_value=0) as execute:
            self.assertEqual(cli.continue_run(args,ROOT),0)
            self.assertEqual(execute.call_args.args[0][-2:],['-resume','session'])

    def test_active_process_cannot_launch_again(self):
        controller=self.root/'.fun_proteomics';controller.mkdir()
        (controller/'run.json').write_text(json.dumps(dict(status='running',pid=os.getpid(),process_identity=cli.process_identity(os.getpid()))))
        args=cli.parser().parse_args(['--continue',str(self.root)])
        with self.assertRaisesRegex(ValueError,'still active'):cli.continue_run(args,ROOT)

    def test_subprocess_failure_is_recorded_and_returned(self):
        state=self.root/'.fun_proteomics/run.json'
        with contextlib.redirect_stdout(io.StringIO()):
            code=cli.execute([sys.executable,'-c','raise SystemExit(7)'],self.root,os.environ.copy(),state,dict(engine='native-counts',outdir=str(self.root)))
        self.assertEqual(code,7)
        self.assertEqual(json.loads(state.read_text())['status'],'failed')

    def test_native_flags_reach_script_help(self):
        with patch.object(cli,'resolve_project',return_value=ROOT), patch.object(cli.subprocess,'run') as run:
            run.return_value.returncode=0
            self.assertEqual(cli.main(['quantify','--help']),0)
            self.assertEqual(run.call_args.args[0][-1],'--help')

    def test_workflow_checksum_detects_damage(self):
        path=self.root/'main.nf';path.write_text('workflow {}')
        manifest={'main.nf':digest(path)}
        verify_files(self.root,manifest)
        path.write_text('modified')
        with self.assertRaisesRegex(ValueError,'checksum mismatch'):verify_files(self.root,manifest)


if __name__=='__main__':unittest.main()
