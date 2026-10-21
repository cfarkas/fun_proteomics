"""Reuse refuses incomplete or changed cached searches before any model runs."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
import importlib.util
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
HAS_PARQUET = importlib.util.find_spec('pyarrow') is not None


@unittest.skipUnless(HAS_PARQUET, 'DIA reuse needs the DIA Parquet environment')
class ReuseTests(unittest.TestCase):
    def setUp(self):
        from reanalyze_dia import cached_report, same_file_content
        self.cached_report=cached_report; self.same_file_content=same_file_content
        temp=tempfile.TemporaryDirectory();self.addCleanup(temp.cleanup);self.root=Path(temp.name)
        self.task=self.root/'work/task';(self.task/'cohort').mkdir(parents=True)
        self.original=self.task/'cohort/report.parquet';self.original.write_bytes(b'fixture report')
        self.report=self.root/'published/report.parquet';self.report.parent.mkdir();self.report.write_bytes(self.original.read_bytes())
        self.report.with_name('report.manifest.txt').write_text(json.dumps({'Main .parquet report':str(self.original)})+'\n')
        self.report.with_name('engine.log').write_text('Finished\n')
        (self.task/'.exitcode').write_text('0\n')

    def test_published_and_cached_bytes_must_match(self):
        self.assertEqual(self.cached_report(self.report),(self.task,self.original))
        self.report.write_bytes(b'changed report')
        with self.assertRaisesRegex(ValueError,'Integrity mismatch'): self.cached_report(self.report)

    def test_failed_or_missing_cached_output_refused(self):
        (self.task/'.exitcode').write_text('1')
        with self.assertRaisesRegex(ValueError,'did not complete'): self.cached_report(self.report)
        self.original.unlink()
        with self.assertRaisesRegex(ValueError,'missing or ambiguous'): self.cached_report(self.report)

    def test_method_or_fasta_changes_refused(self):
        left=self.root/'method.toml';right=self.root/'old.toml'
        left.write_text('fdr=.01');right.write_text('fdr=.05')
        with self.assertRaisesRegex(ValueError,'Integrity mismatch'):self.same_file_content(left,right)
