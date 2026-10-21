"""Result delivery preserves statistical tables, status and historical evidence."""
import csv
import json
from pathlib import Path
import re
import subprocess
import tempfile
import unittest
from unittest.mock import patch
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'pipeline_lib'))
import result_layout as layout


class ResultLayoutTests(unittest.TestCase):
    def setUp(self):
        temp=tempfile.TemporaryDirectory();self.addCleanup(temp.cleanup);self.root=Path(temp.name)
        self.out=self.root/'final_results';self.out.mkdir()
        self.raw=self.root/'analysis_details/reporting/cells/dda';self.raw.mkdir(parents=True)
        self.job={'analysis':'global-label-free','groups':[{'cohort':'cells','mode':'DDA','directory':'../analysis_details/reporting/cells/dda'}],
                  'differential':{'enabled':True,'fdr':.05},'optional_analyses':{'statistics':{'enabled':True,'fdr':.1}}}
        (self.out/'analysis_settings.json').write_text(json.dumps(self.job))
        (self.out/'FINAL_RESULTS_COMPLETE').touch()
        layout.write_rows(self.out/'sample_annotations.csv',['sample_id','condition'],[{'sample_id':'S1','condition':'Control'}])
        layout.write_rows(self.raw/'protein_counts.csv',['FeatureID','S1'],[{'FeatureID':'P1','S1':10}])
        (self.raw/'samples.csv').write_text('sample_id,condition\nS1,Control\n')
        (self.out/'R_session.txt').write_text('session fixture\n')
        self.main=self.table(self.out,feature='protein',pvalue=.001,significant=True)
        self.optional=self.table(self.out/'optional_analyses/statistics',feature='isoform',pvalue='',significant=False)
        retired=self.out/'optional_analyses/open_search/cells/dda';retired.mkdir(parents=True)
        (retired/'global.profile.tsv').write_text('historical discovery\n')
        zscore=self.root/'zscore/cells/dda';zscore.mkdir(parents=True)
        (zscore/'protein_zscore.csv').write_text('FeatureID,S1\nP1,0\n')

    def table(self,base,feature,pvalue,significant):
        path=base/f'differential/edgeR/pairwise/cells/dda/{feature}/edgeR_Control_vs_Treatment.csv'
        fields=['FeatureID','logFC','PValue','FDR','significant','test_status']
        rows=[dict(zip(fields,['P1',2,pvalue,.01 if pvalue else '',str(significant).upper(),'tested' if pvalue else 'insufficient evidence']))]
        layout.write_rows(path,fields,rows)
        sig=path.parent/'significant'/path.name
        layout.write_rows(sig,fields,rows if significant else [])
        index=[{'method':'edgeR','test_type':'pairwise','cohort':'cells','mode':'DDA','feature_type':feature,
                'comparison':'Control_vs_Treatment','tested_features':int(bool(pvalue)),'reported_features':1,
                'significant_features':int(significant),'all_results':str(path.relative_to(base/'differential')),
                'significant_results':str(sig.relative_to(base/'differential'))}]
        layout.write_rows(base/'differential/results_index.csv',list(index[0]),index)
        return path

    def test_final_layout_preserves_tables_and_archives_retired_outputs(self):
        hashes={str(self.main.relative_to(self.out)):layout.digest(self.main),str(self.optional.relative_to(self.out)):layout.digest(self.optional)}
        result=layout.publish_results(self.out,preserve_original=True)
        self.assertEqual(result['tests'],2)
        archive=Path(result['archive']);self.assertTrue((archive/'optional_analyses/open_search/cells/dda/global.profile.tsv').is_file())
        self.assertTrue(archive.with_name(archive.name+'_zscore').is_dir())
        self.assertFalse((self.root/'zscore').exists())
        self.assertFalse((self.out/'differential').exists())
        self.assertFalse((self.out/'optional_analyses/open_search').exists())
        self.assertFalse((self.out/'R_session.txt').exists())
        index=layout.read_rows(self.out/'results_index.csv')
        self.assertEqual({row['section'] for row in index},{'Main results','Optional analyses'})
        self.assertEqual(index[0]['status'],'Significant results')
        self.assertEqual(index[1]['status'],'Insufficient evidence')
        for row in index:
            self.assertTrue((self.out/row['full_results']).is_file())
            self.assertTrue((self.out/row['significant_results']).is_file())
            self.assertIn(layout.digest(self.out/row['full_results']),hashes.values())
        self.assertTrue((self.out/'main_results/counts/cells/dda/raw/protein_counts.csv').is_file())
        audit=Path(result['audit'])
        for row in layout.read_rows(audit/'file_manifest.csv'):
            path=(audit if row['area']=='audit' else self.out)/row['file']
            self.assertEqual(layout.digest(path),row['sha256'])
        self.assertTrue(layout.publish_results(self.out)['already_organized'])
        self.assertFalse(list(self.root.glob('.result_layout_*')))

    def test_missing_input_aborts_without_changing_completed_results(self):
        self.main.unlink()
        with self.assertRaisesRegex(ValueError,'Missing result file'):
            layout.publish_results(self.out,preserve_original=True)
        self.assertTrue((self.out/'analysis_settings.json').is_file())
        self.assertTrue((self.out/'optional_analyses/open_search').is_dir())
        self.assertFalse(list(self.root.glob('.result_layout_*')))

    def test_staged_search_summary_is_delivered_with_final_qc(self):
        summary=self.root/'staged_search_summary';summary.mkdir()
        (summary/'sample_detection_counts.tsv').write_text('sample_id\tpsms\nS1\t10\n')
        layout.publish_results(self.out,search_summary=summary)
        self.assertEqual((self.out/'main_results/quality_control/sample_detection_counts.tsv').read_text(),
                         (summary/'sample_detection_counts.tsv').read_text())

    def test_unfinished_results_are_not_published_as_complete(self):
        (self.out/'FINAL_RESULTS_COMPLETE').unlink()
        with self.assertRaisesRegex(ValueError,'completed'):
            layout.publish_results(self.out)

    def test_result_paths_cannot_escape_their_source_tree(self):
        index=self.out/'differential/results_index.csv';rows=layout.read_rows(index)
        rows[0]['all_results']='../../private.csv';layout.write_rows(index,list(rows[0]),rows)
        with self.assertRaisesRegex(ValueError,'outside'):
            layout.publish_results(self.out)
        self.assertTrue(self.main.is_file())

    def test_empty_and_nonsignificant_tests_remain_distinct(self):
        self.table(self.out,feature='protein',pvalue=.8,significant=False)
        self.optional.write_text(self.optional.read_text().splitlines()[0]+'\n')
        index=self.out/'optional_analyses/statistics/differential/results_index.csv';rows=layout.read_rows(index)
        rows[0].update(reported_features=0,tested_features=0,comparison=Path(rows[0]['all_results']).name);layout.write_rows(index,list(rows[0]),rows)
        layout.publish_results(self.out)
        statuses=[r['status'] for r in layout.read_rows(self.out/'results_index.csv')]
        self.assertEqual(statuses,['No significant results','No eligible evidence'])
        self.assertEqual(layout.read_rows(self.out/'results_index.csv')[1]['comparison'],'Control_vs_Treatment')

    def test_dashboard_payload_links_and_javascript_syntax(self):
        self.job['analysis']='example </script><img src=x onerror=alert(1)>'
        (self.out/'analysis_settings.json').write_text(json.dumps(self.job))
        layout.publish_results(self.out)
        html=(self.out/'START_HERE.html').read_text()
        self.assertNotIn('<img src=x',html)
        payload=re.search(r'<script id="result-data" type="application/json">(.*?)</script>',html,re.S).group(1)
        data=json.loads(payload)
        self.assertEqual(data['analyses'],self.job['analysis'])
        for row in data['results']:
            for key in ('full_results','significant_results'):
                self.assertTrue((self.out/row[key]).is_file())
        for row in data['files']:self.assertTrue((self.out/row['path']).is_file())
        script=re.findall(r'<script>(.*?)</script>',html,re.S)[0]
        path=self.root/'dashboard.js';path.write_text(script)
        result=subprocess.run(['node','--check',str(path)],text=True,capture_output=True)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertNotIn('innerHTML',script)

    def test_gene_results_are_optional_and_keep_main_model_fdr(self):
        self.main=self.table(self.out,feature='gene',pvalue=.001,significant=True)
        layout.write_rows(self.raw/'gene_counts.csv',['FeatureID','S1'],[{'FeatureID':'G1','S1':5}])
        original=layout.digest(self.main)
        layout.publish_results(self.out)
        gene=next(row for row in layout.read_rows(self.out/'results_index.csv') if row['feature_type']=='gene')
        self.assertEqual(gene['section'],'Optional analyses')
        self.assertEqual(float(gene['fdr']),.05)
        self.assertEqual(layout.digest(self.out/gene['full_results']),original)
        self.assertTrue((self.out/'optional_analyses/gene_counts/counts/cells/dda/raw/gene_counts.csv').is_file())
        self.assertFalse((self.out/'main_results/counts/cells/dda/raw/gene_counts.csv').exists())

    def test_gene_outputs_can_be_deselected(self):
        self.main=self.table(self.out,feature='gene',pvalue=.001,significant=True)
        self.job['optional_analyses']['gene_counts']=False
        (self.out/'analysis_settings.json').write_text(json.dumps(self.job))
        layout.write_rows(self.raw/'gene_counts.csv',['FeatureID','S1'],[{'FeatureID':'G1','S1':5}])
        layout.publish_results(self.out)
        self.assertFalse(any(row['feature_type']=='gene' for row in layout.read_rows(self.out/'results_index.csv')))
        self.assertFalse((self.out/'optional_analyses/gene_counts').exists())

    def test_upgrade_existing_layout_moves_gene_files_without_refitting(self):
        self.main=self.table(self.out,feature='gene',pvalue=.001,significant=True)
        layout.publish_results(self.out)
        index=layout.read_rows(self.out/'results_index.csv')
        row=next(row for row in index if row['feature_type']=='gene')
        originals={}
        moves={}
        for key in ('full_results','significant_results'):
            old=row[key];new=old.replace('optional_analyses/gene_counts/comparisons/','main_results/differential_expression/')
            target=self.out/new;target.parent.mkdir(parents=True,exist_ok=True)
            (self.out/old).rename(target);row[key]=new;moves[old]=new;originals[key]=layout.digest(target)
        row.update(section='Main results',analysis='Differential expression')
        layout.write_rows(self.out/'results_index.csv',layout.INDEX_FIELDS,index)
        audit=self.root/'analysis_details/reporting/result_audit/file_manifest.csv'
        manifest=layout.read_rows(audit)
        for entry in manifest:
            if entry['file'] in moves:entry['file']=moves[entry['file']]
        layout.write_rows(audit,list(manifest[0]),manifest)
        result=layout.publish_results(self.out,preserve_original=True)
        self.assertEqual(result['moved_files'],2)
        gene=next(row for row in layout.read_rows(self.out/'results_index.csv') if row['feature_type']=='gene')
        self.assertEqual(gene['section'],'Optional analyses')
        for key in originals:self.assertEqual(layout.digest(self.out/gene[key]),originals[key])
