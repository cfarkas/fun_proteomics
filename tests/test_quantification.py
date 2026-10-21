"""Scientific regression tests for imports, annotations, contrasts and detection."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'pipeline_lib'))
sys.path.insert(0,str(ROOT/'modules/diapasef'))
from quantification import load_matrix, annotate_groups, gene_matrix, benchmark_detection
from quantification_analysis import analyze_matrix, differential_options, validate_design
HAS_PARQUET = importlib.util.find_spec("pyarrow") is not None


class QuantificationTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.root=Path(self.tmp.name)
        self.samples=pd.DataFrame({'sample_id':['S1','S2'],'run':['run1','run2']})

    def test_arbitrary_metadata_windows_paths_and_empty_cells(self):
        p=self.root/'service.csv'
        pd.DataFrame({'Genes':['G1','G2'],'Protein.Group':['P2;P1','P3'],r'C:\raw\run1.d':[10,0],'/raw/run2.d':['',20]}).to_csv(p,index=False)
        matrix,ann,mapping=load_matrix(p,self.samples)
        self.assertEqual(list(matrix.index),['P1;P2','P3']);self.assertTrue(np.isnan(matrix.loc['P1;P2','S2']))
        self.assertEqual(ann.loc['P3','Genes'],'G2');self.assertEqual(list(mapping.sample_id),['S1','S2'])

    def test_fragpipe_explicit_mapping_and_group_canonicalization(self):
        p=self.root/'provider.tsv'
        pd.DataFrame({'Protein':['sp|P1|ONE','rev_sp|P2|TWO'],'Indistinguishable Proteins':['sp|P3|THREE',''],
                      'Gene':['A','B'],'provider1 Spectral Count':[3,99],'provider2 Spectral Count':[0,99]}).to_csv(p,sep='\t',index=False)
        samples=self.samples.assign(service_sample=['provider1','provider2'])
        matrix,_,mapping=load_matrix(p,samples,'spectral_count')
        self.assertEqual(list(matrix.index),['P1;P3']);self.assertEqual(matrix.loc['P1;P3','S1'],3)
        self.assertEqual(mapping.attrs['excluded_decoy_contaminant_rows'],1)
        with self.assertRaisesRegex(ValueError,'matched 0'):load_matrix(p,self.samples)

    def test_invalid_values_groups_and_columns_fail(self):
        p=self.root/'matrix.csv'
        for value in ['-1','inf','oops']:
            p.write_text('Protein.Group,S1,S2\nP1,'+value+',2\n')
            with self.assertRaises(ValueError):load_matrix(p,self.samples)
        p.write_text('Protein.Group,S1,S1\nP1,1,2\n')
        with self.assertRaisesRegex(ValueError,'duplicate column'):load_matrix(p,self.samples)
        p.write_text('Protein.Group,S1,S2\nP1;P2,1,2\nP2;P1,3,4\n')
        with self.assertRaisesRegex(ValueError,'duplicate protein'):load_matrix(p,self.samples)

    def test_gene_fallback_ambiguity_and_missing_rollup(self):
        fasta=self.root/'targets.fasta';fasta.write_text('>sp|P1|ONE GN=A\nMPEPTIDE\n>sp|P2|TWO GN=B\nMPEPTIDE\n')
        report=pd.DataFrame({'Protein.Group':['P1','P1;P2','P9','P1;P8'],'Genes':['','','','']})
        ann=annotate_groups(report,report['Protein.Group'],fasta)
        self.assertEqual(list(ann.gene_mapping_status),['single_gene','ambiguous','unmapped','partial'])
        matrix=pd.DataFrame({'S1':[10,50,20,60],'S2':[np.nan,30,40,60]},index=ann.index)
        genes=gene_matrix(matrix,ann)
        self.assertEqual(list(genes.index),['A']);self.assertEqual(genes.loc['A','S1'],10);self.assertTrue(np.isnan(genes.loc['A','S2']))

    def test_benchmark_gains_losses_and_zero_denominator(self):
        new=pd.DataFrame({'S1':[2,0,3],'S2':[0,1,0]},index=['P1','P2','P3'])
        old=pd.DataFrame({'S1':[4,5],'S2':[0,0]},index=['P1','P4'])
        summary,members=benchmark_detection(new,old);summary=summary.set_index('sample_id')
        self.assertEqual(summary.loc['S1',['shared','pipeline_only','service_only','net_difference']].tolist(),[1,1,1,0])
        self.assertTrue(np.isnan(summary.loc['S2','percent_change']))
        self.assertEqual(summary.loc['__union__','pipeline_detected'],3)
        self.assertEqual(set(members.status),{'shared','pipeline_only','service_only'})

    def test_design_validation(self):
        options=differential_options({'enabled':True})
        with self.assertRaisesRegex(ValueError,'condition and replicate'):validate_design(self.samples,options)
        samples=pd.DataFrame({'sample_id':['A1','A2','B1','B2'],'condition':['A','A','B','B'],'replicate':['1','1','1','2']})
        with self.assertRaisesRegex(ValueError,'unique biological'):validate_design(samples,options)
        with self.assertRaises(ValueError):differential_options({'enabled':True,'min_replicates':1})
        with self.assertRaises(ValueError):differential_options({'enabled':'false'})

    def test_limma_all_pairs_direction_fdr_and_missingness(self):
        rng=np.random.default_rng(7)
        samples=pd.DataFrame({'sample_id':[f'S{i}' for i in range(9)],'condition':['A']*3+['B']*3+['C']*3,'replicate':['1','2','3']*3})
        x=rng.normal(15,.15,(80,9));x[:10,3:6]+=3;x[10:20,6:]+=2
        x[-1,3:]=np.nan
        matrix=pd.DataFrame(2**x,index=pd.Index([f'P{i}' for i in range(80)],name='Protein.Group'),columns=samples.sample_id)
        annotations=pd.DataFrame({'Genes':[f'G{i}' for i in range(80)],'gene_mapping_status':'single_gene'},index=matrix.index)
        output=self.root/'results'
        analyze_matrix(matrix,annotations,samples,output,{'enabled':True,'rscript':os.environ.get('PROTEOMICS_RSCRIPT','Rscript')})
        index=pd.read_csv(output/'differential_expression/protein/results_index.csv')
        self.assertEqual(len(index),3)
        result=pd.read_csv(output/'differential_expression/protein/A_vs_B/all_results.csv').set_index('FeatureID')
        self.assertTrue((result.loc[[f'P{i}' for i in range(10)],'logFC']>2.5).all())
        self.assertTrue((result.loc[[f'P{i}' for i in range(10)],'FDR']<.05).all())
        self.assertTrue(np.isnan(result.loc['P79','PValue']));self.assertEqual(result.loc['P79','test_status'],'insufficient observations')
        tested=result.PValue.dropna().sort_values();expected=np.minimum.accumulate((tested.values*len(tested)/np.arange(1,len(tested)+1))[::-1])[::-1].clip(0,1)
        np.testing.assert_allclose(result.loc[tested.index,'FDR'],expected,rtol=1e-10)
        actual=pd.read_csv(output/'matrices/protein_intensities.csv');self.assertTrue(actual.loc[79,'S8']!=actual.loc[79,'S8'])
        self.assertTrue((output/'differential_expression/gene/A_vs_C/all_results.csv').is_file())
        self.assertTrue((output/'differential_expression/gene/log2_intensities.tsv').is_file())
        self.assertTrue((output/'matrices/protein_intensities.tsv').is_file())
        self.assertEqual(len(pd.read_csv(output/'results_index.csv')), 6)

    @unittest.skipUnless(HAS_PARQUET, "Parquet integration runs in the pinned DIA test environment")
    def test_real_dia_export_and_collect_without_spectrum_search(self):
        from common import filter_report, read_method
        samples=pd.DataFrame({'sample_id':['S1','S2'],'sample_type':['tissue']*2,'run':['run1','run2'],'library_reference':[True,False]})
        sheet=self.root/'samples.tsv';samples.to_csv(sheet,sep='\t',index=False)
        rows=[]
        for i,run in enumerate(samples.run):
            rows.append({'Run':run,'Protein.Group':'P1','Protein.Names':'ONE','Genes':'GENE1','PG.MaxLFQ':100+i,
                         'Precursor.Id':'PEPTIDE2','Stripped.Sequence':'PEPTIDE','Q.Value':.001,'Global.Q.Value':.001,
                         'Lib.Q.Value':.001,'PG.Q.Value':.001,'Global.PG.Q.Value':.001,'Lib.PG.Q.Value':.001})
        report=self.root/'report.parquet';pd.DataFrame(rows).to_parquet(report)
        method=ROOT/'modules/diapasef/method.toml';group=self.root/'tissue'
        subprocess.run([sys.executable,str(ROOT/'modules/diapasef/05_count_proteins.py'),'--samples',str(sheet),'--report',str(report),'--method',str(method),'--outdir',str(group)],check=True,capture_output=True,text=True)
        subprocess.run([sys.executable,str(ROOT/'modules/diapasef/06_collect_results.py'),'--samples',str(sheet),'--groups',str(group)],cwd=self.root,check=True,capture_output=True,text=True)
        matrix=pd.read_csv(self.root/'final_results/matrices/tissue/protein_intensities.csv')
        self.assertEqual(matrix.loc[0,'Genes'],'GENE1');self.assertEqual(matrix.loc[0,'S1'],100)
        self.assertTrue((self.root/'final_results/MATRICES.html').is_file())
        counterpart=pd.read_csv(self.root/'final_results/matrices/tissue/protein_intensities.tsv',sep='\t')
        pd.testing.assert_frame_equal(matrix,counterpart)
        self.assertTrue((group/'protein_matrix_sample_ids.csv').is_file())
        bad=pd.DataFrame(rows);bad.loc[0,'Q.Value']=-.1;bad.loc[1,'Protein.Group']='contam_P2'
        self.assertTrue(filter_report(bad,read_method(method)).empty)


if __name__=='__main__':unittest.main()
