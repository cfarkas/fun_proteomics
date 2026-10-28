"""Scientific regression: FDR choices preserve the fit and identify their null scale."""
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
from quantification_analysis import analyze_matrix
from quantification_config import differential_options


class FdrtoolTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp=tempfile.TemporaryDirectory(prefix='fdrtool-regression-')
        cls.addClassCleanup(cls.temp.cleanup)
        cls.root=Path(cls.temp.name)
        cls.rscript=os.environ.get('PROTEOMICS_RSCRIPT','Rscript')
        subprocess.run([cls.rscript,'--vanilla','-e',
            'stopifnot(requireNamespace("fdrtool",quietly=TRUE))'],check=True)
        samples=pd.DataFrame(dict(sample_id=[f'S{i}' for i in range(12)],
            condition=['A']*4+['B']*4+['C']*4,replicate=[str(i) for i in range(12)]))
        rng=np.random.default_rng(56)
        y=rng.normal(16,.4,(1000,12));y[:30,4:8]+=1.3;y[30:60,8:]-=1
        y[-1,4:]=np.nan;y[-2,:2]=np.nan
        matrix=pd.DataFrame(2**y,index=pd.Index([f'P{i}' for i in range(1000)],name='Protein.Group'),columns=samples.sample_id)
        annotations=pd.DataFrame(dict(Genes=[f'G{i}' for i in range(1000)],gene_mapping_status='single_gene'),index=matrix.index)
        for method in ['BH','fdrtool_pvalue','fdrtool_normal']:
            analyze_matrix(matrix,annotations,samples,cls.root/method,
                dict(enabled=True,gene_level=False,fdr_method=method,rscript=cls.rscript))

    def table(self,method,contrast='A_vs_B'):
        return pd.read_csv(self.root/method/'differential_expression/protein'/contrast/'all_results.csv').set_index('FeatureID')

    def test_all_methods_preserve_effects_original_p_bh_and_missingness(self):
        for contrast in ['A_vs_B','A_vs_C','B_vs_C']:
            native=self.table('BH',contrast)
            np.testing.assert_allclose(native.FDR,native.BH_FDR,equal_nan=True)
            self.assertTrue((native.FDR.dropna()+1e-14>=native.PValue.dropna()).all())
            for method in ['fdrtool_pvalue','fdrtool_normal']:
                result=self.table(method,contrast)
                for column in ['logFC','PValue','BH_FDR','moderated_t','df_total']:
                    np.testing.assert_allclose(result[column],native[column],atol=1e-12,rtol=0,equal_nan=True)
                self.assertEqual(result.test_status.tolist(),native.test_status.tolist())
                self.assertTrue(result.FDR_method.eq(method).all())
                self.assertEqual(result.significant.tolist(),(result.FDR<.05).tolist())
                self.assertTrue(pd.isna(result.loc['P999','FDR']))
        pvalue=self.table('fdrtool_pvalue')
        np.testing.assert_allclose(pvalue.FDR_PValue,pvalue.PValue,atol=1e-14,equal_nan=True)

    def test_q_values_match_independent_calls_and_diagnostics_export_same_fit(self):
        code='''library(fdrtool)
a<-commandArgs(TRUE);d<-read.csv(a[1]);d<-d[d$test_status=='tested',]
x<-if(a[3]=='pvalue')d$PValue else d$moderated_t
f<-fdrtool(x,statistic=a[3],plot=FALSE,verbose=FALSE)
write.csv(data.frame(FeatureID=d$FeatureID,p=f$pval,q=f$qval,lfdr=f$lfdr),a[2],row.names=FALSE)
'''
        script=self.root/'independent.R';script.write_text(code)
        for method,mode in [('fdrtool_pvalue','pvalue'),('fdrtool_normal','normal')]:
            for contrast in ['A_vs_B','A_vs_C','B_vs_C']:
                folder=self.root/method/'differential_expression/protein'/contrast
                reference=self.root/(method+'_'+contrast+'.csv')
                subprocess.run([self.rscript,'--vanilla',str(script),str(folder/'all_results.csv'),str(reference),mode],check=True,capture_output=True,text=True)
                expected=pd.read_csv(reference).set_index('FeatureID')
                result=self.table(method,contrast).loc[expected.index]
                for actual,ref in [('FDR','q'),('FDR_PValue','p'),('local_fdr','lfdr')]:
                    np.testing.assert_allclose(result[actual],expected[ref],atol=1e-12,rtol=0)
                diag=pd.read_csv(folder/'fdrtool_diagnostics.csv').set_index('FeatureID')
                np.testing.assert_allclose(diag.q_value,result.FDR,atol=1e-12)
                self.assertTrue((folder/'fdrtool_diagnostics.pdf').read_bytes().startswith(b'%PDF-'))
                self.assertTrue((folder/'fdrtool_parameters.tsv').is_file())
                self.assertTrue((folder/'fdrtool_warnings.txt').is_file())

    def test_empirical_normal_reports_changed_p_scale_and_method(self):
        from math import erfc,sqrt
        result=self.table('fdrtool_normal').dropna(subset=['PValue'])
        folder=self.root/'fdrtool_normal/differential_expression/protein/A_vs_B'
        parameters=pd.read_csv(folder/'fdrtool_parameters.csv').iloc[0]
        expected=np.array([erfc(abs(t)/(parameters['sd']*sqrt(2))) for t in result.moderated_t])
        np.testing.assert_allclose(result.FDR_PValue,expected,atol=1e-12,rtol=0)
        self.assertGreater(float(abs(result.FDR_PValue-result.PValue).max()),1e-3)
        report=(self.root/'fdrtool_normal/START_HERE.html').read_text()
        self.assertIn('fdrtool_normal',report)
        self.assertIn('not an adjustment of the original PValue',report)

    def test_invalid_method_rejected_before_analysis(self):
        self.assertEqual(differential_options()['fdr_method'],'BH')
        with self.assertRaisesRegex(ValueError,'fdr_method'):
            differential_options(dict(fdr_method='choose_smallest_q'))


if __name__=='__main__':unittest.main()
