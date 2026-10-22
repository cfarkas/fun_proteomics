"""Scientific grouping, conversion corruption and service benchmark regression tests."""
import base64
from pathlib import Path
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'pipeline_lib')]
from baf import validate_mzml, source_files, nonempty_spectra
from quantification import write_benchmark
from fun_proteomics.baf_study import biological_samples
from fun_proteomics.study import normalize


def fixture(path):
    ns='http://psi.hupo.org/ms/mzml'
    def node(parent,tag,**attrs):return ET.SubElement(parent,'{'+ns+'}'+tag,{k:str(v) for k,v in attrs.items()})
    root=ET.Element('{'+ns+'}mzML');run=node(root,'run');spectra=node(run,'spectrumList',count=2)
    for level in (1,2):
        s=node(spectra,'spectrum',id='scan='+str(level),defaultArrayLength=2)
        node(s,'cvParam',accession='MS:1000511',value=level);node(s,'cvParam',accession='MS:1000127')
        node(s,'cvParam',accession='MS:1000016',value=level,unitAccession='UO:0000031')
        if level==2:
            precursor=node(s,'precursor',spectrumRef='scan=1')
            node(precursor,'cvParam',accession='MS:1000744',value=450.25)
        arrays=node(s,'binaryDataArrayList',count=2)
        for accession,values in [('MS:1000514',[450.,451.]),('MS:1000515',[10.,20.])]:
            a=node(arrays,'binaryDataArray');node(a,'cvParam',accession=accession);node(a,'cvParam',accession='MS:1000523')
            node(a,'binary').text=base64.b64encode(np.array(values,dtype='<f8').tobytes()).decode()
    ET.ElementTree(root).write(path,encoding='utf-8')


class BafTests(unittest.TestCase):
    def setUp(self):self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)
    def test_conversion_requires_ms1_ms2_arrays_and_parent(self):
        p=self.root/'run.mzML';fixture(p)
        self.assertEqual(validate_mzml(p),{'ms1_spectra':1,'ms2_spectra':1,'last_rt_seconds':120.})
        original=p.read_text()
        for broken in [original[:-20],original.replace('spectrumRef="scan=1"','spectrumRef="scan=999"'),original.replace('defaultArrayLength="2"','defaultArrayLength="3"'),original.replace('count="2"','count="9"'),original.replace('defaultArrayLength="2"','defaultArrayLength="0"')]:
            p.write_text(broken)
            with self.assertRaises((ValueError,ET.ParseError)):validate_mzml(p)
    def test_baf_preset_is_conventional_lcms_and_preserves_ms1_quantification(self):
        from pipeline import search_properties
        from types import SimpleNamespace
        settings=search_properties(SimpleNamespace(template=ROOT/'assets/fragpipe_baf_lfq.workflow',database=self.root/'db.fasta',cohort='study',mode='DDA',cpus=4,memory=16),{})
        self.assertEqual(settings['workflow.input.data-type.im-ms'],'false')
        self.assertEqual(settings['workflow.input.data-type.regular-ms'],'true')
        self.assertEqual(settings['msbooster.run-msbooster'],'false')
        self.assertEqual(settings['ionquant.run-ionquant'],'true')
    def test_baf_exports_direct_counts_and_excludes_ambiguous_gene_peptides(self):
        from baf_reporting import export_spectral_counts
        from quantification import load_matrix
        run=self.root/'run';run.mkdir()
        pd.DataFrame([
            {'Protein':'P1','Gene':'G1','Indistinguishable Proteins':'P2','S1 Spectral Count':2,'S2 Spectral Count':0},
            {'Protein':'P3','Gene':'G2','Indistinguishable Proteins':'','S1 Spectral Count':4,'S2 Spectral Count':5},
            {'Protein':'contam_P4','Gene':'X','Indistinguishable Proteins':'','S1 Spectral Count':99,'S2 Spectral Count':99},
        ]).to_csv(run/'combined_protein.tsv',sep='\t',index=False)
        pd.DataFrame([
            {'Peptide Sequence':'AAA','Protein':'P1','Gene':'G1','Mapped Genes':'','S1 Spectral Count':2,'S2 Spectral Count':0},
            {'Peptide Sequence':'BBB','Protein':'P3','Gene':'G2','Mapped Genes':'G3','S1 Spectral Count':4,'S2 Spectral Count':5},
            {'Peptide Sequence':'CCC','Protein':'contam_P4','Gene':'X','Mapped Genes':'','S1 Spectral Count':99,'S2 Spectral Count':99},
        ]).to_csv(run/'combined_peptide.tsv',sep='\t',index=False)
        samples=pd.DataFrame({'sample_id':['S1','S2']})
        _,annotations,_=load_matrix(run/'combined_protein.tsv',samples)
        output=self.root/'report';export_spectral_counts(run,samples,annotations,output)
        proteins=pd.read_csv(output/'matrices/protein_spectral_counts.csv')
        self.assertEqual(list(proteins['Protein.Group']),['P1;P2','P3'])
        self.assertEqual(list(proteins.S1),[2,4])
        genes=pd.read_csv(output/'matrices/gene_spectral_counts.csv')
        self.assertEqual(list(genes.Gene),['G1'])
        self.assertEqual(list(genes.S1),[2])
        self.assertEqual(len(pd.read_csv(output/'quality_control/excluded_peptide_gene_assignments.csv')),1)
    def test_vendor_dependencies_and_no_mixed_acquisition(self):
        raw=self.root/'sample.d';raw.mkdir()
        for name in ['analysis.baf','analysis.baf_xtr','Calibrator.ami','analysis.sqlite']:(raw/name).write_text('test')
        self.assertEqual([p.name for p in source_files(raw)],['Calibrator.ami','analysis.baf','analysis.baf_xtr'])
        (raw/'analysis.tdf').touch()
        with self.assertRaisesRegex(ValueError,'Mixed'):source_files(raw)
    def test_technical_injections_do_not_inflate_biological_n(self):
        frame=pd.DataFrame([dict(sample_id=f'{d}_{c}_{i}',biological_sample=f'{d}_{c}',condition=c,donor=d,cohort='study',raw_path=f'{d}_{c}_{i}.d') for d in ['D1','D2'] for c in ['A','B'] for i in range(3)])
        biological=biological_samples(frame,'donor')
        self.assertEqual(len(biological),4);self.assertEqual(biological.donor.nunique(),2)
        frame.loc[1,'condition']='B'
        with self.assertRaisesRegex(ValueError,'disagree'):biological_samples(frame,'donor')
    def test_empty_scans_are_audited_and_parent_evidence_is_not_invented(self):
        rows=[dict(Id=1,AcquisitionKey=1,Parent=None,LineMzId=10,LineIntensityId=11),
              dict(Id=2,AcquisitionKey=2,Parent=1,LineMzId=20,LineIntensityId=21)]
        keys={1:{'MsLevel':0},2:{'MsLevel':1}}
        kept,omitted=nonempty_spectra(rows,keys,lambda key:2 if key<20 else 0)
        self.assertEqual([r['Id'] for r in kept],[1]);self.assertEqual(omitted[0]['spectrum_id'],2)
        with self.assertRaisesRegex(ValueError,'empty MS1'):
            nonempty_spectra(rows,keys,lambda key:0 if key<20 else 2)
        with self.assertRaisesRegex(ValueError,'lengths differ'):
            nonempty_spectra(rows,keys,lambda key:key%2)
    def test_venv_interpreter_is_not_resolved_outside_venv(self):
        binary=self.root/'python';binary.symlink_to(sys.executable)
        config=normalize(dict(project=str(self.root/'study'),mode='baf',python=str(binary),subject_column='donor'))
        self.assertEqual(config['python'],str(binary))
    def test_benchmark_missing_groups_zero_service_and_csv_tsv(self):
        pipeline=pd.DataFrame({'S1':[1.,1.,1.],'S2':[1.,0.,0.]},index=['P1;P2','P3','P4'])
        service=pd.DataFrame({'S1':[1.,1.,1.],'S2':[0.,0.,0.]},index=['P1','P9','P10'])
        service.attrs['incomplete_groups']=['P10']
        out=self.root/'benchmark';summary=write_benchmark(pipeline,service,out)
        overall=summary.set_index('sample_id').loc['__union__']
        self.assertEqual(overall.pipeline_only,3);self.assertEqual(overall.service_only,2)
        accessions=pd.read_csv(out/'accession_overlap_summary.csv').set_index('sample_id')
        self.assertEqual(accessions.loc['__union__','shared'],1)
        self.assertTrue(pd.isna(summary.set_index('sample_id').loc['S2','percent_change']))
        audit=pd.read_csv(out/'group_completeness_audit.csv').set_index('source')
        self.assertEqual(audit.loc['service','observed_groups'],3);self.assertEqual(audit.loc['service','observed_incomplete_groups'],1)
        reported=pd.read_csv(out/'reported_detection_summary.csv').set_index('sample_id')
        self.assertEqual(reported.loc['__union__','service_reported_positive_groups'],3)
        self.assertEqual(reported.loc['__union__','net_reported_difference'],0)
        for csv in out.glob('*.csv'):
            pd.testing.assert_frame_equal(pd.read_csv(csv),pd.read_csv(csv.with_suffix('.tsv'),sep='\t'))


if __name__=='__main__':unittest.main()
