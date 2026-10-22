"""Provider-format imports and donor pairing use synthetic data only."""
import os
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'pipeline_lib'))
if not os.environ.get('FUN_PROTEOMICS_TEST_INSTALLED'): sys.path.insert(0,str(ROOT/'src'))
from quantification import load_matrix, gene_matrix, annotate_groups, benchmark_detection
from quantification_import import read_quantification
from quantification_analysis import analyze_matrix, differential_options, validate_design
from fun_proteomics import study, assignments


class ProviderTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)

    def workbook(self, name='report.xlsx', ambiguous=False):
        p=self.root/name
        table=pd.DataFrame({'Fasta headers':[
            'sp|P1|ONE_HUMAN Protein OS=Human GN=GENE1 PE=1 SV=1',
            'sp|P2|TWO_HUMAN Protein OS=Human GN=GENE2 PE=1 SV=1;tr|P9',
            'sp|P3|THREE_HUMAN Protein OS=Human GN=GENE3 PE=1 SV=1',
            'sp|P3|THREE_HUMAN Protein OS=Human GN=GENE3 PE=1 SV=1'],
            'LFQ intensity D1_A':[10,20,30,0], 'LFQ intensity D1_B':[40,50,60,0]})
        with pd.ExcelWriter(p) as writer:
            pd.DataFrame({'Report':['Provider report']}).to_excel(writer,sheet_name='Cover',index=False)
            table.to_excel(writer,sheet_name='Identified proteins',index=False,startrow=3)
            if ambiguous: table.to_excel(writer,sheet_name='Second quantities',index=False)
        return p

    def test_cover_page_header_offset_gene_annotations_and_empty_duplicate_audit(self):
        p=self.workbook();samples=pd.DataFrame({'sample_id':['D1_A','D1_B']})
        matrix,ann,mapping=load_matrix(p,samples)
        self.assertEqual(matrix.shape,(3,2));self.assertEqual(matrix.loc['P1','D1_B'],40)
        self.assertEqual(mapping.source_column.tolist(),['LFQ intensity D1_A','LFQ intensity D1_B'])
        self.assertEqual(mapping.attrs['import']['worksheet'],'Identified proteins')
        self.assertEqual(mapping.attrs['import']['header_row'],4)
        self.assertEqual(len(mapping.attrs['import']['removed_empty_duplicates']),1)
        self.assertEqual(ann.loc['P2','gene_mapping_status'],'partial')
        self.assertEqual(set(gene_matrix(matrix,ann).index),{'GENE1','GENE3'})
        fasta=self.root/'a.fasta';fasta.write_text('>sp|P2|TWO GN=GENE2\nMPEPTIDE\n')
        supplemented=annotate_groups(ann.reset_index(),matrix.index,fasta)
        self.assertEqual(supplemented.loc['P2','gene_mapping_status'],'partial')
        summary,_=benchmark_detection(matrix,matrix)
        self.assertEqual(summary.iloc[-1].pipeline_detected,2)

    def test_sheet_choice_and_duplicate_observed_groups_are_never_guessed(self):
        p=self.workbook(ambiguous=True)
        with self.assertRaisesRegex(ValueError,'Multiple protein'):read_quantification(p)
        self.assertEqual(len(read_quantification(p,sheet='Second quantities')),3)
        with self.assertRaisesRegex(ValueError,'worksheet not found'):read_quantification(p,sheet='missing')
        p=self.root/'duplicate.tsv'
        pd.DataFrame({'Fasta headers':['sp|P1|ONE GN=G SV=1']*2,'LFQ intensity A':[5,8]}).to_csv(p,sep='\t',index=False)
        with self.assertRaisesRegex(ValueError,'duplicate protein'):load_matrix(p,pd.DataFrame({'sample_id':['A']}))

    def test_maxquant_flags_columns_and_assignment_aliases(self):
        p=self.root/'proteinGroups.txt'
        pd.DataFrame({'Protein IDs':['P1','P2','P3','P4'],'Gene names':['G1','G2','G3','G4'],
                      'Reverse':['','+','',''],'Potential contaminant':['','','+',''],
                      'Only identified by site':['','','','+'],
                      'LFQ intensity A':[5,6,7,8],'LFQ intensity B':[9,10,11,12]}).to_csv(p,sep='\t',index=False)
        cfg=study.normalize({'project':str(self.root/'study'),'matrix':str(p),'differential':False})
        board=assignments.discover(cfg,ROOT)
        self.assertEqual([r['sample_id'] for r in board['rows']],['A','B'])
        board['rows'][0].update(sample_id='Renamed',condition='Control')
        cfg['assignment_board']=board
        matrix,ann,mapping=load_matrix(p,study.samples_for(cfg,ROOT))
        self.assertEqual(matrix.index.tolist(),['P1']);self.assertEqual(matrix.loc['P1','Renamed'],5)
        self.assertEqual(ann.loc['P1','Genes'],'G1');self.assertEqual(mapping.attrs['import']['excluded_provider_flag_rows'],3)

    def test_baf_detection_provides_raw_acquisition_cards(self):
        raw=self.root/'impact.d';raw.mkdir();(raw/'analysis.baf').touch()
        self.assertEqual(study.scan(self.root)[0]['raw_path'],str(raw))

    def test_paired_design_rejects_pseudoreplication_and_unmatched_subjects(self):
        options=differential_options({'enabled':True,'subject_column':'donor'})
        samples=pd.DataFrame({'sample_id':['A1','A2','B1','B2'],'condition':['A','A','B','B'],
                              'replicate':['1','2','1','2'],'donor':['D1','D2','D1','D2']})
        validate_design(samples,options)
        with self.assertRaisesRegex(ValueError,'technical injections'):validate_design(samples.assign(donor=['D1']*4),options)
        with self.assertRaisesRegex(ValueError,'matched donors'):validate_design(samples.assign(donor=['D1','D2','D3','D4']),options)
        with self.assertRaisesRegex(ValueError,'nonempty'):validate_design(samples.assign(donor=['D1','','D1','D2']),options)
        with self.assertRaisesRegex(ValueError,'either independent'):study.normalize({'project':str(self.root),'subject_column':'donor','independent_samples':True})

    def test_paired_limma_all_pairs_matches_within_donor_effects_and_missingness(self):
        rng=np.random.default_rng(123)
        samples=pd.DataFrame([dict(sample_id=f'{donor}_{condition}',condition=condition,replicate=donor,donor=donor)
                              for donor in ['D1','D2','D3','D4'] for condition in ['A','B','C','D']])
        offsets=np.repeat([0,4,8,12],4); effects=np.tile([0,2,-1,1],4)
        x=rng.normal(12,.2,(60,16))+offsets+effects
        x[-1,:]=np.nan;x[-1,[0,4,9,13]]=20 # two per condition, but no matched A/B donor
        x[-2,0]=np.nan # three complete A/B pairs
        matrix=pd.DataFrame(2**x,columns=samples.sample_id,index=pd.Index([f'P{i}' for i in range(60)],name='Protein.Group'))
        ann=pd.DataFrame({'Genes':[f'G{i}' for i in range(60)],'gene_mapping_status':'single_gene'},index=matrix.index)
        output=self.root/'results'
        analyze_matrix(matrix,ann,samples,output,{'enabled':True,'subject_column':'donor','rscript':os.environ.get('PROTEOMICS_RSCRIPT','Rscript')})
        index=pd.read_csv(output/'results_index.csv');self.assertEqual(len(index),12)
        result=pd.read_csv(output/'differential_expression/protein/A_vs_B/all_results.csv').set_index('FeatureID')
        delta=x[:,1::4]-x[:,0::4]
        np.testing.assert_allclose(result.loc['P0':'P58','logFC'],np.nanmean(delta[:-1],axis=1),rtol=1e-10)
        self.assertEqual(result.loc['P58','detected_pairs'],3)
        self.assertEqual(result.loc['P59','detected_pairs'],0);self.assertTrue(np.isnan(result.loc['P59','PValue']))
        self.assertTrue((result.loc['P0':'P57','FDR']<.05).all())
        self.assertEqual(set(index.method),{'paired limma-trend'})
        for stem in ['protein_intensities','gene_intensities']:
            pd.testing.assert_frame_equal(pd.read_csv(output/f'matrices/{stem}.csv'),pd.read_csv(output/f'matrices/{stem}.tsv',sep='\t'))
        self.assertTrue((output/'differential_expression/protein/A_vs_B/paired_samples.tsv').is_file())
        self.assertTrue((output/'differential_expression/protein/A_vs_B/paired_log2_differences.tsv').is_file())
        pairs=pd.read_csv(output/'differential_expression/protein/A_vs_B/paired_samples.csv')
        self.assertEqual(pairs.subject.tolist(),['D1','D2','D3','D4'])


if __name__=='__main__': unittest.main()
