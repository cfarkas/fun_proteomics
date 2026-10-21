#!/usr/bin/env python3
"""Exercise the installed wheel outside its checkout; all data are synthetic."""
import csv
import json
import os
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import numpy as np
import pandas as pd
import fun_proteomics
from fun_proteomics.continuation import verify_delivery
from fun_proteomics.project import digest


def main():
    checkout = Path(__file__).resolve().parents[1]
    assert not Path(fun_proteomics.__file__).resolve().is_relative_to(checkout), 'Install the wheel before running this test'
    with tempfile.TemporaryDirectory(prefix='installed proteomics ') as temp:
        root=Path(temp); foreign=root/'unrelated directory'; foreign.mkdir()
        env=dict(os.environ,FUN_PROTEOMICS_CACHE=str(root/'workflow cache'),MPLCONFIGDIR=str(root/'matplotlib'))
        env.pop('PYTHONPATH',None)
        base=[sys.executable,'-m','fun_proteomics']
        def launch(*arguments, expected=0):
            result=subprocess.run([*base,*map(str,arguments)],cwd=foreign,env=env,capture_output=True,text=True,timeout=180)
            assert result.returncode==expected, result.stdout+'\n'+result.stderr
            return result.stdout+result.stderr
        assert fun_proteomics.__version__ in launch('--version')
        workflow=Path(launch('workflow-path').strip())
        assert workflow.is_relative_to(root/'workflow cache')
        assert (workflow/'pipeline_lib/results_dashboard.html').is_file()
        source=root/'old run';prepared=source/'prepared';prepared.mkdir(parents=True)
        rng=np.random.default_rng(19)
        x=rng.normal(15,.2,(60,9));x[:10,3:6]+=3;x[10:20,6:]+=2
        proteins=[f'P{i}' for i in range(60)]
        fasta=prepared/'search.fasta';fasta.write_text(''.join(f'>sp|{p}|NAME GN=G{i}\nMPEPTIDEK\n' for i,p in enumerate(proteins)))
        samples=[];rows=[]
        for index in range(9):
            sample=f'S{index+1}'; raw=root/'raw'/(sample+'.d');raw.mkdir(parents=True)
            with sqlite3.connect(raw/'analysis.tdf') as db:
                db.execute('CREATE TABLE Frames (MsMsType INTEGER)')
                db.executemany('INSERT INTO Frames VALUES (?)', [(0,),(9,)])
            (raw/'analysis.tdf_bin').write_bytes(b'synthetic acquisition')
            samples.append(dict(sample_id=sample,sample_type='tissue',condition='ABC'[index//3],replicate='1',raw_path=str(raw),run=sample,library_reference=index==0))
            for i,protein in enumerate(proteins):
                row={'Run':sample,'Protein.Group':protein,'Protein.Names':protein,'Genes':'','PG.MaxLFQ':2**x[i,index],
                     'Precursor.Id':f'PEPTIDE{i}2','Stripped.Sequence':f'PEPTIDE{i}'}
                row.update({key:.001 for key in ('Q.Value','Global.Q.Value','Lib.Q.Value','PG.Q.Value','Global.PG.Q.Value','Lib.PG.Q.Value')})
                rows.append(row)
        samples=pd.DataFrame(samples);samples.to_csv(prepared/'samples.tsv',sep='\t',index=False)
        task=source/'work/task';(task/'cohort').mkdir(parents=True)
        shutil.copy2(prepared/'samples.tsv',task/'samples.tsv')
        shutil.copy2(fasta,task/'search.fasta')
        shutil.copy2(workflow/'modules/diapasef/method.toml',task/'method.toml')
        report=task/'cohort/report.parquet';pd.DataFrame(rows).to_parquet(report)
        (task/'.exitcode').write_text('0\n')
        published=source/'analysis_details/search/tissue/cohort';published.mkdir(parents=True)
        shutil.copy2(report,published/'report.parquet')
        (published/'report.manifest.txt').write_text(json.dumps({'Main .parquet report':str(report)})+'\n')
        (published/'engine.log').write_text('Finished\n')
        original_hash=digest(report)
        service=pd.DataFrame(2**x[:50],index=proteins[:50],columns=samples.sample_id)
        service.loc['PROVIDER_ONLY']=100
        service.index.name='Protein.Group';service_path=root/'service quantities.csv';service.to_csv(service_path)
        result=root/'delivery'
        launch('--continue',source,'--outdir',result,'--independent-samples','--hash-raw','--service-matrix','tissue='+str(service_path))
        final=result/'final_results'
        assert verify_delivery(result)==[]
        matrix=pd.read_csv(final/'matrices/tissue/protein_intensities.csv')
        assert matrix.shape[0]==60 and matrix.Genes.notna().all()
        assert len(pd.read_csv(final/'matrices/tissue/gene_intensities.csv'))==60
        assert len(pd.read_csv(final/'results_index.csv'))==6
        for path in final.rglob('*.csv'):
            if path.with_suffix('.tsv').is_file():
                with path.open(newline='') as f: left=list(csv.reader(f))
                with path.with_suffix('.tsv').open(newline='') as f:right=list(csv.reader(f,delimiter='\t'))
                assert left==right,str(path)
        manifest_before=digest(final/'OUTPUT_SHA256SUMS')
        assert 'Completed outputs verified' in launch('--continue',source)
        assert 'Completed outputs verified' in launch('--continue',final)
        assert digest(final/'OUTPUT_SHA256SUMS')==manifest_before
        (final/'matrices/tissue/protein_intensities.csv').write_text('damaged delivery')
        launch('--continue',source)
        state=json.loads((source/'.fun_proteomics/run.json').read_text())
        rebuilt=Path(state['outdir'])
        assert rebuilt!=result and verify_delivery(rebuilt)==[]
        assert digest(report)==original_hash
        assert set(pd.read_csv(rebuilt/'inputs/samples.tsv',sep='\t',dtype=str).source_replicate)=={'1'}
        audit=json.loads((rebuilt/'audit/run.json').read_text())
        assert audit['raw_searches_rerun']==0 and audit['raw_searches_reused']==1
        changed=Path(samples.iloc[0].raw_path)/'analysis.tdf_bin';changed.write_bytes(b'changed acquisition')
        failure=launch('--continue',source,expected=2)
        assert 'Raw acquisition changed' in failure, failure
        # Standalone matrix import also records portable continuation state.
        sheet=root/'independent samples.csv';samples.assign(replicate=samples.sample_id).to_csv(sheet,index=False)
        imported=root/'quantified'
        launch('quantify','--matrix',service_path,'--samples',sheet,'--outdir',imported,'--fasta',fasta,'--differential')
        assert verify_delivery(imported)==[]
        assert 'Completed outputs verified' in launch('--continue',imported)
        (imported/'final_results/MATRICES.html').unlink()
        launch('--continue',imported)
        imported_state=json.loads((imported/'.fun_proteomics/run.json').read_text())
        assert verify_delivery(Path(imported_state['outdir']))==[]
        service_path.write_text(service_path.read_text()+'\n')
        assert 'Analysis input changed' in launch('--continue',imported,expected=2)
        # Count-model continuation uses its original matrices and sample metadata.
        counts_source=root/'delivered counts';raw_counts=counts_source/'main_results/counts/tissue/dda/raw'
        raw_counts.mkdir(parents=True)
        counts=pd.DataFrame(rng.poisson(30,(60,9)),columns=samples.sample_id)
        counts.insert(0,'FeatureID',proteins);counts.to_csv(raw_counts/'protein_counts.csv',index=False)
        metadata=counts_source/'main_results/quality_control/sample_annotations.csv';metadata.parent.mkdir(parents=True)
        samples.assign(cohort='tissue',replicate=samples.sample_id).to_csv(metadata,index=False)
        counted=root/'count analysis'
        launch('counts','--results',counts_source,'--outdir',counted,'--features','protein','--methods','limma')
        assert verify_delivery(counted)==[]
        count_index=pd.read_csv(counted/'final_results/results_index.csv')
        assert len(count_index)>=3
        assert 'Completed outputs verified' in launch('--continue',counted)
        (counted/'final_results/MATRICES.html').unlink()
        launch('--continue',counted)
        counted_state=json.loads((counted/'.fun_proteomics/run.json').read_text())
        assert verify_delivery(Path(counted_state['outdir']))==[]
        metadata.write_text(metadata.read_text()+'\n')
        assert 'Analysis input changed' in launch('--continue',counted,expected=2)
        print('PASS: installed wheel; DDA count modeling/continuation; separate working directory; DIA cache reuse; 6 protein/gene comparisons; gene fallback; CSV/TSV; service import; verified no-op; damaged delivery recovery; changed raw/matrix rejection.')


if __name__=='__main__':main()
