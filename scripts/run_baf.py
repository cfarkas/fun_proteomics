#!/usr/bin/env python3
"""Resume verified BAF conversion, FragPipe DDA search, LFQ and service benchmarking."""
import argparse
import fcntl
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from types import SimpleNamespace
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'pipeline_lib'))
sys.path.insert(0, str(ROOT/'src'))
from pipeline import digest, read, save, prepare_tools, search, search_properties
from search_checkpoint import identity, candidate_matches
from quantification import load_matrix, annotate_groups, file_provenance
from quantification_analysis import analyze_matrix, write_dashboard
from fun_proteomics.baf_study import biological_samples, converter_python


@contextmanager
def directory(path):
    old = Path.cwd(); Path(path).mkdir(parents=True, exist_ok=True); os.chdir(path)
    try: yield
    finally: os.chdir(old)


def quarantine(path):
    if path.exists():
        dest = path.parent/'previous_attempts'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')/path.name
        dest.parent.mkdir(parents=True, exist_ok=True); path.rename(dest)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', required=True)
    p.add_argument('--outdir', required=True)
    args = p.parse_args(); cfg = read(args.config); output = Path(args.outdir).resolve()
    from fun_proteomics.study import samples_for, check
    from fun_proteomics.assignments import service_samples_for
    check(cfg, ROOT)
    frame = samples_for(cfg, ROOT); bio = biological_samples(frame, cfg.get('subject_column'))
    if output.exists() and any(output.iterdir()) and not (output/'inputs/identity.json').exists() and any(p.name!='.fun_proteomics' for p in output.iterdir()):
        raise ValueError('Output folder already contains unrelated files; choose a new project')
    output.mkdir(parents=True, exist_ok=True)
    lock = (output/'.baf.lock').open('a')
    try:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        lock.close()
        raise ValueError('Another BAF process is using this output folder')
    state = output/'inputs/identity.json'
    paths = ['fasta','service_matrix','service_samples','samples','msfragger_zip','ionquant_zip','diatracer_zip','fragpipe_zip']
    wanted = dict(config=cfg, inputs={key:file_provenance(cfg[key]) for key in paths if cfg.get(key)},
                  workflow_sha256=digest(ROOT/'assets/fragpipe_baf_lfq.workflow'))
    if state.is_file() and read(state)!=wanted: raise ValueError('Saved BAF inputs/settings changed; use a new project')
    save(state,wanted)
    template = output/'inputs/search.workflow'
    if not template.is_file() or digest(template)!=wanted['workflow_sha256']:
        shutil.copy2(ROOT/'assets/fragpipe_baf_lfq.workflow',template)
    for suffix, sep in [('csv',','),('tsv','\t')]:
        bio.to_csv(output/('inputs/biological_samples.'+suffix),sep=sep,index=False)
        frame.to_csv(output/('inputs/acquisitions.'+suffix),sep=sep,index=False)
    converted = output/'intermediates/mzML'; converted.mkdir(parents=True,exist_ok=True)
    python = converter_python(cfg, ROOT)
    def convert_one(row):
        target=converted/(row['sample_id']+'.mzML')
        print('Converting/verifying '+row['sample_id']+'; log: '+str(target.with_suffix('.log')),flush=True)
        with (converted/(row['sample_id']+'.log')).open('a') as log:
            subprocess.run([python,str(ROOT/'pipeline_lib/baf.py'),'--raw',row['raw_path'],'--out',str(target)],
                           stdout=log,stderr=subprocess.STDOUT,check=True,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1'))
        audit=read(target.with_suffix('.conversion.json'))
        print('Verified MS1+MS2 conversion: '+row['sample_id'],flush=True)
        return {**row,'raw_path':str(target)}, dict(sample_id=row['sample_id'],**audit)
    # Bound disk traffic/memory independently of search CPU allocation.
    with ThreadPoolExecutor(max_workers=min(8,int(cfg['cpus']),max(1,int(cfg['memory_gib'])//4))) as pool:
        results=list(pool.map(convert_one,frame.to_dict('records')))
    samples,audit=map(list,zip(*results)); save(output/'inputs/conversion_audit.json',audit)
    import pandas as pd
    for suffix,sep in [('csv',','),('tsv','\t')]:
        pd.DataFrame([{k:r[k] for k in ('sample_id','ms1_spectra','ms2_spectra','last_rt_seconds')} for r in audit]).to_csv(output/('inputs/conversion_summary.'+suffix),sep=sep,index=False)
    for item in wanted['inputs'].values():
        if file_provenance(item['path']) != item:
            raise ValueError('Input changed during conversion: '+item['path'])
    if digest(template)!=wanted['workflow_sha256']:
        raise ValueError('Staged workflow changed during conversion')
    plan=dict(samples=samples,raw_audit=audit,tools={k:cfg[k] for k in ('msfragger_zip','ionquant_zip','diatracer_zip','fragpipe_zip') if cfg.get(k)},
        database={'fasta':cfg['fasta']},modes=['DDA'],proteomics={'acquisition':'DDA-BAF','analysis':'global-label-free'})
    planpath=output/'inputs/search_plan.json';save(planpath,plan)
    toolroot=output/'intermediates/tools';toolset=toolroot/'toolset'
    toolmark=toolroot/'integrity.json'
    valid=toolmark.is_file() and all((toolset/k).is_file() and digest(toolset/k)==v for k,v in read(toolmark).items())
    if not valid:
        quarantine(toolset)
        with directory(toolroot): prepare_tools(SimpleNamespace(plan=str(planpath),versions=str(ROOT/'assets/tool_versions.json')))
        save(toolmark,{str(f.relative_to(toolset)):digest(f) for f in toolset.rglob('*') if f.is_file()})
    database=output/'inputs/search.fasta'
    if not database.is_file() or digest(database)!=digest(cfg['fasta']): shutil.copy2(cfg['fasta'],database)
    searchdir=output/'intermediates/search'; searchdir.mkdir(parents=True,exist_ok=True)
    arguments=SimpleNamespace(plan=str(planpath), cohort=str(frame.cohort.iloc[0]),mode='DDA',tools=str(toolset),
        template=str(template),database=str(database),cpus=cfg['cpus'],memory=cfg['memory_gib'],dry_run=False)
    wanted_search=identity(plan,arguments.cohort,'DDA',search_properties(arguments,plan),digest(cfg['fasta']),read(toolset/'software_versions.json'))
    try: reusable=candidate_matches(searchdir/'run',wanted_search)
    except (ValueError,OSError,KeyError): reusable=False
    if not reusable:
        quarantine(searchdir/'run');quarantine(searchdir/'raw')
        print('Starting FragPipe DDA search and MS1 LFQ; log: '+str(searchdir/'run/pipeline.log'),flush=True)
        with directory(searchdir): search(arguments)
    else: print('Verified completed FragPipe search and LFQ; reusing.',flush=True)
    final=output/'final_results'
    # Publication is atomic. A failed report cannot masquerade as a completed run.
    from fun_proteomics.continuation import verify_delivery
    if not verify_delivery(output):
        print('Verified final matrices, comparisons and benchmark.');return
    temporary=output/'report_build';quarantine(temporary)
    pipeline_samples=bio.drop(columns=['service_sample'],errors='ignore')
    quantity = 'maxlfq' if len(bio) > 1 else 'intensity'
    matrix,annotations,mapping=load_matrix(searchdir/'run/combined_protein.tsv',pipeline_samples,quantity)
    annotations=annotate_groups(annotations.reset_index(),matrix.index,cfg['fasta'])
    service=service_map=None
    if cfg.get('service_matrix'):
        service,_,service_map=load_matrix(cfg['service_matrix'],service_samples_for(cfg,bio,ROOT),sheet=cfg.get('service_sheet'))
    options=dict(enabled=cfg['differential'],comparisons='all_pairwise',fdr=cfg['fdr'],min_replicates=cfg['min_replicates'],normalization=cfg['normalization'])
    if cfg.get('subject_column'): options['subject_column']=cfg['subject_column']
    analyze_matrix(matrix,annotations,bio,temporary,options,service,service_map)
    from baf_reporting import export_spectral_counts
    export_spectral_counts(searchdir/'run',pipeline_samples,annotations,temporary)
    with (temporary/'README.md').open('a') as handle:
        handle.write('\nRaw BAF: protein_intensities and gene_intensities contain MS1 LFQ abundance. protein_spectral_counts, peptide_spectral_counts and gene_spectral_counts contain direct spectral evidence, aggregated across technical injections. Gene spectral counts sum only unambiguous peptide-to-gene assignments; exclusions are audited. All-pairs tests use LFQ intensities, not spectral counts. Positive LFQ can include match-between-runs quantification without direct MS2 evidence.\n')
    qc=temporary/'quality_control';mapping.to_csv(qc/'pipeline_sample_mapping.csv',index=False)
    save(qc/'raw_search_provenance.json',dict(inputs=wanted,conversion_audit=audit,
        search_checkpoint=read(searchdir/'run/search_checkpoint.json'),technical_aggregation='FragPipe experiment grouping', quantity=quantity,
        protein_fdr=.01,psm_fdr=.01,mbr=len(bio)>1,msbooster=False))
    if service_map is not None: save(qc/'service_import_audit.json',service_map.attrs)
    from baf_reporting import write_execution_report
    write_execution_report(temporary,frame,bio,audit,read(qc/'raw_search_provenance.json'))
    write_dashboard(temporary)
    (temporary/'FINAL_RESULTS_COMPLETE').write_text('complete\n')
    files=sorted(f for f in temporary.rglob('*') if f.is_file() and f.name!='OUTPUT_SHA256SUMS')
    (temporary/'OUTPUT_SHA256SUMS').write_text(''.join(digest(f)+'  '+f.relative_to(temporary).as_posix()+'\n' for f in files))
    quarantine(final);temporary.rename(final)
    print('Results: '+str(final/'START_HERE.html'),flush=True)


if __name__=='__main__': main()
