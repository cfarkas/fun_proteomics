#!/usr/bin/env python3
"""Run all-pairs DDA count models from delivered matrices, without raw searches."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'pipeline_lib'))
from result_reports import organize
from result_layout import publish_results
from quantification import file_provenance


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results',required=True,help='Existing final_results directory')
    parser.add_argument('--outdir',required=True,help='New empty run directory')
    parser.add_argument('--methods',nargs='+',choices=['edgeR','limma','edgeRun'],default=['edgeR','limma'])
    parser.add_argument('--features',nargs='+',choices=['protein','gene','peptide'],default=['protein','gene','peptide'])
    parser.add_argument('--control',help='Omit for all unordered condition pairs')
    parser.add_argument('--fdr',type=float,default=.05)
    parser.add_argument('--filter-low-counts',action='store_true')
    parser.add_argument('--rscript',default=os.environ.get('PROTEOMICS_RSCRIPT','Rscript'))
    args=parser.parse_args();source=Path(args.results).resolve();out=Path(args.outdir).resolve()
    if out.exists() and any(path.name != '.fun_proteomics' for path in out.iterdir()):
        parser.error('Use a new or empty output directory')
    if not 0<args.fdr<1:parser.error('FDR must be between zero and one')
    samples=pd.read_csv(source/'main_results/quality_control/sample_annotations.csv',dtype=str,keep_default_na=False)
    groups=[];provenance=[file_provenance(source/'main_results/quality_control/sample_annotations.csv')]
    for folder in sorted((source/'main_results/counts').glob('*/*/raw')):
        cohort,mode=folder.parts[-3:-1]
        if mode not in ('dda','ddaplus'):continue
        selected=samples.loc[samples.cohort.eq(cohort)]
        if selected.empty:raise ValueError('No sample metadata for '+cohort)
        if selected.duplicated(['condition','replicate']).any():raise ValueError('Replicate IDs must be independent within each condition')
        if args.control and args.control not in set(selected.condition):raise ValueError('Reference condition not found: '+args.control)
        dest=out/'analysis_details/reporting'/cohort/mode;dest.mkdir(parents=True,exist_ok=True)
        for feature in args.features:
            path=folder/(feature+'_counts.csv') if feature!='gene' else source/'optional_analyses/gene_counts/counts'/cohort/mode/'raw/gene_counts.csv'
            shutil.copy2(path,dest/path.name);provenance.append(file_provenance(path))
        selected.to_csv(dest/'samples.csv',index=False)
        groups.append(dict(cohort=cohort,mode='DDA+' if mode=='ddaplus' else 'DDA',directory=str(dest)))
    if not groups:raise ValueError('No delivered DDA matrices found')
    final=out/'final_results';final.mkdir(exist_ok=True)
    samples.to_csv(final/'sample_annotations.csv',index=False)
    categories={'protein':'proteins','gene':'genes','peptide':'peptides'}
    job=dict(groups=groups,outdir=str(final),feature_categories={key:categories[key] for key in args.features},
             differential=dict(enabled=True,methods=args.methods,comparisons='control' if args.control else 'all_pairwise',
                 control=args.control or '',fdr=args.fdr,filter_low_counts=args.filter_low_counts,workers=1),
             optional_analyses={'gene_counts':'gene' in args.features})
    config=final/'analysis_settings.json';config.write_text(json.dumps(job,indent=2)+'\n')
    env=os.environ.copy();env['R_LIBS']=os.pathsep.join(filter(None,[str(ROOT/'.runtime/R-library'),env.get('R_LIBS')]))
    subprocess.run([args.rscript,str(ROOT/'pipeline_lib/differential.R'),str(config)],env=env,check=True)
    organize(final)
    (final/'FINAL_RESULTS_COMPLETE').write_text('Count models completed from existing matrices.\n')
    (final/'input_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    publish_results(final)
    print('Results: '+str(final/'START_HERE.html'))


if __name__=='__main__':main()
