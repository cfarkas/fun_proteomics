#!/usr/bin/env python3
"""Import service/DIA/FragPipe quantities, benchmark detection and test all pairs."""
import argparse
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'pipeline_lib'))
from quantification import load_matrix, read_table, file_provenance, annotate_groups
from quantification_analysis import analyze_matrix


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--matrix',required=True)
    parser.add_argument('--samples',required=True,help='CSV/TSV with sample_id; differential testing also needs condition and replicate')
    parser.add_argument('--outdir',required=True,help='New or empty output directory')
    parser.add_argument('--quantity',choices=['auto','value','maxlfq','intensity','spectral_count'],default='auto')
    parser.add_argument('--measure',choices=['intensity','spectral_count'],default='intensity')
    parser.add_argument('--id-column')
    parser.add_argument('--fasta')
    parser.add_argument('--service-matrix')
    parser.add_argument('--service-quantity',choices=['auto','value','maxlfq','intensity','spectral_count'],default='auto')
    parser.add_argument('--service-samples',help='Optional separate sample mapping with sample_id and service_sample')
    parser.add_argument('--differential',action='store_true')
    parser.add_argument('--control',help='Omit to compare every unordered condition pair')
    parser.add_argument('--fdr',type=float,default=.05)
    parser.add_argument('--min-replicates',type=int,default=2)
    parser.add_argument('--normalization',choices=['none','median'],default='none')
    parser.add_argument('--rscript',default=None)
    args=parser.parse_args();output=Path(args.outdir)
    if output.exists() and any(output.iterdir()): parser.error('Output directory must be new or empty; existing results are preserved')
    samples=read_table(args.samples)
    if 'sample_type' in samples and samples.sample_type.nunique()>1: parser.error('Analyze sample types separately')
    if 'cohort' in samples and samples.cohort.nunique()>1: parser.error('Analyze cohorts separately')
    matrix,annotations,mapping=load_matrix(args.matrix,samples,args.quantity,args.id_column)
    if args.fasta:
        report=annotations.reset_index();annotations=annotate_groups(report,matrix.index,args.fasta)
    measure='spectral_count' if args.quantity=='spectral_count' or mapping.quantity.eq('spectral_count').any() else args.measure
    service=service_map=None
    if args.service_matrix:
        service_samples=read_table(args.service_samples) if args.service_samples else samples
        service,_,service_map=load_matrix(args.service_matrix,service_samples,args.service_quantity)
    options=dict(enabled=args.differential,comparisons='control' if args.control else 'all_pairwise',control=args.control,
                 fdr=args.fdr,min_replicates=args.min_replicates,normalization=args.normalization)
    if args.rscript: options['rscript']=args.rscript
    analyze_matrix(matrix,annotations,samples,output,options,service,service_map,measure)
    mapping.to_csv(output/'quality_control/input_sample_mapping.csv',index=False)
    provenance={key:file_provenance(path) for key,path in vars(args).items() if key in ('matrix','samples','service_matrix','service_samples','fasta') and path}
    (output/'quality_control/input_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print('Results: '+str((output/'START_HERE.html').resolve()))


if __name__=='__main__': main()
