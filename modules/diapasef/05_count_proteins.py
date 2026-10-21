#!/usr/bin/env python3
"""Export DIA protein/gene quantities, differential results and service benchmarks."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
# Also works when scripts are staged together by Nextflow.
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'pipeline_lib'))
from common import filter_report, load_report, read_method, read_samples
from quantification import canonical_group, load_matrix, annotate_groups, benchmark_detection, file_provenance
from quantification_analysis import analyze_matrix
from matrix_catalog import export_matrix_formats


def count_samples(samples, report, method, service_path=None):
    if set(report.Run) != set(samples.run):
        raise ValueError('Report does not contain exactly the requested runs')
    filtered=filter_report(report,method)
    filtered['group_key']=filtered['Protein.Group'].map(canonical_group)
    fragments=[c for c in filtered if c.startswith('Fr.') and c.endswith('.Quantity')]
    filtered['fragment_signal_exported']=filtered[fragments].gt(0).any(axis=1) if fragments else False
    inconsistent=filtered.groupby(['group_key','Run'])['PG.MaxLFQ'].nunique().gt(1)
    if inconsistent.any(): raise ValueError('Conflicting protein quantities within the same protein group/run')
    matrix=filtered.groupby(['group_key','Run'])['PG.MaxLFQ'].first().unstack('Run').reindex(columns=samples.run)
    matrix.index.name='Protein.Group'
    service=None
    if service_path:
        service,_,_=load_matrix(service_path,samples)
    rows,evidence=[],[]
    annotation=annotate_groups(filtered,matrix.index)
    for sample in samples.itertuples(index=False):
        current=filtered.loc[filtered.Run.eq(sample.run)]
        proteins=current.groupby('group_key',sort=False).agg(
            quantity=('PG.MaxLFQ','first'),distinct_peptides=('Stripped.Sequence','nunique'),
            distinct_precursors=('Precursor.Id','nunique'),best_precursor_q=('Q.Value','min'),
            best_run_protein_q=('PG.Q.Value','min'),precursors_with_fragment_signal=('fragment_signal_exported','sum'))
        for group,protein in proteins.iterrows():
            evidence.append(dict(sample_id=sample.sample_id,sample_type=sample.sample_type,run=sample.run,
                                 protein_group=group,Genes=annotation.at[group,'Genes'],**protein.to_dict()))
        new_count=len(proteins);old=int(service[sample.sample_id].gt(0).sum()) if service is not None else None
        rows.append(dict(sample_id=sample.sample_id,sample_type=sample.sample_type,run=sample.run,
            library_reference=sample.library_reference,detected_protein_groups=new_count,
            protein_groups_with_2_peptides=int(proteins.distinct_peptides.ge(2).sum()),
            detected_peptides=current['Stripped.Sequence'].nunique(),detected_precursors=current['Precursor.Id'].nunique(),
            service_protein_groups=old,count_difference=new_count-old if old is not None else None,
            percent_gain=100*(new_count/old-1) if old else None,beats_service_count=new_count>old if old is not None else None))
    columns=['sample_id','sample_type','run','protein_group','Genes','quantity','distinct_peptides','distinct_precursors',
             'best_precursor_q','best_run_protein_q','precursors_with_fragment_signal']
    return pd.DataFrame(rows),pd.DataFrame(evidence,columns=columns),matrix,filtered


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['samples','report','method','outdir']: parser.add_argument('--'+name,required=True)
    parser.add_argument('--service-matrix');parser.add_argument('--fasta')
    parser.add_argument('--differential',default='{}',help='JSON intensity-model settings')
    args=parser.parse_args();samples=read_samples(args.samples)
    counts,evidence,matrix,filtered=count_samples(samples,load_report(args.report),read_method(args.method),args.service_matrix)
    output=Path(args.outdir);output.mkdir(parents=True,exist_ok=True)
    counts.to_csv(output/'sample_counts.tsv',sep='\t',index=False);counts.to_csv(output/'sample_counts.csv',index=False)
    annotation=annotate_groups(filtered,matrix.index,args.fasta)
    annotation.join(matrix).to_csv(output/'protein_matrix_1pct.tsv',sep='\t')
    named=matrix.rename(columns=dict(zip(samples.run,samples.sample_id)))
    annotation.join(named).to_csv(output/'protein_matrix_sample_ids.tsv',sep='\t')
    for filename in ('protein_matrix_1pct.tsv', 'protein_matrix_sample_ids.tsv'):
        export_matrix_formats(output/filename)
    if not evidence.empty: evidence['Genes']=evidence.protein_group.map(annotation.Genes)
    evidence.to_csv(output/'protein_evidence.tsv',sep='\t',index=False)
    filtered.to_parquet(output/'filtered_precursors.parquet',index=False)
    samples.to_csv(output/'samples.tsv',sep='\t',index=False)
    service=mapping=None
    if args.service_matrix: service,_,mapping=load_matrix(args.service_matrix,samples)
    analyze_matrix(named,annotation,samples,output/'delivery',json.loads(args.differential),service,mapping)
    provenance={key:file_provenance(path) for key,path in vars(args).items() if key in ('samples','report','method','fasta','service_matrix') and path}
    (output/'delivery/quality_control/input_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print(f'Exported {len(counts)} samples and {len(matrix)} protein groups to {output}')


if __name__=='__main__': main()
