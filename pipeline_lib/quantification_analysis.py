"""Post-search reports shared by DIA-NN and imported service quantification."""
import html
import json
import os
from pathlib import Path
import shutil
import subprocess
from urllib.parse import quote
import pandas as pd
from quantification import gene_matrix, write_benchmark
from matrix_catalog import build_matrix_catalog, export_matrix_formats


from quantification_config import differential_options, validate_design_rows


def validate_design(samples, options):
    validate_design_rows(samples.to_dict('records'), options)


def write_dashboard(output, title='Proteomics results'):
    output = Path(output)
    for folder in ('quality_control','differential_expression','benchmark'):
        for path in sorted((output/folder).rglob('*.csv')): export_matrix_formats(path)
    comparisons = []
    for path in sorted((output/'differential_expression').rglob('results_index.csv')):
        frame = pd.read_csv(path, keep_default_na=False)
        parts = path.parent.relative_to(output/'differential_expression').parts
        for row in frame.to_dict('records'):
            row.update(cohort=parts[0] if len(parts)>1 else '', feature=parts[-1])
            for field in ('full_results', 'significant_results'):
                row[field] = (path.parent/row[field]).relative_to(output).as_posix()
            comparisons.append(row)
    sections = []
    if (output/'RAW_EXECUTION.html').is_file():
        sections.append('<h2>Raw search execution</h2><p><a href="RAW_EXECUTION.html">Complete pipeline, sample design, detection totals, comparisons and service benchmark</a></p>')
    if comparisons:
        pd.DataFrame(comparisons).to_csv(output/'results_index.csv', index=False)
        export_matrix_formats(output/'results_index.csv')
        body = ''
        for row in comparisons:
            body += '<tr>' + ''.join('<td>'+html.escape(str(row[key]))+'</td>' for key in
                ('cohort', 'feature', 'baseline', 'comparison_condition', 'tested_features', 'significant_features'))
            body += ''.join('<td><a href="'+quote(row[key])+'">'+label+'</a></td>' for key,label in
                [('full_results','All results'),('significant_results','Significant')]) + '</tr>'
        sections.append('<h2>Differential abundance</h2><p>Positive logFC means comparison minus reference. Each row has its own BH FDR family. '
            'Index: <a href="results_index.csv">CSV</a> · <a href="results_index.tsv">TSV</a>.</p><table><tr><th>Cohort</th><th>Feature</th>'
            '<th>Reference</th><th>Comparison</th><th>Tested</th><th>Significant</th><th>Full table</th><th>Subset</th></tr>'+body+'</table>')
    for path in sorted((output/'benchmark').rglob('detection_summary.csv')):
        frame = pd.read_csv(path)
        union = frame.loc[frame.sample_id.eq('__union__')]
        if not union.empty:
            row = union.iloc[0]
            folder = path.parent.relative_to(output).as_posix()
            audit_path=path.parent/'group_completeness_audit.csv'
            if audit_path.is_file():
                audit=pd.read_csv(audit_path).set_index('source')
                sections.append('<h2>All reported positive groups</h2><p>Before excluding incomplete group identifiers: pipeline <b>'+str(int(audit.loc['pipeline','observed_groups']))+
                    '</b>, service <b>'+str(int(audit.loc['service','observed_groups']))+'</b>. Incomplete positive groups excluded from exact matching: pipeline '+
                    str(int(audit.loc['pipeline','observed_incomplete_groups']))+', service '+str(int(audit.loc['service','observed_incomplete_groups']))+
                    '.</p><p><a href="'+quote(folder+'/reported_detection_summary.csv')+'">All reported detection counts per sample</a> · <a href="'+quote(folder+'/group_completeness_audit.csv')+'">Completeness audit</a> · <a href="'+quote(folder+'/accession_overlap_summary.csv')+'">Supporting accession overlap (not independent protein IDs)</a></p>')

            sections.append('<h2>Service detection benchmark</h2><p>Complete protein groups eligible for exact matching: pipeline <b>'+str(int(row.pipeline_detected))+
                '</b>, service <b>'+str(int(row.service_detected))+'</b>; shared '+str(int(row.shared))+', pipeline-only '+str(int(row.pipeline_only))+
                ', service-only '+str(int(row.service_only))+'. Net difference: <b>'+format(int(row.net_difference),'+d')+'</b>.</p><p><a href="'+
                quote(folder+'/detection_summary.csv')+'">Per-sample counts and differences</a> · <a href="'+quote(folder+'/protein_membership.csv')+
                '">Detected groups gained and lost</a> · <a href="'+quote(folder+'/README.md')+'">Benchmark interpretation</a></p>')
    for folder, label in [('quality_control','Quality control and provenance'),('differential_expression','Model settings and supporting tables'),('optional_analyses','Optional analyses')]:
        links = []
        for path in sorted((output/folder).rglob('*')):
            if path.is_file() and path.suffix in ('.csv','.tsv','.pdf','.png','.md','.json'):
                if folder == 'differential_expression' and path.name in ('all_results.csv','significant_results.csv'):
                    continue
                rel = path.relative_to(output).as_posix()
                links.append('<li><a href="'+quote(rel)+'">'+html.escape(rel)+'</a></li>')
        if links:
            sections.append('<details><summary>'+label+'</summary><ul>'+''.join(links)+'</ul></details>')
    plots = '<p><a href="REPORT.md">Detection report</a></p><img alt="Protein detection by sample" style="max-width:100%" src="protein_counts.png">' if (output/'protein_counts.png').exists() else ''
    (output/'START_HERE.html').write_text('<!doctype html><html><head><meta charset="utf-8"><title>'+html.escape(title)+
        '</title><style>body{font:16px system-ui;max-width:1200px;margin:3rem auto;padding:0 1rem;background:#fafafa;color:#172c35}a{color:#075c8c}li{margin:.5rem 0}h1,h2{color:#15544e}td,th{padding:.6rem;text-align:left;border-bottom:1px solid #ccd}details{margin:1.5rem 0}summary{cursor:pointer;font-weight:bold}.start{padding:1.4rem;background:#e3f2ef;border-radius:8px}</style></head><body><h1>'+html.escape(title)+
        '</h1><div class="start"><h2>Start with the matrices</h2><p><a href="MATRICES.html"><b>Protein and gene matrices — download CSV or TSV</b></a></p>'
        '<p>Primary quantities and transformed matrices are indexed together. Missing intensities stay blank. Gene quantities are sums of unambiguous protein-group assignments.</p></div>'
        '<p><a href="README.md">Methods and interpretation</a></p>'+''.join(sections)+plots+'</body></html>')
    build_matrix_catalog(output)


def analyze_matrix(matrix, annotations, samples, output, options=None, service=None, mapping=None, measure='intensity'):
    output=Path(output); options=differential_options(options)
    validate_design(samples, options)
    if options['enabled'] and measure!='intensity':
        raise ValueError('Intensity models cannot test spectral counts; use the DDA count workflow')
    if list(matrix.columns)!=list(samples.sample_id): raise ValueError('Sample order must match matrix columns')
    output.mkdir(parents=True,exist_ok=True); matrices=output/'matrices';matrices.mkdir(exist_ok=True)
    qc=output/'quality_control';qc.mkdir(exist_ok=True)
    samples.to_csv(qc/'samples.csv',index=False)
    protein=matrices/('protein_intensities.csv' if measure=='intensity' else 'protein_counts.csv')
    annotations.join(matrix).to_csv(protein)
    annotations.to_csv(qc/'protein_annotations.csv')
    genes=gene_matrix(matrix,annotations)
    gene=matrices/('gene_intensities.csv' if measure=='intensity' else 'gene_group_counts.csv')
    genes.to_csv(gene)
    annotations.loc[~annotations.gene_mapping_status.eq('single_gene')].to_csv(qc/'excluded_gene_assignments.csv')
    pd.DataFrame({'sample_id':matrix.columns,'detected_protein_groups':matrix.gt(0).sum().values,
                  'detected_unambiguous_genes':genes.gt(0).sum().values}).to_csv(qc/'sample_detection.csv',index=False)
    if service is not None: write_benchmark(matrix,service,output/'benchmark',mapping)
    if options['enabled']:
        rscript=shutil.which(options['rscript'])
        if not rscript: raise ValueError('Rscript not found. Install limma and jsonlite or set differential.rscript')
        for kind,path in [('protein',protein)]+([('gene',gene)] if options['gene_level'] else []):
            dest=output/'differential_expression'/kind;dest.mkdir(parents=True,exist_ok=True)
            job=dict(options,matrix=str(path.resolve()),samples=str((qc/'samples.csv').resolve()),outdir=str(dest.resolve()))
            config=dest/'settings.json';config.write_text(json.dumps(job,indent=2)+'\n')
            subprocess.run([rscript,str(Path(__file__).with_name('intensity_differential.R')),str(config)],check=True)
    (output/'README.md').write_text('# Quantification results\n\nOpen START_HERE.html. matrices/ contains annotated protein quantities and unambiguous gene rollups. quality_control/ records sample assignments, detection and excluded gene mappings. benchmark/ records matched service detection when supplied.\n\nGene intensity is the sum of positive protein-group quantities assigned to one gene; ambiguous or partial mappings are excluded. It is a protein-group rollup, not DIA-NN Gene.MaxLFQ, RNA expression, or independent gene-level identification FDR. Missing values remain blank. A gene sum can change with the number of observed groups.\n\nDifferential abundance uses limma-trend on log2 positive intensities, with no imputation and no normalization by default (DIA-NN PG.MaxLFQ is already normalized). Optional median centering is recorded. Each pair is fitted separately with at least the configured number of observed biological replicates per condition; each feature level and contrast has its own BH FDR. Full tables retain untestable features with blank p-values and an explicit status. The default model assumes independent replicates. When subject_column is configured, each contrast uses within-donor log2 differences and requires the configured number of complete donor pairs per feature. Technical injections must already be aggregated to one quantity per donor/condition. Paired limma-trend moderates the difference variances using mean abundance across complete pairs; no missing values are imputed. paired_samples.csv and paired_log2_differences.csv record the exact matching. No other batch adjustment is inferred.\n')
    (qc/'analysis_options.json').write_text(json.dumps(dict(options,measure=measure),indent=2)+'\n')
    write_dashboard(output)
