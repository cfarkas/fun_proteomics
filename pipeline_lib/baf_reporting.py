"""Export direct spectral evidence beside the BAF LFQ analysis."""
from pathlib import Path
import numpy as np
import pandas as pd
from final_results import peptide_counts, gene_counts, read, write
from quantification import load_matrix, excluded_group


def export_spectral_counts(run, samples, annotations, output):
    run, output = Path(run), Path(output)
    matrices, qc = output/'matrices', output/'quality_control'
    matrices.mkdir(parents=True, exist_ok=True)
    qc.mkdir(parents=True, exist_ok=True)
    counts, _, _ = load_matrix(run/'combined_protein.tsv', samples, 'spectral_count')
    values = counts.to_numpy(dtype=float)
    if not np.isfinite(values).all() or not np.equal(values, np.floor(values)).all():
        raise ValueError('Direct spectral counts must be nonnegative integers')
    if not counts.index.equals(annotations.index):
        raise ValueError('LFQ and spectral-count protein groups disagree')
    counts = counts.astype('int64')
    annotations.join(counts).to_csv(matrices/'protein_spectral_counts.csv')
    _, groups = read(run/'combined_protein.tsv', '\t')
    proteins = {row['Protein'] for row in groups if not excluded_group(
        row['Protein']+';'+row.get('Indistinguishable Proteins', ''))}
    peptide_path = matrices/'peptide_spectral_counts.csv'
    peptide_counts(run, samples.to_dict('records'), proteins, peptide_path)
    _, peptides = read(peptide_path)
    identifiers = list(samples.sample_id)
    genes, excluded = gene_counts(peptides, identifiers)
    write(matrices/'gene_spectral_counts.csv', ['FeatureID','Gene',*identifiers], genes)
    write(qc/'excluded_peptide_gene_assignments.csv', ['Peptide','Gene','Mapped Genes','Reason'], excluded)
    pd.DataFrame({'sample_id': identifiers,
                  'direct_ms2_protein_groups': counts.gt(0).sum().values}).to_csv(
        qc/'direct_ms2_detection.csv', index=False)


def write_execution_report(output, acquisitions, samples, audit, provenance):
    """Summarize completed raw search evidence without pooling statistical tests."""
    import html
    output=Path(output)
    qc=output/'quality_control'
    annotations=pd.read_csv(qc/'protein_annotations.csv',keep_default_na=False)
    protein=pd.read_csv(output/'matrices/protein_intensities.csv')
    gene=pd.read_csv(output/'matrices/gene_intensities.csv')
    counts=pd.read_csv(output/'matrices/protein_spectral_counts.csv')
    ids=list(samples.sample_id)
    values=dict(acquisitions=len(acquisitions),biological_samples=len(samples),conditions=samples.condition.nunique(),
        source_spectra=sum(row['source_spectra'] for row in audit),
        ms1_spectra=sum(row['ms1_spectra'] for row in audit),ms2_spectra=sum(row['ms2_spectra'] for row in audit),
        omitted_empty_spectra=sum(len(row['omitted_empty_spectra']) for row in audit),
        identified_protein_groups=int(counts[ids].gt(0).any(axis=1).sum()),
        positive_lfq_protein_groups=int(protein[ids].gt(0).any(axis=1).sum()),
        positive_unambiguous_gene_rollups=int(gene[ids].gt(0).any(axis=1).sum()),
        protein_groups_with_gene_annotation=int(annotations.Genes.ne('').sum()))
    cfg=provenance['inputs']['config'];subject=cfg.get('subject_column')
    if subject: values['donors']=int(samples[subject].nunique())
    summary=pd.DataFrame([{'metric':key,'value':value} for key,value in values.items()])
    for suffix,sep in [('csv',','),('tsv','\t')]:summary.to_csv(qc/('execution_summary.'+suffix),sep=sep,index=False)
    def link(path,label):return '<a href="'+html.escape(path)+'">'+html.escape(label)+'</a>'
    displayed=summary.rename(columns={'metric':'Measure','value':'Total'})
    displayed['Measure']=displayed.Measure.str.replace('_',' ').str.capitalize()
    table=displayed.to_html(index=False,border=0,escape=True)
    sections=['<h2>Completed stages and evidence</h2><p>Vendor centroid MS1/MS2 conversion → MSFragger search → Percolator / ProteinProphet filtering → IonQuant MS1 quantification → gene annotation → requested comparisons and service benchmark. Each conversion and the completed search have verified checksums.</p>'+table]
    design='Paired contrasts by '+subject if subject else 'Independent biological samples'
    sections.append('<h2>Study design</h2><p>'+html.escape(design)+'. Technical injections share one FragPipe experiment and do not increase the biological sample count.</p>'+samples.groupby('condition',sort=True).size().rename('biological_samples').reset_index().to_html(index=False,border=0))
    index=output/'results_index.csv'
    if index.is_file():
        results=pd.read_csv(index)
        columns=['feature','baseline','comparison_condition','tested_features','significant_features']
        sections.append('<h2>All condition comparisons</h2><p>Positive logFC means comparison minus reference. BH FDR is controlled separately for each feature level and contrast. Untestable features retain blank p-values; they are not evidence of no difference.</p>'+results[columns].rename(columns={'feature':'Feature','baseline':'Reference','comparison_condition':'Comparison','tested_features':'Tested','significant_features':'Significant'}).to_html(index=False,border=0)+link('results_index.csv','CSV index')+' · '+link('results_index.tsv','TSV index'))
    benchmark=output/'benchmark/reported_detection_summary.csv'
    if benchmark.is_file():
        reported=pd.read_csv(benchmark);reported=reported[reported.sample_id.eq('__union__')]
        exact=pd.read_csv(output/'benchmark/detection_summary.csv');exact=exact[exact.sample_id.eq('__union__')]
        reported=reported.drop(columns='sample_id').rename(columns={'pipeline_reported_positive_groups':'Pipeline groups','service_reported_positive_groups':'Service groups','net_reported_difference':'Net change','percent_change':'Change (%)'})
        exact=exact.drop(columns='sample_id').rename(columns={key:key.replace('_',' ').capitalize() for key in exact.columns})
        sections.append('<h2>Service benchmark</h2><p>All positive groups as reported by each analysis:</p>'+reported.to_html(index=False,border=0)+'<p>Exact matching of complete groups:</p>'+exact.to_html(index=False,border=0)+'<p>'+link('benchmark/group_completeness_audit.csv','Incomplete-group audit')+' · '+link('benchmark/protein_membership.csv','Shared, additional and missing groups')+' · '+link('benchmark/README.md','Interpretation')+'</p><p>Incomplete service identifiers contribute to reported counts but cannot establish exact overlap. Databases, protein inference, FDR and transfer rules can differ; more positive quantities alone do not prove greater accuracy.</p>')
    sections.append('<h2>Methods and provenance</h2><p>Protein and PSM identification filtering: 1% FDR. Conventional LC–MS mode; MSBooster disabled. IonQuant match-between-runs uses 1% ion-transfer FDR when enabled. Positive LFQ can lack direct MS2 evidence in an individual sample. Missing intensities are not imputed. Gene LFQ is a sum of unambiguous protein-group assignments; gene spectral counts use unambiguous peptides.</p><p>Empty centroid scans are omitted and audited; no spectral peaks are invented. Original BAF components and calibration inputs are hashed and read from isolated copies.</p><p>'+link('quality_control/raw_search_provenance.json','Full input, conversion and search provenance')+' · '+link('../inputs/biological_samples.csv','Biological sample mapping')+' · '+link('../inputs/acquisitions.csv','Acquisition mapping')+' · '+link('../intermediates/search/run/pipeline.log','Complete search log')+' · '+link('../intermediates/search/run/fragpipe.workflow','Executed workflow')+' · '+link('../intermediates/search/run/fragger.params','Optimized search parameters')+'</p>')
    sections.append('<h2>Find every matrix</h2><p>'+link('MATRICES.html','Protein, peptide and gene matrices: CSV and TSV')+'. LFQ intensities and direct spectral counts are labeled separately. Statistical tests use LFQ, not spectral counts.</p>')
    (output/'RAW_EXECUTION.html').write_text('<!doctype html><meta charset="utf-8"><title>Raw BAF execution report</title><style>body{font:16px system-ui;color:#173d45;max-width:1250px;margin:3rem auto;padding:0 1rem}p{line-height:1.6}table{border-collapse:collapse}th,td{padding:.55rem;border-bottom:1px solid #ccd;text-align:left}a{color:#075c8c}h2{margin-top:2rem}</style><h1>Raw BAF execution report</h1><p>'+link('START_HERE.html','Result dashboard')+' · '+link('quality_control/execution_summary.csv','Execution summary CSV')+' · '+link('quality_control/execution_summary.tsv','Execution summary TSV')+'</p>'+''.join(sections))
