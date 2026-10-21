"""Consolidate published count matrices without combining search modes or cohorts."""
import argparse,csv,json,shutil,subprocess,os,re
from collections import defaultdict
from pathlib import Path

def read(path,delimiter=','):
    with Path(path).open() as f:
        r=csv.DictReader(f,delimiter=delimiter);return r.fieldnames or [],list(r)
def write(path,fields,rows):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
def peptide_counts(run,samples,proteins,out):
    fields,rows=read(Path(run)/'combined_peptide.tsv','\t')
    ids=[s['sample_id'] for s in samples]
    needed=['Peptide Sequence','Protein','Gene']+[sid+' Spectral Count' for sid in ids]
    if not set(needed)<=set(fields):raise ValueError('Missing peptide spectral-count columns')
    output=[];seen=set()
    for r in rows:
        if r['Protein'] not in proteins:continue
        key=r['Peptide Sequence']
        if key in seen:raise ValueError('Duplicate combined peptide sequence: '+key)
        seen.add(key)
        output.append({'Peptide':key,'Protein':r['Protein'],'Gene':r['Gene'],'Mapped Genes':r.get('Mapped Genes',''),**{sid:r[sid+' Spectral Count'] for sid in ids}})
    write(out,['Peptide','Protein','Gene','Mapped Genes']+ids,output)

def gene_counts(rows,ids):
    totals=defaultdict(lambda:{sid:0 for sid in ids});excluded=[]
    for row in rows:
        if 'Mapped Genes' not in row:raise ValueError('Peptide report lacks Mapped Genes. Regenerate peptide_spectral_counts.csv from combined_peptide.tsv before reporting gene counts.')
        genes={v for v in re.split(r'[;,\s]+',row.get('Gene','')+' '+row.get('Mapped Genes','')) if v and v not in {'NA','-'}}
        if len(genes)!=1:
            excluded.append({'Peptide':row['Peptide'],'Gene':row.get('Gene',''),'Mapped Genes':row.get('Mapped Genes',''),'Reason':'ambiguous genes' if genes else 'no gene annotation'})
            continue
        gene=next(iter(genes))
        for sid in ids:totals[gene][sid]+=int(row[sid])
    return [{'FeatureID':gene,'Gene':gene,**totals[gene]} for gene in sorted(totals)],excluded

def assemble(plan,reports,out,rscript='Rscript',rlib=None,search_summary=None):
    from optional_analyses import validate_options, export_isoforms
    from pipeline import digest
    validate_options(plan)
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    (out/'FINAL_RESULTS_COMPLETE').unlink(missing_ok=True)
    sample_ids=[s['sample_id'] for s in plan['samples']]
    hla=plan.get('proteomics',{}).get('analysis')=='immunopeptidomics'
    if len(sample_ids)!=len(set(sample_ids)):raise ValueError('Sample IDs must be globally unique')
    combined={k:[] for k in ('protein','peptide','modification','gene','modification_gene','modification_protein')};annotations={}
    manifest=[]
    for cohort in sorted({s['cohort'] for s in plan['samples']}):
        samples=[s for s in plan['samples'] if s['cohort']==cohort];ids=[s['sample_id'] for s in samples]
        for mode in plan['modes']:
            slug='ddaplus' if mode=='DDA+' else 'dda'
            name=f'tables_{cohort}_{slug}'
            matches=[Path(r) for r in reports if Path(r).name==name]
            if len(matches)!=1:raise ValueError('Expected one report directory: '+name)
            report=matches[0];dest=out.parent/'analysis_details/reporting'/cohort/slug;dest.mkdir(parents=True,exist_ok=True)
            if (report/'STUB_NOT_SCIENTIFIC_RESULTS').exists():raise ValueError('Final scientific results cannot use synthetic inputs')
            if hla:
                shutil.copytree(report/'immunopeptidomics',out/'immunopeptidomics'/cohort/slug,dirs_exist_ok=True)
            if plan.get('optional_analyses',{}).get('peptide_isoforms'):
                optional_dest=out/'optional_analyses/peptide_isoforms'/cohort/slug
                existing=report/'optional_analyses/peptide_isoforms'
                settings=json.loads((existing/'settings.json').read_text()) if (existing/'settings.json').exists() else {}
                if (settings.get('database_sha256')==digest(plan['database']['fasta']) and
                    settings.get('peptide_counts_sha256')==digest(report/'peptide_spectral_counts.csv')):
                    shutil.copytree(existing,optional_dest,dirs_exist_ok=True)
                else:export_isoforms(plan,report,cohort,optional_dest)
            ptm_report=dest/'gene_and_ptm_counts'
            shutil.copytree(report/'ptm_counts',ptm_report,dirs_exist_ok=True)
            if not (ptm_report/'protein_modification_peptide_counts.csv').is_file():
                raw_root=Path(plan['outdir'])/'analysis_details'/cohort/slug/'run'
                if not (raw_root/'SEARCH_COMPLETE').is_file():
                    raise ValueError('Rebuilding modification counts requires completed PSM files at '+str(raw_root))
                from ptm_counts import export as export_ptm_counts
                export_ptm_counts(plan,raw_root,cohort,ptm_report)
            sources={'protein':(report/'protein_spectral_counts.tsv','\t','Protein'),
                     'peptide':(report/'peptide_spectral_counts.csv',',','Peptide'),
                     'modification':(ptm_report/'modified_peptide_psm_counts.csv',',','Peptide'),
                     'modification_gene':(ptm_report/'gene_modification_peptide_counts.csv',',','Gene'),
                     'modification_protein':(ptm_report/'protein_modification_peptide_counts.csv',',','Protein')}
            for kind,(source,sep,idcol) in sources.items():
                fields,rows=read(source,sep)
                if not set(ids)<=set(fields):raise ValueError('Missing sample columns: '+str(source))
                ann=[x for x in fields if x not in ids];annotations[kind]=list(dict.fromkeys(annotations.get(kind,[])+ann))
                matrix=[];seen=set()
                for row in rows:
                    feature=row[idcol]+(' | '+row['Assigned Modifications'] if kind=='modification' else '')
                    if kind in ('modification_gene','modification_protein'):feature=row[idcol]+' | '+row['Modification']
                    if feature in seen:raise ValueError('Duplicate feature: '+feature)
                    seen.add(feature)
                    values={sid:int(row[sid]) for sid in ids}
                    if any(v<0 or str(v)!=str(row[sid]).strip() for sid,v in values.items()):raise ValueError('Expected nonnegative integer counts')
                    matrix.append({'FeatureID':feature,**{k:row[k] for k in ann},**values})
                    combined[kind].append({'Cohort':cohort,'Mode':mode,'FeatureID':feature,**{k:row[k] for k in ann},**{sid:values.get(sid,'') for sid in sample_ids}})
                write(dest/(kind+'_counts.csv'),['FeatureID']+ann+ids,matrix)
                if kind=='peptide':
                    genes,excluded=gene_counts(rows,ids)
                    annotations['gene']=['Gene']
                    write(dest/'gene_counts.csv',['FeatureID','Gene']+ids,genes)
                    write(dest/'excluded_gene_assignments.csv',['Peptide','Gene','Mapped Genes','Reason'],excluded)
                    combined['gene'].extend({'Cohort':cohort,'Mode':mode,**row,**{sid:row.get(sid,'') for sid in sample_ids}} for row in genes)
            for f in ('protein_intensity.tsv','protein_maxlfq_intensity.tsv'):
                if (report/f).exists():shutil.copy2(report/f,dest/f)
            write(dest/'samples.csv',['sample_id','condition','replicate','cohort'],[{k:s[k] for k in ('sample_id','condition','replicate','cohort')} for s in samples])
            manifest.append({'cohort':cohort,'mode':mode,'directory':str(dest.resolve())})
    for kind,rows in combined.items():
        write(out/(kind+'_counts.csv'),['Cohort','Mode','FeatureID']+annotations[kind]+sample_ids,rows)
    write(out/'sample_annotations.csv',['sample_id','condition','replicate','cohort'],[{k:s[k] for k in ('sample_id','condition','replicate','cohort')} for s in plan['samples']])
    options=dict(plan.get('differential',{'enabled':False}));options['workers']=min(8,options.get('workers',1),plan.get('resources',{}).get('cpus',1))
    job={'groups':manifest,'differential':options,'outdir':str(out.resolve()),
         'resources':plan.get('resources',{}),
         'analysis':plan.get('proteomics',{}).get('analysis','global-label-free'),
         'optional_analyses':plan.get('optional_analyses',{})}
    if hla:job['immunopeptidomics']=plan['immunopeptidomics']
    jobfile=out/'analysis_settings.json';jobfile.write_text(json.dumps(job,indent=2))
    env=os.environ.copy()
    if rlib:env['R_LIBS']=str(Path(rlib).resolve())
    subprocess.run([rscript,str(Path(__file__).with_name('differential.R')),str(jobfile)],check=True,env=env)
    from optional_statistics import run_statistics
    run_statistics(plan,out,rscript,env)
    if job['optional_analyses'].get('drugs',{}).get('interaction',{}).get('enabled'):
        subprocess.run([rscript,str(Path(__file__).with_name('drug_interaction.R')),str(jobfile)],check=True,env=env)
    try:
        from .result_reports import organize
    except ImportError:
        from result_reports import organize
    organize(out)
    job['outdir']='.'
    for group in job['groups']:group['directory']=os.path.relpath(group['directory'],out.resolve())
    jobfile.write_text(json.dumps(job,indent=2))
    comparisons=read(out/'comparisons.csv')[1] if (out/'comparisons.csv').exists() else []
    comparison_names=', '.join(dict.fromkeys(row['contrast'] for row in comparisons)) or 'Differential testing disabled.'
    methods=', '.join(options.get('methods',['edgeR','limma'])) if options.get('enabled') else 'Disabled'
    status=read(out/'analysis_summary.csv')[1] if (out/'analysis_summary.csv').exists() else []
    skipped=[f"{row['Cohort']} / {row['Mode']} / {row['FeatureType']}: {row['status']}" for row in status if row['status']!='tested']
    if (out/'method_notes.csv').exists():
        skipped.extend(f"{row['Cohort']} / {row['Mode']} / {row['FeatureType']} / {row['method']}: {row['status']}" for row in read(out/'method_notes.csv')[1])
    filtering=(f"Enabled: library-adjusted minimum {options.get('min_count',10)} and total minimum {options.get('min_total_count',15)} for proteins, genes and peptides; modification minima {options.get('modification_min_count',1)} and {options.get('modification_min_total_count',3)}." if options.get('filter_low_counts',False) else 'Disabled: every feature with a nonzero total is eligible for testing. No minimum-count or replicate-detection threshold is applied.')
    (out/'README.md').write_text(f"""# Results for this run

Open `START_HERE.html` in a browser for explanations and direct links to the results.

Methods: {methods}. Comparisons: {comparison_names}.

- `protein_counts.csv`, `gene_counts.csv`, `peptide_counts.csv`: raw spectral counts. Gene counts use unambiguous peptide-to-gene assignments.
- `modification_gene_counts.csv`, `modification_protein_counts.csv`: distinct modified peptide sequences per gene/protein **and modification type**. Repeated spectra, charge states and sites of the same sequence count once per type. Shared entity assignments are excluded and audited.
- `modification_counts.csv`: supporting modified-peptidoform PSM counts; this is a different counting unit.
- `differential/<method>/<pairwise|omnibus>/<cohort>/<mode>/<feature_type>/`: named contrasts for proteins, genes, peptides and optional modifications by gene/protein. For example, `genes/edgeR_CTRL_vs_Treatment.csv`. Each test has a matching CSV under `significant/`. `differential/results_index.csv` lists all tests, without combining pairwise and omnibus discoveries.
- `matrices/<cohort>/<mode>/`: normalized CPM and log2 CPM. `../zscore/<cohort>/<mode>/`: row z-scores of log2 CPM for all six count matrices; constant rows are zero, empty libraries are missing.
- `analysis_summary.csv`: tested-feature counts, sparsity, filtering and skipped analyses. `../analysis_details/reporting/`: complete filtering masks, model inputs, assignment audits and detailed counts.

In `REFERENCE_vs_COMPARISON`, positive logFC means higher in COMPARISON. Tables record both actual condition names and the direction. `FDR` is the only adjusted-p-value column and uses Benjamini–Hochberg within each method, feature type, contrast, cohort and search mode; significant means adjusted p < {options.get('fdr',0.05)}. Omnibus tables test any condition difference and contain condition averages instead of ambiguous fold-change columns. Every table includes `mean_count_GROUP`, `mean_normalized_CPM_GROUP`, `detected_samples_GROUP` and `total_samples_GROUP`. These describe evidence; use pairwise tests to identify significant group comparisons.

Sparse modification analysis uses TMMwsp normalization and the edgeR v4 small-count quasi-likelihood GLM (`legacy=FALSE`, robust moderation; edgeR ≥ 4.4.2). Low-count filtering: {filtering} All-zero features are excluded; model failures or empty libraries are reported separately. edgeRun (summation bound {options.get('edgeRun_upper',100000)}) and limma-voom are secondary comparisons for sparse counts. Zeros are not imputed and pseudocounts are not added to model input. See the [edgeR manual](https://bioconductor.org/packages/release/bioc/manuals/edgeR/man/edgeR.pdf).

Distinct-peptide counts describe detection diversity, not PTM occupancy or protein-abundance-adjusted modification changes. Gene-level tests use proteomics-derived counts, not RNA expression. Small counts and few replicates limit power. Low coverage does not fail statistical reporting: full result tables retain nonsignificant rows, and untestable rows carry an insufficient-evidence test_status with blank p-values/FDR and significant=false. Untestable rows are not evidence of no effect. Independent biological replicates and comparable sampling are assumed. DDA and DDA+ remain separate searches, never replicates. Blank cells denote other cohorts; zero denotes no retained evidence.

Skipped analyses: {'; '.join(skipped) if skipped else 'None.'}
""")
    if hla:
        note=('HLA-enriched immunopeptidomics: start with `immunopeptidomics/<cohort>/dda/peptides.csv`. '
              'DEPs tests peptide spectral counts only. Protein/gene/modification matrices are supporting source annotations, '
              'not source-protein expression. Identification uses 1% PSM FDR without a 1% peptide/protein FDR claim; '
              'this differs from the adjusted p-values of differential tests. MBR is off. '
              'Binding input lists are exported, but binding prediction is not run. See the HLA report README for settings and limitations.\n\n')
        text=(out/'README.md').read_text()
        (out/'README.md').write_text(text.replace('# Results for this run\n\n','# Results for this run\n\n'+note,1))
    if plan.get('optional_analyses',{}).get('peptide_isoforms') or job['optional_analyses'].get('drugs',{}).get('interaction',{}).get('enabled'):
        with (out/'README.md').open('a') as handle:
            handle.write('\nOptional analyses are under `optional_analyses/`: peptide/isoform sequence evidence and/or drug-combination interaction, according to the saved settings. Each analysis has its own README and settings; direct links are in `START_HERE.html`.\n')
    if job['optional_analyses'].get('statistics',{}).get('enabled'):
        with (out/'README.md').open('a') as handle:
            handle.write('\nStatistical comparisons for optional evidence are under `optional_analyses/statistics/START_HERE.html`. Their methods, reference/all-pairs comparisons, FDR and optional filtering are independent of DEPs and drug interaction. The statistics README explains peptide and isoform count definitions and exclusions.\n')
    (out/'FINAL_RESULTS_COMPLETE').write_text('Consolidated matrices and requested analyses completed.\n')
    from result_layout import publish_results
    publish_results(out,search_summary=search_summary)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',required=True);p.add_argument('--reports',nargs='+',required=True);p.add_argument('--out',default='final_results');p.add_argument('--rscript',default='Rscript');p.add_argument('--rlib');p.add_argument('--search-summary');a=p.parse_args()
    assemble(json.loads(Path(a.plan).read_text()),a.reports,a.out,a.rscript,a.rlib,a.search_summary)
