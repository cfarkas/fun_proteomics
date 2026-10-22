"""Dependency-free preflight checks, before the DIA Python environment exists."""
import csv
import math
import os
from pathlib import Path
import re


def differential_options(options=None):
    if options is None: options={}
    allowed={'enabled','comparisons','control','fdr','min_replicates','normalization','gene_level','rscript','subject_column'}
    if not isinstance(options,dict) or set(options)-allowed: raise ValueError('Unknown intensity differential options')
    result=dict(enabled=False,comparisons='all_pairwise',fdr=.05,min_replicates=2,normalization='none',gene_level=True,rscript=os.environ.get('PROTEOMICS_RSCRIPT','Rscript'))
    result.update(options)
    for key in ('enabled','gene_level'):
        if not isinstance(result[key],bool): raise ValueError(key+' must be boolean')
    if result['comparisons'] not in ('all_pairwise','control'): raise ValueError('comparisons must be all_pairwise or control')
    if result['comparisons']=='control' and not result.get('control'): raise ValueError('A reference condition is required')
    if isinstance(result['fdr'],bool) or not 0<float(result['fdr'])<1: raise ValueError('FDR must be between zero and one')
    result['fdr']=float(result['fdr']);value=result['min_replicates']
    if isinstance(value,bool) or not str(value).isdigit() or int(value)<2: raise ValueError('min_replicates must be an integer >= 2')
    result['min_replicates']=int(value)
    if result['normalization'] not in ('none','median'): raise ValueError('normalization must be none or median')
    if result.get('subject_column') is not None:
        column=result['subject_column']
        if not isinstance(column,str) or not column.strip() or column != column.strip():
            raise ValueError('subject_column must name the donor/subject column in the sample table')
    return result


def validate_design_rows(samples,options):
    if not options['enabled']: return
    groups={};seen=set()
    for row in samples:
        if not {'sample_id','condition','replicate'}.issubset(row): raise ValueError('Differential analysis requires condition and replicate in the sample sheet')
        key=(str(row['condition']),str(row['replicate']))
        if not all(key) or key in seen: raise ValueError('Use nonempty conditions and unique biological replicate IDs within each condition')
        seen.add(key);groups[key[0]]=groups.get(key[0],0)+1
    if len(groups)<2 or min(groups.values())<options['min_replicates']: raise ValueError('Differential analysis requires at least two conditions with independent biological replicates')
    if options['comparisons']=='control' and options['control'] not in groups: raise ValueError('Reference condition not found in this sample type')
    column=options.get('subject_column')
    if column:
        from itertools import combinations
        subjects={}; seen_subjects=set()
        for row in samples:
            subject=str(row.get(column, '')).strip(); condition=str(row['condition'])
            if not subject: raise ValueError('Paired analysis needs a nonempty '+column+' for every sample')
            key=(condition,subject)
            if key in seen_subjects: raise ValueError('Repeated donor within a condition: aggregate technical injections before paired analysis')
            seen_subjects.add(key); subjects.setdefault(condition,set()).add(subject)
        pairs=list(combinations(groups,2)) if options['comparisons']=='all_pairwise' else [(options['control'],g) for g in groups if g!=options['control']]
        for a,b in pairs:
            if len(subjects[a] & subjects[b])<options['min_replicates']:
                raise ValueError('Paired comparison '+a+' vs '+b+' needs at least '+str(options['min_replicates'])+' matched donors')


def audit_service(path,samples):
    """Preflight CSV/TSV service tables with the same explicit column matching."""
    with Path(path).open('rb') as handle:
        if handle.read(4).startswith(b'PK'):
            raise ValueError('Native DIA service tables require CSV/TSV. Import Excel with analyze_quantification.py or export it as TSV.')
    with Path(path).open(encoding='utf-8-sig',newline='') as h:
        first=h.readline();h.seek(0)
        reader=csv.DictReader(h,delimiter='\t' if '\t' in first else ',');header=reader.fieldnames or []
        if not header or len(set(header))!=len(header) or any(not v for v in header): raise ValueError('Invalid service matrix header')
        identifier=next((c for c in ('Protein.Group','Protein','FeatureID','Protein ID','protein_group','group_key') if c in header),None)
        if identifier is None: raise ValueError('service matrix needs Protein.Group, Protein, FeatureID or Protein ID')
        suffix=next((s for s in (' MaxLFQ Intensity',' Intensity',' Spectral Count') if any(c.endswith(s) for c in header)), '')
        def normalize(value): return re.sub(r'(?i)\.(d|raw|mzml)$','',str(value).replace('\\','/').rsplit('/',1)[-1])
        lookup={}
        for column in header:
            if suffix and not column.endswith(suffix): continue
            if suffix==' Spectral Count' and column.endswith((' Unique Spectral Count',' Total Spectral Count')): continue
            lookup.setdefault(normalize(column[:-len(suffix)] if suffix else column),[]).append(column)
        selected=[]
        for sample in samples:
            aliases=[sample['service_sample']] if sample.get('service_sample') else [sample['sample_id'],Path(sample['raw_path']).stem]
            matches={c for alias in aliases for c in lookup.get(normalize(alias),[])}
            if len(matches)!=1: raise ValueError('service matrix sample match is missing or ambiguous: '+sample['sample_id'])
            column=matches.pop()
            if column in selected: raise ValueError('Duplicate service sample column mapping')
            selected.append(column)
        keys=set()
        for row in reader:
            raw=row.get(identifier,'') or ''
            if identifier=='Protein': raw+=';'+(row.get('Indistinguishable Proteins','') or '')
            parts=[p.strip() for p in re.split('[;,]',raw) if p.strip()]
            if any(re.match(r'(?i)^(contam|contaminant|rev|reverse|decoy|xxx)[_|-]',v) for p in parts for v in p.split('|')): continue
            key=';'.join(sorted({p.split('|')[1] if p.count('|')>=2 else p for p in parts}))
            if not key or key in keys: raise ValueError('Empty or duplicate service protein groups')
            keys.add(key)
            for c in selected:
                value=row.get(c)
                if value in ('','NA','NaN','nan','N/A','null'):continue
                try: numeric=float(value)
                except (ValueError,TypeError) as error: raise ValueError('Non-numeric service quantity') from error
                if not math.isfinite(numeric) or numeric<0: raise ValueError('Invalid service quantity')
