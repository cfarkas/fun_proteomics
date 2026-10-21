#!/usr/bin/env python3
"""Count assigned variable modifications in completed FragPipe PSM outputs."""
import argparse
import collections
import csv
import json
from pathlib import Path
import re

MOD = re.compile(r'([^,;]+?)\(([+-]?\d+(?:\.\d+)?)\)')

def modifications(text):
    found = MOD.findall(text)
    if MOD.sub('', text).strip(' ,;'):
        raise ValueError(f'Unrecognized Assigned Modifications: {text}')
    return [(site.strip(), float(mass)) for site, mass in found]

def kind(site, mass):
    residue = re.sub(r'\d', '', site).strip()
    if residue == 'C' and abs(mass-57.02146) < .002:
        return None  # Fixed alkylation is not counted as a variable biological PTM.
    if residue == 'M' and abs(mass-15.9949) < .002:
        return 'Oxidation (M)'
    if 'term' in residue.lower() and abs(mass-42.0106) < .002:
        return 'N-terminal acetylation'
    if residue in ('S','T','Y') and abs(mass-79.96633) < .002:
        return 'Phosphorylation (STY)'
    return f'{residue} ({mass:+.4f} Da)'

def write(path, fields, rows):
    with path.open('w', newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(rows)

def export(config, run, cohort, out):
    if not (run/'SEARCH_COMPLETE').is_file():
        raise ValueError('Search is not complete: missing SEARCH_COMPLETE. Retry after this cohort/mode finishes.')
    if (run/'STUB_NOT_SCIENTIFIC_RESULTS').exists():
        raise ValueError('Refusing synthetic search results')
    samples=[s for s in config['samples'] if s['cohort']==cohort]
    if not samples: raise ValueError(f'No samples in cohort {cohort}')
    # Match the organism filter via the inferred target protein groups. Include
    # indistinguishable members but do not multiply shared peptides across proteins.
    allowed=None
    if config.get('organism'):
        with (run/'combined_protein.tsv').open() as f:
            reader=csv.DictReader(f,delimiter='\t')
            if not {'Organism','Protein'} <= set(reader.fieldnames or []):
                raise ValueError('Cannot apply organism filter: missing combined protein annotations')
            allowed=set()
            for r in reader:
                if r['Organism']==config['organism'] and not r['Protein'].startswith(('rev_','contam_')):
                    allowed.add(r['Protein'])
                    allowed.update(x.strip() for x in r.get('Indistinguishable Proteins','').split(',') if x.strip())
    by_type=collections.defaultdict(collections.Counter)
    forms=collections.defaultdict(collections.Counter)
    summaries=[]
    gene_peptides=collections.defaultdict(lambda: collections.defaultdict(set))
    gene_forms=collections.defaultdict(lambda: collections.defaultdict(set))
    gene_all=collections.defaultdict(lambda: collections.defaultdict(set))
    gene_types=collections.defaultdict(lambda: collections.defaultdict(set))
    protein_types=collections.defaultdict(lambda: collections.defaultdict(set))
    protein_all=collections.defaultdict(lambda: collections.defaultdict(set))
    protein_modified=collections.defaultdict(lambda: collections.defaultdict(set))
    ambiguous_proteins=[]
    ambiguous=[]
    for s in samples:
        sid=s['sample_id']; seen=set(); modified=0; retained=0; ambiguous_count=0; missing_gene=0
        with (run/sid/'psm.tsv').open() as f:
            reader=csv.DictReader(f,delimiter='\t')
            required={'Spectrum','Peptide','Assigned Modifications','Protein','Is Decoy','Is Contaminant','Gene','Mapped Genes','Mapped Proteins'}
            if not required <= set(reader.fieldnames or []):
                raise ValueError(f'{sid}: missing PSM columns {required-set(reader.fieldnames or [])}')
            for r in reader:
                if any(r[k].lower() in ('true','1') for k in ('Is Decoy','Is Contaminant')) or r['Protein'].startswith(('rev_','contam_')): continue
                if allowed is not None and r['Protein'] not in allowed: continue
                mods=modifications(r['Assigned Modifications'])
                variable=sorted((site,mass) for site,mass in mods if kind(site,mass))
                assignment='; '.join(f'{site}({mass:.4f})' for site,mass in variable)
                # DDA+ can have multiple peptide IDs for a spectrum. Retain each
                # distinct spectrum/peptidoform, collapse duplicate protein mappings.
                key=(r['Spectrum'],r['Peptide'],assignment,r.get('Charge',''))
                if key in seen: continue
                seen.add(key); retained+=1
                genes={x for value in (r['Gene'],r['Mapped Genes']) for x in re.split(r'[,;\s]+',value.strip()) if x and x != '-'}
                if len(genes)==1:
                    gene=next(iter(genes))
                    gene_all[gene][sid].add(r['Peptide'])
                    if variable:
                        gene_peptides[gene][sid].add(r['Peptide'])
                        gene_forms[gene][sid].add((r['Peptide'],assignment))
                        for label in {kind(site,mass) for site,mass in variable}:
                            gene_types[(gene,label)][sid].add(r['Peptide'])
                elif variable:
                    if genes: ambiguous_count+=1
                    else: missing_gene+=1
                    ambiguous.append({'sample_id':sid,'Spectrum':r['Spectrum'],'Peptide':r['Peptide'],'Assigned Modifications':assignment,'Genes':'; '.join(sorted(genes)),'reason':'multiple genes' if genes else 'no gene annotation'})
                proteins={x.strip() for value in (r['Protein'],r['Mapped Proteins']) for x in re.split(r'[,;]',value) if x.strip() and x.strip()!='-'}
                if len(proteins)==1:
                    protein=next(iter(proteins))
                    protein_all[protein][sid].add(r['Peptide'])
                    if variable:
                        protein_modified[protein][sid].add(r['Peptide'])
                        for label in {kind(site,mass) for site,mass in variable}:
                            protein_types[(protein,label)][sid].add(r['Peptide'])
                elif variable:
                    ambiguous_proteins.append({'sample_id':sid,'Spectrum':r['Spectrum'],'Peptide':r['Peptide'],'Proteins':'; '.join(sorted(proteins)),'reason':'multiple proteins' if proteins else 'no protein annotation'})
                if not variable: continue
                modified+=1
                for label in {kind(site,mass) for site,mass in variable}: by_type[label][sid]+=1
                forms[(r['Peptide'],assignment)][sid]+=1
        summaries.append({**{k:s[k] for k in ('sample_id','condition','replicate')},'retained_psms':retained,'variable_modified_psms':modified,'ambiguous_gene_modified_psms':ambiguous_count,'missing_gene_modified_psms':missing_gene})
    out.mkdir(parents=True,exist_ok=True)
    ids=[s['sample_id'] for s in samples]
    # Include all identified unambiguous genes so zeros are interpretable.
    for filename,values in [('gene_modified_peptide_counts',gene_peptides),('gene_modified_peptidoform_counts',gene_forms),('gene_total_identified_peptide_counts',gene_all)]:
        write(out/(filename+'.csv'),['Gene']+ids,[{'Gene':g,**{sid:len(values[g][sid]) for sid in ids}} for g in sorted(gene_all)])
    write(out/'gene_modification_peptide_counts.csv',['Gene','Modification']+ids,[{'Gene':g,'Modification':label,**{sid:len(values[sid]) for sid in ids}} for (g,label),values in sorted(gene_types.items())])
    write(out/'protein_modification_peptide_counts.csv',['Protein','Modification']+ids,[{'Protein':protein,'Modification':label,**{sid:len(values[sid]) for sid in ids}} for (protein,label),values in sorted(protein_types.items())])
    for filename,values in [('protein_modified_peptide_counts',protein_modified),('protein_total_identified_peptide_counts',protein_all)]:
        write(out/(filename+'.csv'),['Protein']+ids,[{'Protein':protein,**{sid:len(values[protein][sid]) for sid in ids}} for protein in sorted(protein_all)])
    write(out/'ambiguous_protein_modified_peptides.csv',['sample_id','Spectrum','Peptide','Proteins','reason'],ambiguous_proteins)
    write(out/'ambiguous_modified_peptides.csv',['sample_id','Spectrum','Peptide','Assigned Modifications','Genes','reason'],ambiguous)
    write(out/'ptm_psm_counts.csv',['Modification']+ids,[{'Modification':label,**{sid:c[sid] for sid in ids}} for label,c in sorted(by_type.items())])
    write(out/'modified_peptide_psm_counts.csv',['Peptide','Assigned Modifications']+ids,[{'Peptide':key[0],'Assigned Modifications':key[1],**{sid:c[sid] for sid in ids}} for key,c in sorted(forms.items())])
    write(out/'ptm_sample_summary.csv',['sample_id','condition','replicate','retained_psms','variable_modified_psms','ambiguous_gene_modified_psms','missing_gene_modified_psms'],summaries)
    (out/'README.txt').write_text('PRIMARY TABLE: gene_modified_peptide_counts.csv = distinct stripped peptide sequences bearing at least one variable modification, per gene and sample. Different spectra, charge states, sites and modification types of the same sequence count once.\nOnly unambiguous gene annotations (Gene plus Mapped Genes has exactly one gene) enter gene tables; ambiguous or missing genes are listed separately. No razor assignment to a single gene is invented.\nGene total identified peptides supplies a detection denominator, not a PTM occupancy denominator. Peptidoform and gene-by-modification matrices are separate.\nCounts of distinct FDR-filtered spectrum/peptidoform/charge matches, not intensity or occupancy.\nOne PSM is counted once per modification type even if multiple residues carry that type.\nA PSM with two modification types contributes to both rows; rows must not be summed as unique PSMs.\nFixed C +57.02146 is excluded. Other unrecognized assigned modifications retain residue and mass labels.\nPeptide assignments are not confidence-filtered PTM sites; site localization scoring is disabled in the default search.\nShared peptides are counted once, not allocated across proteins. Species filtering uses retained combined protein assignments.\nZeros indicate no retained evidence, not absence of the modification. DDA and DDA+ must remain separate.\n')
    if config.get('proteomics',{}).get('analysis')=='phosphoproteomics':
        readme=out/'README.txt'
        readme.write_text(readme.read_text().replace('site localization scoring is disabled in the default search.', 'PTMProphet localization is enabled, but these peptide counts are not filtered by localization probability.'))
    with (out/'README.txt').open('a') as handle:
        handle.write('Protein-by-modification counts use distinct peptide sequences assigned to exactly one protein across Protein and Mapped Proteins. Shared protein assignments are excluded and audited; gene and protein counts therefore need not sum to the same total. Different variable modification types are separate rows.\n')
    return {'samples':len(samples),'modification_types':len(by_type),'modified_peptides':len(forms),'genes':len(gene_all),'output':str(out)}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',type=Path,required=True)
    p.add_argument('--run',type=Path,required=True,help='Completed cohort/mode run directory containing SEARCH_COMPLETE')
    p.add_argument('--cohort',default='study')
    p.add_argument('--outdir',type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(export(json.loads(a.config.read_text()),a.run,a.cohort,a.outdir),indent=2))
if __name__=='__main__': main()
