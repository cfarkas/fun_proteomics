"""Report HLA-enriched peptides from filtered PSMs without a protein-FDR gate."""
import collections
import csv
import json
import math
from pathlib import Path
import re

REFERENCE = 'https://doi.org/10.1016/j.mcpro.2025.101492'
FDR_DEFINITION = ('1% PSM FDR after MSBooster/Percolator; peptide, ion and protein '
                  'thresholds are 100%. No claim of 1% peptide or protein FDR. MBR is disabled.')


def enabled(plan):
    return plan.get('proteomics', {}).get('analysis') == 'immunopeptidomics'


def validate_workflow(properties):
    required = {'msfragger.num_enzyme_termini': '0', 'msfragger.search_enzyme_name_1': 'nonspecific',
                'msfragger.digest_min_length': '7', 'msfragger.digest_max_length': '25',
                'msfragger.output_report_topN': '1', 'ionquant.mbr': '0',
                'quantitation.run-label-free-quant': 'true', 'ionquant.run-ionquant': 'true',
                'msbooster.run-msbooster': 'true', 'percolator.run-percolator': 'true',
                'workflow.input.data-type.im-ms': 'true',
                'phi-report.filter': '--psm 0.01 --pep 1 --ion 1 --prot=1'}
    for key, value in required.items():
        if properties.get(key) != value:
            raise ValueError('Invalid immunopeptidomics workflow setting: ' + key)
    fixed = properties.get('msfragger.table.fix-mods', '')
    if not fixed or any(float(entry.strip().split(',')[0]) != 0 for entry in fixed.split(';')):
        raise ValueError('The immunopeptidomics preset must have no fixed modifications.')


def write(path, fields, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def organism_proteins(plan):
    """Read species from FASTA annotations, independently of protein inference."""
    organism = plan.get('organism')
    if not organism:
        return None
    allowed = set()
    annotated = False
    with Path(plan['database']['fasta']).open() as handle:
        for line in handle:
            if not line.startswith('>'):
                continue
            identifier = line[1:].split()[0]
            if identifier.startswith(('rev_', 'contam_')):
                continue
            match = re.search(r'\bOS=(.*?)(?=\s(?:OX|GN|PE|SV)=|$)', line.strip())
            annotated |= match is not None
            if match and match[1] == organism:
                allowed.add(identifier)
                parts = identifier.split('|')
                if len(parts) > 1:
                    allowed.add(parts[1])
    if not annotated:
        raise ValueError('HLA organism filtering requires OS= annotations in the FASTA. '
                         'Use an annotated FASTA or leave the organism filter blank.')
    if not allowed:
        raise ValueError('No FASTA targets match the HLA organism filter: ' + organism)
    return allowed


def export(plan, run, cohort, out):
    run, out = Path(run), Path(out)
    if not enabled(plan) or plan.get('modes') != ['DDA']:
        raise ValueError('HLA reporting requires the dedicated immunopeptidomics DDA workflow.')
    if not (run/'SEARCH_COMPLETE').is_file() or (run/'STUB_NOT_SCIENTIFIC_RESULTS').exists():
        raise ValueError('HLA reporting requires a completed scientific search.')
    properties = dict(line.split('=', 1) for line in (run/'fragpipe.workflow').read_text().splitlines()
                      if line and not line.startswith('#') and '=' in line)
    validate_workflow(properties)
    samples = [s for s in plan['samples'] if s['cohort'] == cohort]
    if not samples:
        raise ValueError('No HLA samples for cohort ' + cohort)
    ids = [s['sample_id'] for s in samples]
    hla_class = plan['immunopeptidomics']['hla_class']
    allowed = organism_proteins(plan)
    counts = collections.defaultdict(collections.Counter)
    proteins = collections.defaultdict(set)
    genes = collections.defaultdict(set)
    best = {}
    summaries, lengths = [], []
    for sample in samples:
        sid = sample['sample_id']
        seen, retained = {}, set()
        reasons = collections.Counter()
        length_counts = collections.Counter()
        with (run/sid/'psm.tsv').open(newline='') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            required = {'Spectrum', 'Peptide', 'Protein', 'Gene', 'Mapped Genes', 'Mapped Proteins',
                        'Assigned Modifications', 'Is Decoy', 'Is Contaminant'}
            if not required <= set(reader.fieldnames or []):
                raise ValueError(f'{sid}: missing HLA PSM columns {sorted(required-set(reader.fieldnames or []))}')
            for row in reader:
                if any(row[k].lower() in ('true', '1') for k in ('Is Decoy', 'Is Contaminant')) or row['Protein'].startswith(('rev_', 'contam_')):
                    reasons['decoy_or_contaminant_rows'] += 1
                    continue
                mappings = {p.strip() for value in (row['Protein'], row['Mapped Proteins'])
                            for p in re.split('[,;]', value) if p.strip() and p.strip() != '-'}
                mappings = {p for p in mappings if not p.startswith(('rev_', 'contam_'))}
                if allowed is not None and not mappings.intersection(allowed):
                    reasons['other_organism_rows'] += 1
                    continue
                peptide = row['Peptide']
                if not re.fullmatch('[A-Z]{7,25}', peptide):
                    raise ValueError(f'{sid}: unexpected peptide outside the 7–25 residue search: {peptide}')
                spectrum = (row.get('Spectrum File', ''), row['Spectrum'])
                if not row['Spectrum']:
                    raise ValueError(f'{sid}: missing spectrum identifier')
                assignment = (peptide, row['Assigned Modifications'], row.get('Charge', ''))
                if spectrum in seen and seen[spectrum] != assignment:
                    raise ValueError(f'{sid}: multiple assignments for one spectrum; expected a top-1 DDA search.')
                proteins[peptide].update(mappings)
                genes[peptide].update(g for value in (row['Gene'], row['Mapped Genes'])
                                      for g in re.split(r'[,;\s]+', value) if g and g not in ('-', 'NA'))
                score = row.get('Probability', row.get('PeptideProphet Probability', ''))
                if score:
                    probability = float(score)
                    if not math.isfinite(probability) or not 0 <= probability <= 1:
                        raise ValueError(f'{sid}: invalid PSM probability')
                    best[peptide] = max(best.get(peptide, 0), probability)
                if spectrum in seen:
                    reasons['duplicate_mapping_rows'] += 1
                    continue
                seen[spectrum] = assignment
                counts[peptide][sid] += 1
                length_counts[len(peptide)] += 1
                retained.add(peptide)
        summaries.append({'sample_id': sid, 'condition': sample['condition'],
                          'replicate': sample['replicate'], 'retained_psms': len(seen),
                          'unique_peptides': len(retained),
                          **{key: reasons[key] for key in ('decoy_or_contaminant_rows', 'other_organism_rows', 'duplicate_mapping_rows')}})
        lengths.extend({'sample_id': sid, 'length': length, 'psms': length_counts[length],
                        'unique_peptides': sum(len(p) == length for p in retained)} for length in range(7, 26))
    matrix, candidates = [], []
    for peptide in sorted(counts):
        mapped_genes = sorted(genes[peptide])
        values = {sid: counts[peptide][sid] for sid in ids}
        matrix.append({'Peptide': peptide, 'Protein': '; '.join(sorted(proteins[peptide])),
                       'Gene': mapped_genes[0] if mapped_genes else '',
                       'Mapped Genes': '; '.join(mapped_genes[1:]), **values})
        candidates.append({'Peptide': peptide, 'Length': len(peptide), 'Enriched HLA class': hla_class,
                           'Proteins': '; '.join(sorted(proteins[peptide])), 'Genes': '; '.join(mapped_genes),
                           'Best PSM probability': best.get(peptide, ''),
                           'Detected samples': sum(v > 0 for v in values.values()), **values})
    write(out.parent/'peptide_spectral_counts.csv', ['Peptide', 'Protein', 'Gene', 'Mapped Genes']+ids, matrix)
    write(out/'peptides.csv', ['Peptide', 'Length', 'Enriched HLA class', 'Proteins', 'Genes',
                             'Best PSM probability', 'Detected samples']+ids, candidates)
    write(out/'sample_summary.csv', ['sample_id', 'condition', 'replicate', 'retained_psms', 'unique_peptides',
                                   'decoy_or_contaminant_rows', 'other_organism_rows', 'duplicate_mapping_rows'], summaries)
    write(out/'length_distribution.csv', ['sample_id', 'length', 'psms', 'unique_peptides'], lengths)
    binding_min, binding_max = (8, 13) if hla_class == 'I' else (9, 22)
    (out/'binding_input.txt').write_text(''.join(p+'\n' for p in sorted(counts) if binding_min <= len(p) <= binding_max))
    (out/'settings.json').write_text(json.dumps({'hla_class': hla_class, 'search_length': [7, 25],
        'binding_input_length': [binding_min, binding_max], 'fdr': FDR_DEFINITION,
        'organism': plan.get('organism'), 'reference': REFERENCE, 'fragpipe_version': '24.0'}, indent=2)+'\n')
    (out/'README.md').write_text(f'''# HLA-enriched immunopeptidomics

HLA class {hla_class} is the sample enrichment selected by the user, not an allele or binding prediction.

- `peptides.csv`: sequences (7–25 residues), all source mappings, best PSM probability and spectral counts per sample. Shared mappings do not multiply counts. Modified forms are combined by stripped sequence.
- `sample_summary.csv`: retained PSMs, sequence counts and excluded/duplicate rows.
- `length_distribution.csv`: counts by peptide length and sample.
- `binding_input.txt`: unique unmodified sequences of length {binding_min}–{binding_max}, ready for external binding prediction. PTMs are not represented in this sequence list. No binding software is run here; no binding rank, allele assignment or immunogenicity is inferred.

{FDR_DEFINITION} PSM probability is a confidence score, not an adjusted p-value. Counts come directly from filtered `psm.tsv` files, without gating on combined protein groups. Decoys and contaminants are excluded; an optional organism filter uses FASTA OS annotations. Empty libraries and unsupported differential models are reported separately.

The DEPs settings apply to peptide spectral counts only in this workflow. Protein, gene and modification matrices are supporting source annotations; they do not measure source-protein expression. HLA enrichment and independent biological replicates are assumed.

Search preset: nonspecific digestion, 7–25 residues, charge 1–4, mass 200–5000 Da, precursor ±15 ppm, fragments 0.03 Da, no fixed modifications; variable oxidation (M), protein N-terminal acetylation and cysteinylation (C), at most three modifications per peptide. MSBooster/Percolator and IonQuant remain enabled; MBR is off.

Adapted from [Chen et al., Benchmarking Software for DDA-PASEF Immunopeptidomics]({REFERENCE}), using the pinned FragPipe 24.0 Nonspecific-HLA workflow. This is not an exact reproduction of the benchmark's software versions. Validation with representative HLA-enriched acquisitions is still required before scientific use.
''')
