"""Validate optional analyses and export sequence-supported peptide/isoform evidence."""
import collections
import csv
import json
import math
import re
from pathlib import Path

ROLES = ('control', 'drug_a', 'drug_b', 'combination')


def migrate_options(plan):
    """Drop retired discovery settings when loading configurations from older versions."""
    options = plan.setdefault('optional_analyses', {})
    # Preserve existing studies' gene analyses; new studies explicitly default off.
    options.setdefault('gene_counts', plan.get('proteomics', {}).get('analysis') != 'immunopeptidomics')
    options.pop('open_search', None)
    settings = options.get('statistics', {})
    if isinstance(settings, dict) and isinstance(settings.get('features'), list) and 'mass_shift' in settings['features']:
        settings['features'] = [kind for kind in settings['features'] if kind != 'mass_shift']
        if not settings['features']:
            settings['enabled'] = False
            settings.pop('features')
    return plan


def validate_options(plan):
    migrate_options(plan)
    from optional_statistics import validate_statistics
    validate_statistics(plan)
    options = plan.get('optional_analyses', {})
    for key in ('hla_enriched', 'peptide_isoforms', 'gene_counts'):
        if key in options and not isinstance(options[key], bool):
            raise ValueError(key + ' must be true or false.')
    hla = plan.get('proteomics', {}).get('analysis') == 'immunopeptidomics'
    if hla and options.get('gene_counts'):
        raise ValueError('Gene-count comparisons require whole-cell proteomics; HLA evidence supports peptide tests only.')
    if 'hla_enriched' in options and options['hla_enriched'] != hla:
        raise ValueError('HLA-enriched samples must use the immunopeptidomics workflow; select No for global samples.')
    drugs = options.get('drugs', {})
    interaction = drugs.get('interaction', {})
    if 'applied' in drugs and not isinstance(drugs['applied'], bool):
        raise ValueError('Drugs applied must be true or false.')
    if 'enabled' in interaction and not isinstance(interaction['enabled'], bool):
        raise ValueError('Drug-combination interaction enabled must be true or false.')
    if 'matched_design' in interaction and not isinstance(interaction['matched_design'], bool):
        raise ValueError('Dose and treatment-time matching annotation must be true or false.')
    if not drugs.get('applied'):
        if drugs.get('study_type') or interaction.get('enabled'):
            raise ValueError('Drug analyses require Drugs applied = Yes.')
        return
    if drugs.get('study_type') not in ('synergy', 'single'):
        raise ValueError('Choose whether you are testing synergy or single drugs.')
    if not interaction.get('enabled'):
        return
    if drugs['study_type'] != 'synergy':
        raise ValueError('Drug-combination interaction requires the synergy study option.')
    methods = interaction.get('methods', ['edgeR'])
    if (not isinstance(methods, list) or not methods
            or any(method not in ('edgeR', 'limma', 'monte_carlo') for method in methods)
            or len(set(methods)) != len(methods)):
        raise ValueError('Choose one or more different interaction methods: edgeR, limma-voom or Monte Carlo.')
    fdr = interaction.get('fdr', 0.05)
    if isinstance(fdr, bool) or not isinstance(fdr, (int, float)) or not math.isfinite(fdr) or not 0 < fdr <= 1:
        raise ValueError('Interaction FDR must be greater than 0 and at most 1.')
    if 'monte_carlo' in methods:
        monte_carlo = interaction.get('monte_carlo', {})
        if not isinstance(monte_carlo, dict):
            raise ValueError('Monte Carlo settings must be an object.')
        for key, default, low, high in (('simulations', 9999, 99, 1000000),
                                      ('seed', 1, 0, 2147483646), ('workers', 1, 1, 8)):
            value = monte_carlo.get(key, default)
            if isinstance(value, bool) or not isinstance(value, int) or not low <= value <= high:
                raise ValueError(f'Monte Carlo {key} must be an integer between {low} and {high}.')
    conditions = [interaction.get(role) for role in ROLES]
    if not all(conditions) or len(set(conditions)) != 4:
        raise ValueError('Select four different conditions: control, drug A, drug B and combination.')
    samples = [s for s in plan['samples'] if s['cohort'] == interaction.get('cohort')]
    if not samples:
        raise ValueError('Select an existing cohort for the drug-combination analysis.')
    for condition in conditions:
        replicates = [s['replicate'] for s in samples if s['condition'] == condition]
        if len(replicates) < 2 or len(replicates) != len(set(replicates)):
            raise ValueError('Drug-combination interaction needs at least two independently numbered biological replicates per selected condition: ' + condition)


def read_csv(path):
    with Path(path).open(newline='') as handle:
        return list(csv.DictReader(handle))


def write_csv(path, fields, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def export_isoforms(plan, report, cohort, out):
    """Match every observed peptide against all FASTA targets, treating I/L as equal."""
    from pipeline import fasta_records, accession, digest
    out = Path(out)
    samples = [s for s in plan['samples'] if s['cohort'] == cohort]
    ids = [s['sample_id'] for s in samples]
    rows = read_csv(Path(report)/'peptide_spectral_counts.csv')
    peptides = {}
    for row in rows:
        peptide = row['Peptide']
        if peptide in peptides or not re.fullmatch('[A-Z]+', peptide):
            raise ValueError('Duplicate or invalid peptide sequence: ' + peptide)
        values = {sid: int(row[sid]) for sid in ids}
        if any(v < 0 for v in values.values()):
            raise ValueError('Negative peptide spectral count')
        peptides[peptide] = values
    # A short seed avoids a peptide-by-protein Cartesian scan, without extra dependencies.
    seed_length = min(5, min(map(len, peptides), default=5))
    seeds = collections.defaultdict(list)
    for peptide in peptides:
        normalized = peptide.replace('I', 'L')
        seeds[normalized[:seed_length]].append((peptide, normalized))
    records, families = {}, collections.defaultdict(set)
    for header, sequence in fasta_records(Path(plan['database']['fasta']).read_text()):
        identifier = header.split()[0]
        if identifier.startswith(('rev_', 'contam_')):
            continue
        acc = accession(header)
        if acc in records:
            raise ValueError('Ambiguous duplicate FASTA accession for isoform mapping: ' + acc)
        gene = re.search(r'\bGN=([^\s]+)', header)
        family = re.sub(r'-\d+$', '', acc)
        records[acc] = {'sequence': sequence, 'gene': gene[1] if gene else '', 'family': family}
        families[family].add(acc)
    matches = collections.defaultdict(dict)
    for acc, record in records.items():
        sequence = record['sequence'].replace('I', 'L')
        for pos in range(len(sequence)-seed_length+1):
            for peptide, normalized in seeds.get(sequence[pos:pos+seed_length], ()):
                if sequence.startswith(normalized, pos):
                    matches[peptide].setdefault(acc, []).append(pos)
    evidence, summaries = [], []
    unique_counts = collections.defaultdict(collections.Counter)
    unique_peptides = collections.defaultdict(set)
    isoform_peptides = collections.defaultdict(set)
    covered = collections.defaultdict(set)
    for peptide, values in sorted(peptides.items()):
        mapped = matches[peptide]
        unique = next(iter(mapped)) if len(mapped) == 1 else ''
        detected = sum(v > 0 for v in values.values())
        discriminating = bool(detected and unique and len(families[records[unique]['family']]) > 1)
        if unique and detected:
            unique_counts[unique].update(values)
            unique_peptides[unique].add(peptide)
            if discriminating:
                isoform_peptides[unique].add(peptide)
        for acc, positions in mapped.items():
            if not detected:
                continue
            for pos in positions:
                covered[acc].update(range(pos, pos+len(peptide)))
        evidence.append({'Peptide': peptide, 'Length': len(peptide), 'Detected_samples': detected,
            'Matching_accessions': '; '.join(sorted(mapped)), 'Matching_proteins': len(mapped),
            'Unique_accession': unique, 'Isoform_discriminating': discriminating,
            'Mapping_status': 'unique' if unique else 'shared' if mapped else 'unmapped',
            'Positions_1_based': '; '.join(acc+':'+','.join(str(p+1) for p in positions) for acc, positions in sorted(mapped.items())),
            **values})
        for condition in dict.fromkeys(s['condition'] for s in samples):
            selected = [values[s['sample_id']] for s in samples if s['condition'] == condition]
            summaries.append({'Peptide': peptide, 'condition': condition,
                'mean_spectral_count': sum(selected)/len(selected),
                'detected_samples': sum(v > 0 for v in selected), 'total_samples': len(selected)})
    protein_rows = []
    for acc in sorted(covered):
        record = records[acc]
        protein_rows.append({'Accession': acc, 'Gene': record['gene'], 'Family': record['family'],
            'Family_entries_in_FASTA': len(families[record['family']]),
            'Unique_peptides': len(unique_peptides[acc]),
            'Isoform_discriminating_peptides': len(isoform_peptides[acc]),
            'Sequence_coverage_fraction': len(covered[acc])/len(record['sequence']),
            **{sid: unique_counts[acc][sid] for sid in ids}})
    write_csv(out/'peptide_evidence.csv', ['Peptide', 'Length', 'Detected_samples', 'Matching_accessions', 'Matching_proteins',
        'Unique_accession', 'Isoform_discriminating', 'Mapping_status', 'Positions_1_based']+ids, evidence)
    write_csv(out/'peptide_condition_summary.csv', ['Peptide', 'condition', 'mean_spectral_count', 'detected_samples', 'total_samples'], summaries)
    write_csv(out/'isoform_evidence.csv', ['Accession', 'Gene', 'Family', 'Family_entries_in_FASTA', 'Unique_peptides',
        'Isoform_discriminating_peptides', 'Sequence_coverage_fraction']+ids, protein_rows)
    settings = {'database_sha256': digest(plan['database']['fasta']), 'target_entries': len(records),
        'peptide_counts_sha256': digest(Path(report)/'peptide_spectral_counts.csv'),
        'families_with_alternatives': sum(len(v) > 1 for v in families.values()), 'il_equivalent': True,
        'peptides': len(peptides), 'unique_peptides': sum(bool(r['Unique_accession']) and r['Detected_samples']>0 for r in evidence),
        'isoform_discriminating_peptides': sum(r['Isoform_discriminating'] for r in evidence)}
    (out/'settings.json').write_text(json.dumps(settings, indent=2)+'\n')
    (out/'README.md').write_text('''# Peptide-level and isoform evidence

`peptide_evidence.csv` remaps observed sequences to every target in the searched FASTA, with I/L treated as indistinguishable. All matching proteins and sequence positions are retained. Shared peptides are never assigned to one isoform. `peptide_condition_summary.csv` describes detection and mean spectral counts per condition; zeros do not establish biological absence.

`isoform_evidence.csv` reports observed sequence coverage (including shared peptides), unique peptide counts and per-sample spectral counts from unique peptides only. An isoform-discriminating peptide must map to exactly one target across the entire FASTA, whose accession family contains alternatives (UniProt accession and -N suffix convention). Uniqueness is relative to this database, not all possible sequences. Supply an isoform-inclusive FASTA in the Database step before searching; canonical-only databases cannot establish isoform specificity. Non-UniProt identifiers remain protein evidence without an inferred isoform family.

These are sequence-supported candidates, not a separate isoform-identification FDR claim or reconstructed full proteoforms. Optional analyses → Statistics can compare peptide spectral counts and counts from isoform-discriminating peptides between conditions. The resulting BH FDR concerns condition differences, not isoform identification or switching. Peptide differential testing is also available through DEPs. HLA peptides describe source-sequence evidence, not source-protein abundance or isoform expression. Review the identification FDR of the parent search. No spectrum is searched by this reporting step.
''')
