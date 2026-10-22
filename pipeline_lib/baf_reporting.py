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
