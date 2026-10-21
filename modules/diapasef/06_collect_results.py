#!/usr/bin/env python3
"""Validate all sample-type results and write one combined count table and report."""

import argparse
from pathlib import Path
import os
import shutil
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "pipeline_lib"))
from quantification_analysis import write_dashboard
from matrix_catalog import build_matrix_catalog
import re

import pandas as pd


def merge_counts(frames, expected):
    counts = pd.concat(frames, ignore_index=True)
    keys = ['sample_type', 'sample_id', 'run']
    if counts.duplicated(['sample_type', 'sample_id']).any():
        raise ValueError('Duplicate sample in output counts')
    if set(map(tuple, counts[keys].values)) != set(map(tuple, expected[keys].values)):
        raise ValueError('Final count table must contain every requested sample exactly once')
    def natural(value):
        return tuple(int(part) if part.isdigit() else part for part in re.split(r'(\d+)', value))
    order = sorted(range(len(counts)), key=lambda i: (counts.iloc[i].sample_type, natural(counts.iloc[i].sample_id)))
    return counts.iloc[order].reset_index(drop=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samples', required=True)
    parser.add_argument('--groups', nargs='+', required=True)
    args = parser.parse_args()
    types = {'sample_id': str, 'sample_type': str}
    frames = [pd.read_csv(Path(path) / 'sample_counts.tsv', sep='\t', dtype=types) for path in args.groups]
    expected = pd.read_csv(args.samples, sep='\t', dtype=types)
    counts = merge_counts(frames, expected)
    output = Path('final_results')
    output.mkdir(exist_ok=True)
    counts.to_csv(output / 'sample_counts.tsv', sep='\t', index=False)
    counts.to_csv(output / 'sample_counts.csv', index=False)
    expected.to_csv(output / 'samples.tsv', sep='\t', index=False)
    os.environ.setdefault('MPLCONFIGDIR', str(Path('.cache/matplotlib').resolve()))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    groups = list(counts.groupby('sample_type', sort=True))
    labels_by_group, label_rows = {}, []
    for sample_type, group in groups:
        ids = group.sample_id.astype(str).tolist()
        prefix = re.sub(r'[^_ /-]*$', '', os.path.commonprefix(ids))
        if len(prefix) < 16:
            prefix = ''
        labels = [value[len(prefix):] for value in ids]
        labels_by_group[sample_type] = labels
        label_rows.extend(dict(sample_type=sample_type, plot_label=label, sample_id=sid,
                               omitted_shared_prefix=prefix) for sid, label in zip(ids, labels))
    height = sum(4 + .065 * max(map(len, labels)) for labels in labels_by_group.values())
    fig, axes = plt.subplots(len(groups), 1, squeeze=False,
                            figsize=(max(10, max(len(g) for _, g in groups) * .3), height), constrained_layout=True)
    lines = ['# Protein quantification', '', f'Completed {len(counts)} requested samples across {len(groups)} sample types.', '',
             'Counts use positive protein quantities and 1% run, global and library confidence filters. '
             'Protein groups are not expanded into individual accessions. Missing quantities are not imputed.', '',
             '| Sample type | Samples | Median protein groups | Samples above service |',
             '| --- | ---: | ---: | ---: |']
    summaries = []
    for (sample_type, group), axis in zip(groups, axes[:, 0]):
        x = np.arange(len(group))
        benchmark = group.service_protein_groups.notna().all()
        if benchmark:
            axis.bar(x - .2, group.service_protein_groups, .4, label='Service', color='#aab2bd')
            axis.bar(x + .2, group.detected_protein_groups, .4, label='New: 1% filters', color='#207f8c')
            wins = int(group.count_difference.gt(0).sum())
        else:
            axis.bar(x, group.detected_protein_groups, color='#207f8c', label='1% filters')
            wins = None
        axis.set_xticks(x, labels_by_group[sample_type], rotation=90, fontsize=8)
        axis.set(title=sample_type, xlabel='Sample ID', ylabel='Quantified protein groups')
        axis.legend(frameon=False)
        median = float(group.detected_protein_groups.median())
        lines.append(f'| {sample_type} | {len(group)} | {median:g} | {str(wins) if wins is not None else "No benchmark supplied"} |')
        summaries.append({'sample_type': sample_type, 'samples': len(group), 'median_protein_groups': median,
                          'samples_above_service': wins, 'median_percent_gain': group.percent_gain.median() if benchmark else None})
    pd.DataFrame(summaries).to_csv(output / 'group_summary.tsv', sep='\t', index=False)
    for extension in ['png', 'pdf']:
        fig.savefig(output / f'protein_counts.{extension}', dpi=180)
    plt.close(fig)
    lines += ['', 'Service comparison is optional. Empty benchmark cells mean no service matrix was supplied; '
              'percentage gain is undefined when the service count is zero. Service settings may differ, '
              'so count differences alone do not establish matched empirical false-discovery rates.', '',
              'Protein and gene matrices are in matrices/<sample_type>/. Differential results and service benchmarks have separate folders. Detailed evidence is retained under analysis_details/groups/.']
    (output / 'REPORT.md').write_text('\n'.join(lines) + '\n')
    for path in args.groups:
        folder = Path(path)
        group = pd.read_csv(folder / 'sample_counts.tsv', sep='\t', dtype=str).sample_type.iloc[0]
        for name in ('matrices', 'differential_expression', 'benchmark', 'quality_control'):
            source = folder / 'delivery' / name
            if source.is_dir(): shutil.copytree(source, output / name / group)
        shutil.copy2(folder / 'delivery/README.md', output / ('README_' + group + '.md'))
    (output / 'quality_control').mkdir(exist_ok=True)
    labels_path = output / 'quality_control/plot_sample_labels.csv'
    pd.DataFrame(label_rows).to_csv(labels_path, index=False)
    from matrix_catalog import export_matrix_formats
    export_matrix_formats(labels_path)
    (output / 'README.md').write_text('Plot labels omit only the shared sample-ID prefix; quality_control/plot_sample_labels.csv maps every label to the full sample ID. Open START_HERE.html. Matrices are in matrices/<sample_type>/. See REPORT.md and README_<sample_type>.md for methods and interpretation.\n')
    write_dashboard(output, 'DIA-PASEF results')
    build_matrix_catalog(output)
    (output / 'FINAL_RESULTS_COMPLETE').write_text('DIA quantification and selected analyses complete.\n')
    print(f'Completed count table: {len(counts)} samples across {len(groups)} sample types')


if __name__ == '__main__':
    main()
