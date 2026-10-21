"""Publish lossless CSV/TSV companions and a single index of delivered matrices."""
import argparse
import csv
import hashlib
import html
from pathlib import Path
import os
import shutil
import tempfile
from urllib.parse import quote


def export_matrix_formats(path):
    """Change only the delimiter; preserve IDs, blanks and numeric text exactly."""
    path = Path(path)
    if path.suffix not in {'.csv', '.tsv'}:
        raise ValueError('Matrix format must be CSV or TSV: ' + str(path))
    destination = path.with_suffix('.tsv' if path.suffix == '.csv' else '.csv')
    temporary = None
    try:
        with path.open(encoding='utf-8-sig', newline='') as source, tempfile.NamedTemporaryFile(
                mode='w', encoding='utf-8', newline='', dir=path.parent, delete=False) as target:
            temporary = Path(target.name)
            reader = csv.reader(source, delimiter=',' if path.suffix == '.csv' else '\t', strict=True)
            writer = csv.writer(target, delimiter='\t' if destination.suffix == '.tsv' else ',', lineterminator='\n')
            width = None
            for row in reader:
                if width is None:
                    width = len(row)
                if not width or len(row) != width:
                    raise ValueError('Ragged or empty matrix: ' + str(path))
                writer.writerow(row)
            if width is None:
                raise ValueError('Empty matrix: ' + str(path))
        # Existing companions are refreshed, never used as a source of values.
        os.chmod(temporary, path.stat().st_mode & 0o777)
        os.replace(temporary, destination)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()
    return path.with_suffix('.csv'), path.with_suffix('.tsv')


def matrix_files(base):
    """One authoritative file per stem; CSV wins when both formats exist."""
    paths = sorted(p for p in base.rglob('*') if p.is_file() and p.suffix in {'.csv', '.tsv'})
    return [p for p in paths if p.suffix == '.csv' or not p.with_suffix('.csv').exists()]


def build_matrix_catalog(output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    rows, seen = [], set()

    def add(path, cohort, mode, feature, measure, source=None):
        csv_path, tsv_path = export_matrix_formats(path)
        rel = csv_path.relative_to(output).as_posix()
        if rel in seen:
            return
        seen.add(rel)
        rows.append(dict(cohort=cohort, mode=mode, feature=feature, measure=measure,
                         path=rel, source=(source or path).relative_to(output).as_posix(),
                         csv_path=rel, tsv_path=tsv_path.relative_to(output).as_posix(),
                         csv_sha256=hashlib.sha256(csv_path.read_bytes()).hexdigest(),
                         tsv_sha256=hashlib.sha256(tsv_path.read_bytes()).hexdigest()))

    bases = [output/'main_results/counts'] + sorted((output/'optional_analyses').glob('*/counts'))
    for base in bases:
        for path in matrix_files(base):
            parts = path.relative_to(base).parts
            if len(parts) != 4:
                continue
            cohort, mode, transformation, name = parts
            export_matrix_formats(path)
            destination = path
            if transformation == 'raw' and name in {'protein_counts.csv', 'gene_counts.csv', 'peptide_counts.csv'}:
                destination = output/'matrices'/cohort/mode/name
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(path, destination)
                if hashlib.sha256(path.read_bytes()).digest() != hashlib.sha256(destination.read_bytes()).digest():
                    raise ValueError('Matrix copy checksum mismatch: ' + str(path))
            add(destination, cohort, mode, name.split('_')[0],
                'spectral_count' if transformation == 'raw' else transformation, path)
    for path in matrix_files(output/'matrices'):
        parts = path.relative_to(output/'matrices').parts
        add(path, parts[0] if len(parts) > 1 else '',
            parts[1] if len(parts) > 2 and parts[1] in ('dda', 'ddaplus') else 'DIA/import',
            path.name.split('_')[0], 'intensity' if 'intensit' in path.name else 'spectral_count')
    for path in sorted((output/'differential_expression').rglob('log2_intensities.csv')):
        parts = path.relative_to(output/'differential_expression').parts
        add(path, parts[0] if len(parts) > 2 else '', 'DIA/import', parts[-2], 'log2_intensity')
    for path in sorted((output/'benchmark').rglob('service_quantification.csv')):
        parts = path.relative_to(output/'benchmark').parts
        add(path, parts[0] if len(parts) > 1 else '', 'Service', 'protein', 'service_quantity')
    fields = ['cohort', 'mode', 'feature', 'measure', 'path', 'source', 'csv_path', 'tsv_path', 'csv_sha256', 'tsv_sha256']
    with (output/'MATRICES.csv').open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    export_matrix_formats(output/'MATRICES.csv')
    body = ''.join('<tr>' + ''.join('<td>' + html.escape(r[k]) + '</td>' for k in fields[:4]) +
                   ''.join(f'<td><a href="{quote(r[k])}">{label}</a></td>' for k, label in [('csv_path', 'CSV'), ('tsv_path', 'TSV')]) +
                   '<td><small>' + html.escape(r['path']) + '</small></td></tr>' for r in rows)
    (output/'MATRICES.html').write_text('<!doctype html><meta charset="utf-8"><title>Count and abundance matrices</title>'
        '<style>body{font:16px system-ui;margin:3rem;color:#173d45}td,th{padding:.7rem;border-bottom:1px solid #ccd;text-align:left}a{color:#075c8c}</style>'
        '<h1>Count and abundance matrices</h1><p>Every matrix has matching <b>CSV and TSV</b> versions. Identifiers, values and missing cells are identical; only the delimiter changes. '
        'Primary matrices are in <b>matrices/</b>. DDA counts are discrete spectral counts; DIA intensities are continuous quantities. Missing intensities stay blank.</p>'
        '<p>Index: <a href="MATRICES.csv">CSV</a> · <a href="MATRICES.tsv">TSV</a> · <a href="START_HERE.html">All results</a></p>'
        '<table><thead><tr><th>Cohort</th><th>Search</th><th>Feature</th><th>Measure</th><th>CSV</th><th>TSV</th><th>Location</th></tr></thead><tbody>' + body + '</tbody></table>')
    dashboard = output/'START_HERE.html'
    if dashboard.is_file():
        text = dashboard.read_text()
        if 'href="MATRICES.html"' not in text:
            banner = '<p style="padding:1rem;background:#e3f2ef"><a href="MATRICES.html"><b>Find matrices — CSV and TSV</b></a> · <a href="MATRICES.csv">Index CSV</a> · <a href="MATRICES.tsv">Index TSV</a></p>'
            text = text.replace('<body>', '<body>' + banner, 1) if '<body>' in text else text + banner
            dashboard.write_text(text)
    return rows


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results')
    args = parser.parse_args()
    print('Indexed', len(build_matrix_catalog(args.results)), 'matrices in CSV and TSV')
