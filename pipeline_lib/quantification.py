"""Import quantification tables and compare detection on explicit sample matches."""
import hashlib
import csv
from pathlib import Path
import re
import numpy as np
import pandas as pd


def tokens(value):
    if value is None or pd.isna(value):
        return []
    return sorted({v.strip() for v in re.split(r'[;,]', str(value)) if v.strip()})


def canonical_group(value):
    return ';'.join(sorted({t.split('|')[1] if t.count('|') >= 2 else t for t in tokens(value)}))


def excluded_group(value):
    return any(re.match(r'(?i)^(?:contam|contaminant|rev|reverse|decoy|xxx)[_|-]', v)
               for token in tokens(value) for v in token.split('|'))


def read_table(path):
    """Check original headers before pandas can silently rename duplicates."""
    path = Path(path)
    with path.open('rb') as h:
        excel = h.read(4).startswith(b'PK')
    if excel:
        try:
            table = pd.read_excel(path, dtype=str, keep_default_na=False)
            header = pd.read_excel(path, header=None, nrows=1, keep_default_na=False).iloc[0].astype(str).tolist()
        except ImportError as error:
            raise ValueError('Excel input requires openpyxl; install it or export CSV/TSV.') from error
    else:
        with path.open(encoding='utf-8-sig', newline='') as h:
            first = h.readline(); delimiter = '\t' if '\t' in first else ','; h.seek(0)
            header = next(csv.reader(h, delimiter=delimiter), [])
        table = pd.read_csv(path, sep=delimiter, dtype=str, keep_default_na=False, encoding='utf-8-sig')
    if not header or len(set(header)) != len(header) or any(not v.strip() for v in header):
        raise ValueError('Quantification table has empty or duplicate column names')
    return table


def run_name(value):
    return re.sub(r'(?i)\.(d|raw|mzml)$', '', str(value).replace('\\', '/').rsplit('/', 1)[-1])


SUFFIXES = {'maxlfq': ' MaxLFQ Intensity', 'intensity': ' Intensity',
            'spectral_count': ' Spectral Count', 'value': ''}


def load_matrix(path, samples, quantity='auto', id_column=None):
    """Return values, annotations and explicit sample-column mapping.

    service_sample in the sample sheet overrides provider column matching.
    No fuzzy matching, duplicate collapsing, or implicit zero imputation.
    """
    table = read_table(path)
    if 'sample_id' not in samples or not samples.sample_id.is_unique or samples.sample_id.eq('').any():
        raise ValueError('Sample IDs must be nonempty and unique within this matrix')
    id_column = id_column or next((v for v in ('Protein.Group', 'Protein', 'FeatureID', 'Protein ID', 'protein_group', 'group_key') if v in table), None)
    if id_column not in table:
        raise ValueError('Quantification/service matrix needs Protein.Group, Protein, FeatureID or Protein ID')
    ids = table[id_column].copy()
    if id_column == 'Protein' and 'Indistinguishable Proteins' in table:
        ids = ids + ';' + table['Indistinguishable Proteins']
    excluded = ids.map(excluded_group)
    table = table.loc[~excluded].copy(); keys = ids.loc[~excluded].map(canonical_group)
    if keys.eq('').any() or keys.duplicated().any():
        raise ValueError('Empty or duplicate protein groups in quantification/service matrix')
    if quantity == 'auto':
        quantity = next((k for k in ('maxlfq', 'intensity', 'spectral_count') if any(c.endswith(SUFFIXES[k]) for c in table)), 'value')
    if quantity not in SUFFIXES:
        raise ValueError('Unknown quantity type: ' + quantity)
    suffix = SUFFIXES[quantity]; lookup = {}
    for column in table:
        if suffix and not column.endswith(suffix): continue
        if quantity == 'intensity' and column.endswith(' MaxLFQ Intensity'): continue
        if quantity == 'spectral_count' and column.endswith((' Unique Spectral Count', ' Total Spectral Count')): continue
        name = run_name(column[:-len(suffix)] if suffix else column)
        lookup.setdefault(name, []).append(column)
    selected, mapping = {}, []
    for sample in samples.to_dict('records'):
        explicit = sample.get('service_sample', '')
        aliases = [explicit] if explicit else [sample['sample_id'], sample.get('run', '')]
        matches = set(c for alias in aliases if alias for c in lookup.get(run_name(alias), []))
        if len(matches) != 1:
            raise ValueError(f"Quantification/service matrix: sample {sample['sample_id']!r} matched {len(matches)} columns. "
                             'Use service_sample in the sample sheet (column name without measurement suffix).')
        column = matches.pop()
        if column in selected.values(): raise ValueError('Multiple samples map to the same column: ' + column)
        selected[sample['sample_id']] = column
        mapping.append({'sample_id': sample['sample_id'], 'source_column': column, 'quantity': quantity})
    values = table[list(selected.values())].replace({'': np.nan, 'NA': np.nan, 'NaN': np.nan, 'nan': np.nan, 'N/A': np.nan, 'null': np.nan})
    try: values = values.apply(pd.to_numeric, errors='raise')
    except (ValueError, TypeError) as error: raise ValueError('Quantities must be numeric or blank/NA') from error
    array = values.to_numpy(dtype=float)
    if np.isinf(array).any() or (array < 0).any(): raise ValueError('Quantities must be finite nonnegative values or missing')
    values.columns = list(selected); values.index = pd.Index(keys, name='Protein.Group')
    annotations = pd.DataFrame(index=values.index)
    for dest, candidates in {'Genes': ('Genes', 'Gene'), 'Protein.Names': ('Protein.Names', 'Entry Name'), 'Description': ('Description', 'Protein.Description')}.items():
        source = next((v for v in candidates if v in table), None)
        annotations[dest] = table[source].to_numpy() if source else ''
    annotations['gene_mapping_status'] = annotations.Genes.map(lambda v: 'unmapped' if not tokens(v) else 'single_gene' if len(tokens(v)) == 1 else 'ambiguous')
    mapping = pd.DataFrame(mapping); mapping.attrs['excluded_decoy_contaminant_rows'] = int(excluded.sum())
    return values, annotations, mapping


def annotate_groups(report, groups, fasta=None):
    """Retain reported genes; supplement missing symbols from FASTA GN= fields."""
    genes, names, fallback = {}, {}, {}
    if fasta:
        with Path(fasta).open() as h:
            for line in h:
                if not line.startswith('>'): continue
                identifier = canonical_group(line[1:].split()[0]); match = re.search(r'(?:^|\s)GN=([^\s]+)', line)
                fallback.setdefault(identifier, set()).update(tokens(match.group(1)) if match else [])
    fields = [name for name in ('Protein.Group', 'Genes', 'Protein.Names') if name in report]
    for row in report[fields].drop_duplicates().to_dict('records'):
        key = canonical_group(row['Protein.Group'])
        genes.setdefault(key, set()).update(tokens(row.get('Genes', '')))
        names.setdefault(key, set()).update(tokens(row.get('Protein.Names', '')))
    records = []
    for key in groups:
        reported = genes.get(key, set()); mapped = set().union(*(fallback.get(a, set()) for a in key.split(';')))
        symbols = reported | mapped
        source = 'report+fasta' if reported and mapped else 'report' if reported else 'fasta' if mapped else 'unmapped'
        complete = (all(fallback.get(a) for a in key.split(';')) or (';' not in key and bool(reported))) if fasta else bool(reported)
        status = 'unmapped' if not symbols else 'ambiguous' if len(symbols) > 1 else 'single_gene' if complete else 'partial'
        records.append({'Protein.Group': key, 'Genes': ';'.join(sorted(symbols)), 'Protein.Names': ';'.join(sorted(names.get(key, set()))),
                        'gene_mapping_status': status, 'gene_annotation_source': source})
    return pd.DataFrame(records, columns=['Protein.Group','Genes','Protein.Names','gene_mapping_status','gene_annotation_source']).set_index('Protein.Group')


def gene_matrix(matrix, annotations):
    keep = annotations.gene_mapping_status.eq('single_gene')
    selected = matrix.loc[keep].copy(); selected.index = annotations.loc[keep, 'Genes']
    result = selected.groupby(level=0).sum(min_count=1); result.index.name = 'Gene'
    return result


def benchmark_detection(matrix, service):
    if set(matrix.columns) != set(service.columns): raise ValueError('Benchmark matrices require the same matched samples')
    if '__union__' in matrix.columns: raise ValueError('Sample name __union__ is reserved for benchmark summaries')
    rows, membership = [], []
    for sample in [*matrix.columns, '__union__']:
        new = set(matrix.index[matrix.gt(0).any(axis=1) if sample == '__union__' else matrix[sample].gt(0)])
        old = set(service.index[service.gt(0).any(axis=1) if sample == '__union__' else service[sample].gt(0)])
        shared, gained, lost = new & old, new - old, old - new
        rows.append({'sample_id': sample, 'pipeline_detected': len(new), 'service_detected': len(old), 'shared': len(shared),
                     'pipeline_only': len(gained), 'service_only': len(lost), 'net_difference': len(new)-len(old),
                     'percent_change': 100*(len(new)-len(old))/len(old) if old else np.nan})
        for group in sorted(new | old):
            membership.append({'sample_id': sample, 'Protein.Group': group, 'status': 'shared' if group in shared else 'pipeline_only' if group in gained else 'service_only'})
    return pd.DataFrame(rows), pd.DataFrame(membership, columns=['sample_id','Protein.Group','status'])


def write_benchmark(matrix, service, output, mapping=None):
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    summary, membership = benchmark_detection(matrix, service)
    summary.to_csv(output/'detection_summary.csv', index=False); membership.to_csv(output/'protein_membership.csv', index=False)
    service.to_csv(output/'service_quantification.csv')
    if mapping is not None: mapping.to_csv(output/'sample_mapping.csv', index=False)
    (output/'README.md').write_text('# Service detection benchmark\n\nPositive values count as detected. Rows are exact canonical protein groups; groups are not expanded into accessions. pipeline_only and service_only list gains and losses; net_difference is pipeline minus service. The union counts each group once across matched samples. Percentage change is blank when the service count is zero. Decoys/contaminants are excluded on import.\n\nSearch databases, protein grouping, identification FDR, transfer-between-runs and quantification methods may differ. More detections alone do not establish better accuracy. Counts are not differential abundance. Review the recorded sample mapping and input provenance.\n')
    return summary


def file_provenance(path):
    path = Path(path); digest = hashlib.sha256()
    with path.open('rb') as h:
        for chunk in iter(lambda: h.read(1024*1024), b''): digest.update(chunk)
    return {'path': str(path.resolve()), 'sha256': digest.hexdigest(), 'bytes': path.stat().st_size}
