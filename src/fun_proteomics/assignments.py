"""Discover sample cards and retain an explicit, auditable condition assignment."""
from pathlib import Path

from .project import digest

SOURCE_FIELDS = ('mode', 'matrix', 'samples', 'raw_root', 'quantity')
METADATA = {'protein.group', 'protein', 'featureid', 'protein id', 'protein_group',
            'group_key', 'genes', 'gene', 'protein.names', 'entry name', 'description',
            'protein.description', 'indistinguishable proteins', 'protein length',
            'length', 'molecular weight', 'protein probability', 'top peptide probability',
            'combined total peptides', 'combined unique peptides', 'combined total spectral count',
            'combined unique spectral count', 'combined spectral count', 'number of proteins',
            'q-value', 'q.value', 'pg.q.value', 'global.pg.q.value'}


def source(cfg):
    value = {key: cfg.get(key, '') for key in SOURCE_FIELDS}
    for key in ('matrix', 'samples'):
        if cfg.get(key):
            value[key + '_sha256'] = digest(cfg[key])
    return value


def discover(cfg, workflow):
    from . import study
    cfg = study.normalize(cfg)
    if cfg['mode'] == 'continue' or cfg.get('study_config'):
        raise ValueError('An existing run or native configuration supplies its own assignments')
    study.imports(workflow)
    import pandas as pd
    from quantification import read_table, load_matrix, SUFFIXES, run_name
    before = source(cfg)
    if cfg.get('samples'):
        frame = read_table(cfg['samples'])
        if 'sample_id' not in frame or frame.empty or frame.sample_id.eq('').any() or not frame.sample_id.is_unique:
            raise ValueError('The sample table needs nonempty, unique sample_id values')
        if 'condition' not in frame:
            frame['condition'] = ''
    elif cfg['mode'] in ('dia', 'dda'):
        frame = pd.DataFrame(study.scan(cfg.get('raw_root', '')))
    else:
        if not cfg.get('matrix'):
            raise ValueError('Select a protein matrix before discovering samples')
        table = read_table(cfg['matrix'])
        quantity = cfg['quantity']
        if quantity == 'auto':
            quantity = next((k for k in ('maxlfq', 'intensity', 'spectral_count')
                             if any(c.endswith(SUFFIXES[k]) for c in table)), 'value')
        suffix = SUFFIXES[quantity]
        columns = []
        for column in table:
            if suffix:
                if not column.endswith(suffix): continue
                if quantity == 'intensity' and column.endswith(' MaxLFQ Intensity'): continue
                if quantity == 'spectral_count' and column.endswith((' Unique Spectral Count', ' Total Spectral Count')): continue
            else:
                if column.lower() in METADATA: continue
                values = table[column].replace(['', 'NA', 'NaN', 'nan', 'N/A', 'null'], float('nan')).dropna()
                if values.empty or pd.to_numeric(values, errors='coerce').isna().any(): continue
            columns.append(column[:-len(suffix)] if suffix else column)
        if not columns:
            raise ValueError('No sample quantity columns detected. Supply a sample table with explicit sample_id / service_sample mappings.')
        frame = pd.DataFrame([dict(sample_id=run_name(c), run=c, condition='') for c in columns])
    if cfg['mode'] == 'matrix':
        _, _, mapping = load_matrix(cfg['matrix'], frame, cfg['quantity'])
        # Keep the source alias stable when the display/analysis name is edited.
        for index, item in enumerate(mapping.to_dict('records')):
            suffix = SUFFIXES[item['quantity']]
            alias = item['source_column'][:-len(suffix)] if suffix else item['source_column']
            frame.loc[frame.index[index], 'run'] = alias
            frame.loc[frame.index[index], 'service_sample'] = alias
    else:
        if 'raw_path' not in frame:
            raise ValueError('Raw sample tables need raw_path')
        base = Path(cfg.get('raw_root') or Path(cfg.get('samples', cfg['project'])).parent)
        frame['raw_path'] = [str((base / Path(p).expanduser()).resolve()) for p in frame.raw_path]
    if len(frame) > 2000:
        raise ValueError('Select at most 2000 samples')
    if 'detected_sample' not in frame: frame['detected_sample'] = frame.sample_id
    rows = frame.fillna('').astype(str).to_dict('records')
    if before != source(cfg):
        raise ValueError('Inputs changed during discovery. Discover samples again.')
    conditions = list(dict.fromkeys(row['condition'].strip() for row in rows if row['condition'].strip()))
    for row in rows: row['condition'] = row['condition'].strip()
    return dict(source=before, conditions=conditions or ['Control', 'Treatment'], rows=rows)


def rows_for(cfg):
    board = cfg['assignment_board']
    if not isinstance(board, dict) or set(board) != {'source', 'conditions', 'rows'}:
        raise ValueError('Invalid sample assignment board; discover samples again')
    if board['source'] != source(cfg):
        raise ValueError('Input files or quantity settings changed after sample discovery. Discover samples again before saving.')
    conditions = board['conditions']
    if (not isinstance(conditions, list) or any(not isinstance(c, str) or not c.strip() or c != c.strip()
            or any(char in c for char in '\r\n\t') for c in conditions)
            or len(set(conditions)) != len(conditions)):
        raise ValueError('Condition names must be nonempty and unique')
    rows = board['rows']
    if not isinstance(rows, list) or not rows or len(rows) > 2000:
        raise ValueError('Discover samples before assigning conditions')
    ids = set()
    for row in rows:
        if not isinstance(row, dict) or any(not isinstance(k, str) or not isinstance(v, str) for k, v in row.items()):
            raise ValueError('Sample assignment fields must be text')
        name = row.get('sample_id', '')
        if not name.strip() or name != name.strip() or any(c in name for c in '\r\n\t') or name in ids:
            raise ValueError('Sample IDs must be nonempty and unique, including unassigned samples')
        ids.add(name)
        if row.get('condition', '') and row['condition'] not in conditions:
            raise ValueError('A sample refers to an unknown condition')
    selected = [dict(row) for row in rows if row.get('condition')]
    if not selected:
        raise ValueError('Assign at least one sample to a condition. Unassigned samples are excluded.')
    return selected


def service_samples_for(cfg, samples, workflow):
    """Align a separate provider mapping to the selected (possibly renamed) cards."""
    from . import study
    study.imports(workflow)
    from quantification import read_table
    import pandas as pd
    if not cfg.get('service_samples'):
        return samples
    provider = read_table(cfg['service_samples'])
    if 'sample_id' not in provider or not provider.sample_id.is_unique:
        raise ValueError('Service sample mapping needs unique sample_id values')
    if 'assignment_board' not in cfg:
        return provider
    indexed = provider.set_index('sample_id', drop=False)
    rows = []
    for sample in samples.to_dict('records'):
        original = sample.get('detected_sample', sample['sample_id'])
        if original not in indexed.index:
            raise ValueError('No service mapping for detected sample ' + original)
        row = indexed.loc[original].to_dict()
        # Preserve the provider alias even when no explicit service_sample was needed.
        row['service_sample'] = row.get('service_sample') or row.get('run') or original
        row['sample_id'] = sample['sample_id']
        rows.append(row)
    return pd.DataFrame(rows)
