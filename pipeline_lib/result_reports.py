"""Publish explicit test families and condition summaries without refitting models."""
import json, re, shutil, tempfile
from datetime import datetime
from pathlib import Path
import pandas as pd
from html import escape
from urllib.parse import quote

CATEGORIES = dict(protein="proteins", gene="genes", peptide="peptides", modification_gene="modifications_by_gene", modification_protein="modifications_by_protein")

def organize(out, preserve_original=False):
    out = Path(out).resolve()
    source = out / "differential"
    if not source.exists():
        write_guide(out, [], json.loads((out / "analysis_settings.json").read_text()))
        return
    if (source / "results_index.csv").exists():
        return
    job = json.loads((out / "analysis_settings.json").read_text())
    stage = Path(tempfile.mkdtemp(prefix=".differential_reports_", dir=out))
    index, mapping = [], []
    try:
        for group in job["groups"]:
            cohort, mode = group["cohort"], group["mode"]
            slug = "ddaplus" if mode == "DDA+" else "dda"
            raw = (out / group["directory"]).resolve()
            samples = pd.read_csv(raw / "samples.csv", keep_default_na=False)
            conditions = list(dict.fromkeys(samples.condition))
            labels, used = {}, set()
            for condition in conditions:
                label = job["differential"].get("condition_labels", {}).get(condition) or condition
                base = re.sub(r"[^A-Za-z0-9_-]", "_", label)
                label, number = base, 1
                while label in used:
                    label = f"{base}_{number}"
                    number += 1
                labels[condition] = label
                used.add(label)
                mapping.append(dict(Cohort=cohort, Mode=mode, column_label=label, condition=condition, samples=int((samples.condition==condition).sum())))
            for kind, category in job.get("feature_categories", CATEGORIES).items():
                counts = pd.read_csv(raw / f"{kind}_counts.csv", keep_default_na=False).set_index("FeatureID")
                normalized_path = raw / "models" / category / "counts_model_normalized_cpm.csv"
                if not normalized_path.exists():
                    continue
                normalized = pd.read_csv(normalized_path, keep_default_na=False).set_index("FeatureID")
                means = pd.DataFrame(index=counts.index)
                for condition, label in labels.items():
                    ids = samples.loc[samples.condition==condition, "sample_id"].tolist()
                    values = counts[ids].apply(pd.to_numeric)
                    means[f"mean_count_{label}"] = values.mean(axis=1)
                    means[f"mean_normalized_CPM_{label}"] = normalized[ids].apply(pd.to_numeric, errors="coerce").mean(axis=1)
                    means[f"detected_samples_{label}"] = (values>0).sum(axis=1)
                    means[f"total_samples_{label}"] = len(ids)
                for method in ("edgeR", "edgeRun", "limma"):
                    original = source / method / cohort / slug
                    for path in sorted(original.glob(f"*_{category}.csv")):
                        table = pd.read_csv(path, keep_default_na=False)
                        if "PValue" not in table:
                            continue
                        omnibus = "omnibus" in path.stem
                        family = "omnibus" if omnibus else "pairwise"
                        if omnibus:
                            annotation = [c for c in counts.columns if c not in samples.sample_id.tolist()]
                            allowed = set(annotation + ["FeatureID", "logCPM", "AveExpr", "F", "t", "B", "PValue", "adjusted_p_value", "FDR", "significant", "method", "contrast", "test_status"])
                            table = table[[c for c in table if c in allowed and not c.lower().startswith("logfc")]]
                        if "adjusted_p_value" in table:
                            table["FDR"] = table["adjusted_p_value"]
                        elif "adj.P.Val" in table:
                            table["FDR"] = table["adj.P.Val"]
                        if "FDR" not in table:
                            raise ValueError(f"Missing adjusted p-values: {path}")
                        table = table.drop(columns=[c for c in ("adjusted_p_value", "adj.P.Val", "qvalue", "q_value", "P.Value") if c in table])
                        table["method"] = method
                        table["test_type"] = family
                        sparse = kind.startswith("modification_") or kind in job.get("sparse_kinds", [])
                        table["model"] = ("edgeR v4 small-count quasi-likelihood GLM" if sparse else "edgeR quasi-likelihood GLM") if method=="edgeR" else "unconditional exact test" if method=="edgeRun" else "limma-voom"
                        table["Cohort"], table["Mode"], table["FeatureType"] = cohort, mode, kind
                        table["count_unit"] = job.get("count_units", {}).get(kind, "distinct modified peptide sequences" if kind.startswith("modification_") else "spectral counts")
                        table["normalization"] = "TMMwsp" if sparse else "TMM"
                        table = table.merge(means, left_on="FeatureID", right_index=True, how="left", validate="one_to_one")
                        if table.filter(like="mean_count_").isna().any().any():
                            raise ValueError(f"Missing condition counts: {path}")
                        filename = path.stem.removesuffix("_"+category)+".csv"
                        relative = Path(method) / family / cohort / slug / category / filename
                        dest = stage / relative
                        dest.parent.mkdir(parents=True, exist_ok=True)
                        table.to_csv(dest, index=False)
                        significant = table[table.significant.astype(str).str.lower()=="true"]
                        sigrelative = relative.parent / "significant" / filename
                        (stage/sigrelative).parent.mkdir(parents=True, exist_ok=True)
                        significant.to_csv(stage/sigrelative, index=False)
                        index.append(dict(method=method, test_type=family, cohort=cohort, mode=mode, feature_type=kind, comparison="any_condition_difference" if omnibus else table.contrast.iloc[0] if len(table) else filename, tested_features=int(pd.to_numeric(table.PValue,errors="coerce").notna().sum()), reported_features=len(table), significant_features=len(significant), all_results=str(relative), significant_results=str(sigrelative)))
                    for note in original.glob(f"*_{category}.txt"):
                        dest=stage/method/"notes"/cohort/slug/category/note.name
                        dest.parent.mkdir(parents=True,exist_ok=True)
                        shutil.copy2(note,dest)
        if not index:
            shutil.rmtree(stage)
            write_guide(out, [], job)
            return
        pd.DataFrame(index).to_csv(stage/"results_index.csv", index=False)
        pd.DataFrame(mapping).to_csv(stage/"condition_columns.csv", index=False)
        (stage/"README.md").write_text("""# Differential results

Choose a method, then `pairwise/` (named group comparisons) or `omnibus/` (any condition difference), then cohort, search mode and feature type. Each CSV contains tested features or explicitly labeled untestable rows; its `significant/` counterpart contains only FDR-significant rows. `results_index.csv` lists every test and its files. There is no combined discovery table mixing test types.

Low coverage does not fail the study. Full tables retain nonsignificant results. Rows that cannot be tested have an explicit `test_status`, blank p-values/FDR and `significant=false`; they are not evidence of no effect.

`PValue` is the unadjusted p-value. **FDR is the only adjusted-p-value column**, using Benjamini–Hochberg separately for each table. `significant` uses the configured FDR cutoff. Pairwise `logFC` is comparison minus reference. Omnibus tests have no single direction: use named pairwise tests to establish which comparisons are significant.

`mean_count_GROUP` is the arithmetic mean of raw counts across all samples of that condition, including zeros. `mean_normalized_CPM_GROUP` is the mean model-normalized CPM. `detected_samples_GROUP` reports samples with counts greater than zero, out of `total_samples_GROUP`. `condition_columns.csv` maps labels to the full condition names.

Modification counts are distinct modified peptide sequences per gene/protein and modification type. Detection does not establish biological absence in other groups, site occupancy or a significant pairwise difference. Condition averages are descriptive, not additional statistical tests.
""")
        if preserve_original:
            backup = out.parent/"analysis_details"/("previous_differential_layout_"+datetime.now().strftime("%Y%m%d_%H%M%S_%f"))
        else:
            backup = out/(stage.name+"_original")
        backup.parent.mkdir(parents=True, exist_ok=True)
        source.rename(backup)
        try:
            stage.rename(source)
        except Exception:
            backup.rename(source)
            raise
        if not preserve_original:
            shutil.rmtree(backup)
        write_guide(out, index, job)
        return len(index)
    except Exception:
        if stage.exists():
            shutil.rmtree(stage)
        raise


def write_guide(out, index, job):
    def link(path, text):
        return '<a href="'+quote(str(path), safe="/")+'">'+escape(str(text))+'</a>'
    names = dict(protein="Proteins",gene="Genes (proteomics)",peptide="Peptides",modification_gene="Modified peptides per gene",modification_protein="Modified peptides per protein")
    names.update(isoform="Isoform-discriminating evidence")
    rows=[]
    for entry in index:
        cells=[entry["method"], "Overall condition test" if entry["test_type"]=="omnibus" else "Group comparison", entry["cohort"], entry["mode"], names[entry["feature_type"]], "Any condition difference" if entry["test_type"]=="omnibus" else entry["comparison"], entry["tested_features"], entry["significant_features"]]
        rows.append('<tr>'+''.join('<td>'+escape(str(x))+'</td>' for x in cells)+'<td>'+link('differential/'+entry['all_results'],'Full result table')+' · '+link('differential/'+entry['significant_results'],'Significant features')+'</td></tr>')
    counts=''.join('<li>'+link(kind+'_counts.csv',label)+'</li>' for kind,label in [('protein','Protein counts'),('gene','Gene counts'),('peptide','Peptide counts'),('modification_gene','Distinct modified peptides per gene'),('modification_protein','Distinct modified peptides per protein'),('modification','Supporting modification PSM counts')])
    if job.get('feature_categories'):
        counts=''.join('<li>'+link(kind+'_counts.csv',names[kind])+'</li>' for kind in job['feature_categories'])
    page='<!doctype html><html lang="en"><meta charset="utf-8"><title>Start here — Proteomics results</title><style>body{font:16px/1.5 system-ui;margin:30px;color:#20392f}h1,h2{color:#173e32}a{color:#17684b}section{background:#f1f6f2;padding:18px;border-radius:10px;margin:18px 0}table{border-collapse:collapse;width:100%;font-size:14px}td,th{padding:9px;text-align:left;border-bottom:1px solid #d3ded6}th{background:#e6eee8}input{padding:12px;width:min(95%,600px)}.scroll{overflow:auto}</style><h1>Your proteomics results</h1><p>Start with significant features for your chosen method and group comparison. Use all tested features to inspect a particular gene or protein, even when it is not significant.</p><section><h2>1. Find differences between groups</h2><p><b>Group comparison:</b> tests a named pair of conditions. For <code>CTRL_vs_TENI</code>, positive logFC means higher in TENI; negative means higher in CTRL.</p><p><b>Overall condition test:</b> asks whether any condition differs. It does not identify which pair differs. Check the named group comparisons for that.</p><p><b>FDR</b> is the single Benjamini–Hochberg adjusted p-value column. The significance cutoff for this run is FDR &lt; @FDR@. Low coverage is not a run failure. Full tables retain nonsignificant results. Rows that cannot be tested have an insufficient-evidence status, blank p-values/FDR and significant=false; these are not evidence of no effect.</p><p>Keep DDA and DDA+ separate: they are different searches, not biological replicates.</p><label>Find a method, comparison or feature type<br><input id="search" type="search" placeholder="For example: edgeR CTRL_vs_NILO genes" oninput="const terms=this.value.toLowerCase().split(/\\s+/);document.querySelectorAll(\'tbody tr\').forEach(r=>r.hidden=!terms.every(t=>r.textContent.toLowerCase().includes(t)))"></label></section><div class="scroll"><table><thead><tr><th>Method</th><th>Test</th><th>Cohort</th><th>Search</th><th>Features</th><th>Comparison</th><th>Tested</th><th>Significant</th><th>Open CSV</th></tr></thead><tbody>@ROWS@</tbody></table></div><section><h2>2. See which conditions contain the evidence</h2><p>Every test table includes these labeled columns for each condition:</p><ul><li><b>mean_count_GROUP:</b> average raw count, including zeros.</li><li><b>mean_normalized_CPM_GROUP:</b> average count per million after model normalization.</li><li><b>detected_samples_GROUP / total_samples_GROUP:</b> how many samples have at least one count.</li></ul><p>For example, mean_count_CTRL = 0 and detected_samples_TENI = 3 out of total_samples_TENI = 4 means the feature was detected in three TENI samples and not detected in CTRL. This alone does not establish a significant difference or biological absence.</p><p>Modification counts measure distinct modified peptide sequences for each gene/protein and modification type. They do not measure modification occupancy. Gene results are derived from proteomics, not RNA expression.</p><p><a href="differential/condition_columns.csv">Full condition names and column labels</a></p></section><section><h2>3. Counts, normalized matrices and z-scores</h2><ul>@COUNTS@</ul><p><a href="matrices/">Normalized CPM and log2 CPM matrices</a> · <a href="../zscore/">Z-score matrices</a></p><p>Z-scores show relative patterns across samples within each feature. They are not p-values or differential tests.</p></section><p><a href="README.md">Run README and limitations</a> · <a href="analysis_summary.csv">Testing and filtering summary</a> · <a href="differential/results_index.csv">Complete test index</a> · <a href="differential/README.md">Detailed column definitions</a></p></html>'
    page=page.replace('@FDR@',escape(str(job['differential'].get('fdr',0.05)))).replace('@ROWS@',''.join(rows)).replace('@COUNTS@',counts)
    if job.get('feature_categories'):
        page=page.replace('Your proteomics results','Your optional-analysis statistics').replace('href="../zscore/"','href="zscore/"')
        page=page.replace('<h1>Your optional-analysis statistics</h1>','<h1>Your optional-analysis statistics</h1><p>These exploratory tests compare peptide counts or isoform-discriminating evidence, as selected. Isoform-associated counts do not establish isoform switching. HLA counts describe peptide evidence, not source-protein expression. The README records the precise count definitions and exclusions.</p>')
    if not index:
        page=page.replace('<h2>1. Find differences between groups</h2>','<h2>1. Differential testing was not performed</h2><p>Review the run README for disabled or skipped analyses. Count matrices are available below.</p>')
        page=re.sub(r'<a href="(?:differential/[^"]+|analysis_summary.csv)">([^<]+)</a>',r'\1 (not available)',page)
    if job.get('analysis') == 'immunopeptidomics':
        links=[]
        for group in job['groups']:
            folder='immunopeptidomics/'+group['cohort']+'/dda/'
            for filename,label in [('peptides.csv','Peptides'),('sample_summary.csv','Sample summary'),
                                   ('length_distribution.csv','Peptide lengths'),('binding_input.txt','Binding input'),('README.md','HLA report guide')]:
                links.append('<li>'+link(folder+filename,group['cohort']+' · '+label)+'</li>')
        section='<section><h2>HLA-enriched immunopeptidomics</h2><p>Start with the peptide evidence and sample summaries. Differential tests apply to peptide spectral counts only. Other count matrices describe source assignments, not source-protein expression.</p><p>Identification: 1% PSM FDR; no 1% peptide/protein FDR claim. This is separate from differential-test FDR. Match-between-runs is off. HLA binding is not predicted here.</p><ul>'+''.join(links)+'</ul></section>'
        page=page.replace('<h1>Your proteomics results</h1>','<h1>Your HLA peptide results</h1>'+section)
    optional=out/'optional_analyses'
    if optional.exists():
        links=[]
        if (optional/'statistics/START_HERE.html').is_file():
            links.append('<li>'+link('optional_analyses/statistics/START_HERE.html','Statistical comparisons for optional analyses — methods, contrasts and full tables')+'</li>')
        for path in sorted(optional.rglob('*')):
            if path.is_file() and path.suffix in ('.csv','.tsv','.md'):
                links.append('<li>'+link(path.relative_to(out),path.relative_to(optional))+'</li>')
        if links:
            section='<section><h2>Optional analyses</h2><p>Sequence-supported peptide and isoform evidence, and drug-combination interactions are reported separately. Each analysis includes its settings and interpretation guide.</p><ul>'+''.join(links)+'</ul></section>'
            page=page.replace('</html>',section+'</html>')
    (out/'START_HERE.html').write_text(page)
