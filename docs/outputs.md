# Find and understand your results

Start with `PROJECT/START_HERE.html`. Use **`PROJECT/MATRICES.html`** for direct CSV and TSV links. These two shortcuts point to the current delivery, including after successful recovery into a new result generation.

```text
my-study/
  config/
    run.json                 saved analysis choices
    samples.csv              normalized sample assignments for matrix imports
  START_HERE.html             open the current report
  MATRICES.html               open the matrix catalogue
  results/
    final_results/           reports, matrices, comparisons, QC, benchmarks
    .fun_proteomics/          execution state, logs and integrity records
  logs/                      browser-launched analysis log
  .fun_proteomics/            project lock and launch provenance
```

The `results/` directory is the run root. Its internal layout depends on the analysis type. Raw workflows also retain Nextflow work/cache and supporting analysis details. Avoid deleting these if you need continuation.

| You need | Where to look inside `results/final_results/` |
| --- | --- |
| All matrices | `MATRICES.html`, `MATRICES.csv`, `MATRICES.tsv` |
| Imported protein quantities | `matrices/protein_intensities.csv` and `.tsv` |
| Imported gene quantities | `matrices/gene_intensities.csv` and `.tsv` |
| Raw DIA quantities | `matrices/<sample_type>/protein_intensities.*` and `gene_intensities.*` |
| Raw DDA count matrix | `matrices/<cohort>/<dda\|ddaplus>/protein_counts.*` |
| Every condition comparison | `results_index.csv` / `.tsv`, linked from the report |
| Imported matrix differential tables | `differential_expression/<protein\|gene>/<contrast>/` |
| Imported matrix detection benchmark | `benchmark/detection_summary.csv` and `protein_membership.csv` |
| Excel sheet, source column mapping, removed empty duplicates | `quality_control/import_audit.json` and `input_sample_mapping.csv/.tsv` |
| Sample assignments and gene mapping audit | `quality_control/` for imported matrices; corresponding cohort reports for raw workflows |

Both matrix formats preserve identical identifiers, values and missing cells. DIA and service intensity matrices contain abundances, not spectral counts. DDA protein counts are a different measurement and use count-based statistics. [Full output definitions](../OUTPUTS.md).

## Table previews

The following snippets illustrate formats and selected columns; they are not study results. The generated files include additional annotations and audit fields. CSV and TSV contain the same table with different separators.

**`matrices/protein_intensities.csv`** — protein-group abundance; selected columns:

```csv
Protein.Group,Genes,C1,C2,T1,T2
P04637,TP53,120000,135000,260000,245000
P60709,ACTB,480000,505000,490000,510000
P02768,ALB,90000,85000,,100000
```

**`matrices/gene_intensities.csv`** — gene abundance; selected columns:

```csv
Gene,C1,C2,T1,T2
TP53,120000,135000,260000,245000
ACTB,480000,505000,490000,510000
```

**`protein_counts.csv / protein_spectral_counts.csv`** — integer spectral counts; selected columns:

```csv
FeatureID,C1,C2,T1,T2
P04637,12,15,28,30
P60709,45,50,47,49
```

**`all_results.csv`** — one contrast; selected columns:

```csv
FeatureID,Genes,logFC,PValue,FDR,test_status,significant
P04637,TP53,1.0,0.01,0.04,tested,True
P60709,ACTB,0.02,0.8,0.9,tested,False
```

`significant_results.csv` has the same columns and retains only rows passing the configured FDR threshold.

**`results_index.csv`** — links to comparisons; selected columns:

```csv
feature,baseline,comparison_condition,tested_features,significant_features,full_results
protein,Control,Treatment,100,4,differential_expression/protein/Control_vs_Treatment/all_results.csv
gene,Control,Treatment,90,3,differential_expression/gene/Control_vs_Treatment/all_results.csv
```

**`MATRICES.csv`** — matrix download paths; selected columns:

```csv
feature,measure,csv_path,tsv_path
protein,intensity,matrices/protein_intensities.csv,matrices/protein_intensities.tsv
gene,intensity,matrices/gene_intensities.csv,matrices/gene_intensities.tsv
```

## Comparisons

Each unordered condition pair is tested once. The table names identify baseline and comparison condition; positive logFC means higher abundance in the second condition. Three groups yield three protein contrasts and three gene contrasts when unambiguous gene mappings are available.

Intensity analysis uses limma-trend on log2 positive intensities, no imputation, and at least two observed biological replicates in each condition by default. Missing or untestable features retain their status and blank p-values. BH FDR is calculated separately within each contrast and feature level. With an explicit `subject_column`, imported matrix comparisons use within-donor log2 differences. `detected_pairs` reports the number of complete observed pairs; unmatched values cannot satisfy the replicate threshold. Per-contrast `paired_samples.csv/.tsv` and `paired_log2_differences.csv/.tsv` show the matching and model input. No other batch adjustment is inferred.

Gene matrices sum positive protein-group quantities assigned unambiguously to one gene. Ambiguous or incomplete mappings are reported and excluded from the rollup. This is not RNA expression or a separately controlled gene-identification FDR.

## Service benchmark

A protein group is detected when its quantity is positive. The benchmark matches exact canonical groups and the same samples. Incomplete groups reconstructed from clipped FASTA headers are retained in reported counts and excluded from exact-group identity comparisons. It reports shared groups, pipeline-only groups, service-only groups, and pipeline minus service totals. The union row counts each group once across samples. More detections alone do not demonstrate better accuracy; database, grouping, transfer and filtering settings can differ.

For saved configuration and sample-file previews, see [input files](input_files.md). For the service worksheets, sample mappings and benchmark result example, see [benchmarking](benchmarking.md).
