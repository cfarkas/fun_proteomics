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
| Sample assignments and gene mapping audit | `quality_control/` for imported matrices; corresponding cohort reports for raw workflows |

Both matrix formats preserve identical identifiers, values and missing cells. DIA and service intensity matrices contain abundances, not spectral counts. DDA protein counts are a different measurement and use count-based statistics. [Full output definitions](../OUTPUTS.md).

## Comparisons

Each unordered condition pair is tested once. The table names identify baseline and comparison condition; positive logFC means higher abundance in the second condition. Three groups yield three protein contrasts and three gene contrasts when unambiguous gene mappings are available.

Intensity analysis uses limma-trend on log2 positive intensities, no imputation, and at least two observed biological replicates in each condition by default. Missing or untestable features retain their status and blank p-values. BH FDR is calculated separately within each contrast and feature level. These models do not automatically adjust for repeated subjects or batches.

Gene matrices sum positive protein-group quantities assigned unambiguously to one gene. Ambiguous or incomplete mappings are reported and excluded from the rollup. This is not RNA expression or a separately controlled gene-identification FDR.

## Service benchmark

A protein group is detected when its quantity is positive. The benchmark matches exact canonical groups and the same samples. It reports shared groups, pipeline-only groups, service-only groups, and pipeline minus service totals. The union row counts each group once across samples. More detections alone do not demonstrate better accuracy; database, grouping, transfer and filtering settings can differ.
