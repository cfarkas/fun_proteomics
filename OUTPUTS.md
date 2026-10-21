# Results directory

## Matrix entry point and DIA outputs

Every delivered matrix is available in both **CSV and TSV**, with identical identifiers, values and missing cells. `MATRICES.html` has separate links for both formats; the index is also supplied as `MATRICES.csv` and `MATRICES.tsv`.

Start with `final_results/MATRICES.html` or its machine-readable `MATRICES.csv`.
For DDA, `matrices/<cohort>/<dda|ddaplus>/protein_counts.csv` is a checksum-verified copy
of the existing annotated raw count matrix. Peptide and selected gene matrices are beside it.
The historical detailed folders remain available. `python pipeline_lib/matrix_catalog.py /data/run/final_results`
adds this index to completed results without statistical refitting.

For DIA, `matrices/<sample_type>/protein_intensities.csv` contains `Protein.Group`, `Genes`,
protein names, gene mapping status/source and sample intensities. `gene_intensities.csv` sums quantities
only for unambiguous single-gene groups, with missingness retained and no duplication of ambiguous groups.
These are protein-derived gene rollups, not RNA counts or independently tested gene identification FDR.
Unmapped/ambiguous/partial assignments are listed under `quality_control/<sample_type>/`.

DIA comparisons use `differential_expression/<sample_type>/<protein|gene>/<reference>_vs_<comparison>/`.
`all_results.csv` retains insufficient-observation rows; `significant_results.csv` uses the selected BH cutoff.
Positive logFC means comparison minus reference. `results_index.csv`, model settings, log2 matrices,
normalization offsets, design and R session document each analysis. Every contrast and feature level has
its own FDR family. No imputation or subject/batch adjustment is inferred.

`benchmark/<sample_type>/detection_summary.csv` reports pipeline/service totals, shared, pipeline-only,
service-only, net difference and percentage change for each matched sample and the union. `protein_membership.csv`
lists each category. Zero baseline produces a blank percentage. Exact canonical group membership is required;
protein groups are never expanded into independent proteins. Search/grouping/FDR differences limit interpretation.
Native DIA evidence is under `analysis_details/`; raw acquisition files are never copied into deliverables.


Open **`<outdir>/final_results/START_HERE.html`**. The dashboard works offline and separates the study overview, main results, optional analyses, counts/QC and interpretation guidance. Statistical tables can be filtered by analysis, search mode, method and result status. Each test links to its full results and significant subset.

```text
<outdir>/
├── final_results/
│   ├── START_HERE.html
│   ├── README.md
│   ├── results_index.csv
│   ├── FINAL_RESULTS_COMPLETE
│   ├── main_results/
│   │   ├── differential_expression/<cohort>/<dda|ddaplus>/<comparison>/<feature_type>/
│   │   │   ├── <method>_all_results.csv
│   │   │   └── <method>_significant_results.csv
│   │   ├── counts/<cohort>/<dda|ddaplus>/
│   │   │   ├── raw/*_counts.csv
│   │   │   ├── normalized/*_normalized_cpm.csv
│   │   │   ├── normalized/*_log2_cpm.csv
│   │   │   └── zscore/*_zscore.csv
│   │   └── quality_control/
│   │       ├── sample_annotations.csv
│   │       ├── analysis_summary.csv
│   │       ├── comparisons.csv
│   │       └── condition_columns.csv
│   └── optional_analyses/
│       ├── gene_counts/{counts,comparisons}/  # optional gene and gene-modification results
│       ├── statistics/
│       │   ├── comparisons/<cohort>/<mode>/<comparison>/<peptide|isoform>/
│       │   │   ├── <method>_all_results.csv
│       │   │   └── <method>_significant_results.csv
│       │   └── counts/<cohort>/<mode>/<raw|matrices|zscore>/
│       ├── drug_interaction/<cohort>/<mode>/<protein|gene|peptide>/
│       │   ├── <method>_all_results.csv
│       │   └── <method>_significant_results.csv
│       ├── peptide_isoforms/<cohort>/<mode>/
│       │   ├── peptide_evidence.csv
│       │   ├── peptide_condition_summary.csv
│       │   └── isoform_evidence.csv
│       └── immunopeptidomics/<cohort>/dda/
│           ├── peptides.csv
│           ├── sample_summary.csv
│           ├── length_distribution.csv
│           └── binding_input.txt
├── analysis_details/
│   ├── <cohort>/<dda|ddaplus>/run/         # original search files
│   ├── <cohort>/<dda|ddaplus>/tables_*/    # detailed search counts
│   ├── reporting/<cohort>/<mode>/         # main model inputs, filters, assignments
│   ├── reporting/result_audit/           # settings, optional model diagnostics, file checksums
│   └── final/                           # search summary and workbook
├── logs/
└── database/search.fasta
```

Only generated analyses have result folders. Existing search-level coverage tables are included under `main_results/quality_control/` when available. `method_notes.csv` is present when a method cannot be estimated. The overall condition test uses `overall_condition_test/` in place of a named comparison. Methods stay together under the biological comparison, while cohorts and DDA/DDA+ remain separate.

`results_index.csv` records section, analysis, cohort, mode, evidence type, method, test type, comparison, reported/tested/significant feature counts, status, FDR cutoff and both result paths. It describes separate statistical families; it is not a merged list of discoveries. Dashboard totals count tables, not unique significant features across methods.

## Main analyses

Protein differential abundance is primary; supporting peptide and selected protein-modification comparisons follow. Gene-level results appear only under Optional analyses.

**Counts.** Protein and peptide matrices contain spectral counts. Gene counts sum peptide spectral counts only for unambiguous gene assignments, using Gene plus Mapped Genes. Gene results describe proteomics-derived evidence, not RNA expression. Modification-by-gene and modification-by-protein matrices count distinct modified peptide sequences per entity and modification type: repeated spectra, charges and sites do not multiply a sequence within one type. A peptide with two modification types can contribute once to each type. Ambiguous gene/protein assignments are excluded and audited. Supporting `modification_counts.csv` contains peptidoform PSM counts instead.

**Comparisons.** Positive pairwise `logFC` means comparison minus reference. Overall tests ask whether any condition differs and have no single fold-change direction. `PValue` is unadjusted; `FDR` is the Benjamini–Hochberg adjusted p-value within each method, feature type, contrast, cohort and search mode. Do not select the smallest p-value across methods as a combined test. Condition means and detection counts are descriptive; use the named comparisons to assess differences.

**Normalization and low counts.** Standard counts use TMM. Sparse modification tests use TMMwsp and edgeR's small-count quasi-likelihood GLM (`legacy=FALSE`, robust moderation); edgeRun and limma-voom are complementary methods. Filtering is optional, with saved thresholds in the model records. No imputation or pseudocounts are added to the count-model inputs. Modification counts do not measure occupancy, localization-filtered sites or protein-adjusted modification changes.

**Matrices.** Counts are grouped by cohort and search mode rather than combined into a matrix with blank columns from other cohorts. Normalization is separate for each evidence type. Z-scores standardize rows of log2 CPM; constant rows are zero and empty libraries are missing. They are visualization values, not statistical tests. DDA and DDA+ must not be summed or treated as biological replicates. Independent replication and comparable sampling are assumed; batch and subject effects are not inferred.

## Optional analyses

**Gene counts and comparisons.** Select `optional_analyses.gene_counts=true` in Optional analyses → Gene counts. `gene_counts/counts/<cohort>/<mode>/` contains raw, normalized and z-score matrices; `gene_counts/comparisons/` contains named condition and overall tests, including selected gene-modification tests. These use the DEPs methods, reference condition, FDR and count filters; testing requires DEPs enabled. Shared-gene peptides are excluded from this aggregate only. Protein counts and their separately corrected statistics are unaffected. The same Gene counts option includes gene-level drug-interaction tests when that analysis is enabled; those tests use its own methods and FDR. Existing studies preserve gene analyses; new studies default off.

**Peptide and isoform evidence.** Peptides are remapped against all searched FASTA targets, with I/L indistinguishable. Shared sequences remain ambiguous. Isoform specificity requires a uniquely matching accession with alternative isoforms in the searched FASTA. The evidence tables include coverage, detection by condition and unique-peptide counts; they do not establish full proteoforms or isoform-identification FDR.

**Statistical comparisons.** Selected peptide spectral counts and isoform-discriminating counts are compared with their own methods, reference/all-pairs design, FDR and filtering settings. edgeR's small-count model and optional limma-voom use TMMwsp, retaining full-library offsets through filtering. Canonical-only or ambiguous isoform assignments are excluded with reasons in the audit directory. The FDR concerns condition differences; it does not establish isoform identity or switching.

**Drug interaction.** `interaction_log2FC` is combination minus drug A minus drug B plus control on a log2 scale. The null is multiplicative on the original abundance scale. Its sign describes proteomic nonadditivity, not beneficial synergy. Tables record `dose_time_matching_confirmed` and `analysis_interpretation`; unconfirmed matching is an exploratory combination contrast. This annotation does not change the test. Methods have separate results and FDR: edgeR quasi-likelihood, limma-voom, and/or negative-binomial Monte Carlo simulation.

Limma interaction tables include moderated t-statistics and 95% effect confidence intervals. Monte Carlo tables include the likelihood-ratio statistic, empirical p-value, FDR, exceedance/valid simulation counts, seed, minimum p-value, simulation standard error and exact binomial interval for the conditional tail probability. Dispersion and normalization stay fixed during simulation, so calibration is approximate. Its interval describes simulation uncertainty, not biological effect uncertainty. The p-value floor `1/(simulations+1)` can limit BH discoveries. Voom weights, null fits, designs, offsets and simulation settings are supporting model records under `analysis_details/reporting/result_audit/drug_interaction/`.

**HLA-enriched samples.** Dedicated reports contain peptide evidence, sample summaries, lengths and stripped sequences for external binding prediction. Identification applies 1% PSM FDR without a 1% peptide/protein FDR claim; this is separate from differential-test FDR. MBR is off. Differential tests use peptide counts; supporting gene/protein annotations do not measure source-protein abundance. No binding predictor or allele assignment is run.

## Status and provenance

**No significant results** means tests ran and no features passed FDR. **Insufficient evidence** means tests could not be estimated. **No eligible evidence** means there are no eligible rows. Full tables preserve reported nonsignificant/untestable features; untestable rows have blank p-values/FDR and an explicit status. Empty significant subsets are valid and remain downloadable. Counts excluded from main fits remain in count matrices and filter reports.

Final result tables are copied without changing statistics. `--continue` reuses this study’s matching completed Nextflow tasks and verified completed FragPipe searches; preserve `work/`, `.nextflow/` and the output directory. Search checkpoints and provenance remain with completed search files. `analysis_details/reporting/result_audit/file_manifest.csv` records source/destination paths and verified SHA-256 hashes. The audit directory holds run settings, R session details, model inputs, exclusion/filter reports and normalization factors. Historical indices/settings there describe the original computation; use the final `results_index.csv` for current result locations.

Open-search discovery is no longer part of the pipeline or final results. Older configurations drop that option and mass-shift comparisons on loading. An existing completed run can be reformatted without searching or fitting models:

```bash
.runtime/env/bin/python pipeline_lib/result_layout.py /path/to/output/final_results --preserve-original
```

The prior layout, including retired discovery results, is preserved under `analysis_details/previous_result_layout_*`. Future runs generate the same curated final layout automatically.
