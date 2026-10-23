# Advanced workflow reference

For routine analysis, start with [setup](setup.md), [input previews](input_files.md), or the [single BAF chapter](raw_baf.md). This page covers specialist native configurations.

A Nextflow pipeline with browser-based setup for native Bruker DDA-PASEF and DIA-PASEF acquisitions. The DDA branch uses FragPipe for DDA/DDA+ searches, protein and peptide count tables, modification counts, and optional differential analyses. The DIA branch uses DIA-NN for identification and quantification.

Raw PASEF analyses start with `nextflow run main.nf`. Existing count or quantification tables can also be analyzed directly using the scripts below. Use `--interactively` to configure either acquisition type. Direct DDA runs use `--config`; direct DIA runs use `--dia_config`.

The workflow is independent of cell line, treatment, and organism. Supply your own acquisitions, sample assignments, and matching FASTA. No study data, results, registered software, or machine-specific configuration is included.

## Install and launch from any directory

```bash
python -m pip install fun_proteomics
python -m fun_proteomics --continue /data/my_study
```

`fun_proteomics --continue /data/my_study` is the equivalent installed command. It detects saved Nextflow runs or completed DIA searches, validates reusable results and preserves prior attempts. Completed native results are verified and reused; incomplete native reporting creates a new generation under `continuations/`. Use `--outdir /data/new_delivery` for an explicit new destination. Nextflow runs retain their recorded work directory and session. See [PACKAGE_README.md](../PACKAGE_README.md) for dependencies and detailed continuation behavior.

```bash
fun_proteomics run --config /data/dda.json --outdir /data/results
fun_proteomics run --dia-config /data/dia.json --outdir /data/results_dia
nextflow run main.nf --continue /data/my_study
```

**`dda.json`** — selected fields; start from the complete examples/run.json template:

```json
{
  "outdir": "results/my_study",
  "modes": [
    "DDA"
  ],
  "resources": {
    "cpus": 8,
    "memory_gib": 32
  },
  "database": {
    "fasta": "/data/reference/search.fasta",
    "contains_decoys": true
  }
}
```

**`dia.json`** — selected fields; start from the complete examples/diapasef.json template:

```json
{
  "acquisition": "DIA-PASEF",
  "samplesheet": "diapasef_samples.tsv",
  "raw_root": "/data/raw",
  "canonical_fasta": "/data/targets.fasta",
  "outdir": "results/dia"
}
```

The Python package includes workflow code and reporting assets. It excludes study data, results, credentials and registered binaries. Raw PASEF workflows use Nextflow; conventional [raw BAF](raw_baf.md), existing quantification and verified-search reporting run in native Python. The project is distributed under the [MIT license](../LICENSE).

## Find matrices and compare all conditions

Every delivered matrix is available in both **CSV and TSV**, with identical identifiers, values and missing cells. `MATRICES.html` has separate links for both formats; the index is also supplied as `MATRICES.csv` and `MATRICES.tsv`.

Open `<outdir>/final_results/MATRICES.html` for direct links to every matrix.
The main DDA count matrix is `final_results/matrices/<cohort>/<dda|ddaplus>/protein_counts.csv`.
The DIA abundance matrix is `final_results/matrices/<sample_type>/protein_intensities.csv`,
with `Genes`, protein names and annotation status beside the sample columns.
DIA quantities are intensities, not spectral counts. Gene intensity rollups are in `gene_intensities.csv`.

For DDA, enable DEPs and choose **No control — compare all pairs**.
For DIA, enable **Differential abundance** in DEPs; **All conditions versus all conditions** is the default.
Four conditions produce six unique comparisons. Positive logFC means the second named condition exceeds the first.
DDA uses count models; DIA uses limma-trend on log2 intensities with no missing-value imputation.

Recompute all-pairs DDA statistics from completed results without repeating searches:

```bash
python -m fun_proteomics counts --results /data/run/final_results \
  --outdir /data/run_all_pairs --methods edgeR limma
```

Use a new output directory. This tests proteins, genes and peptides; use `--features protein peptide`
when gene matrices were not requested. Add `--control Control` for reference-only comparisons.
The existing R dependency installation is reused. See [OUTPUTS.md](../OUTPUTS.md) for matrix definitions.

## Reuse completed DIA searches with integrity checks

Rebuild organized matrices, gene annotations, service benchmarks and all-pairs protein/gene models from a completed run:

```bash
python scripts/reanalyze_dia.py --source-run /data/completed_run \
  --outdir /data/rerun_all_pairs --hash-raw \
  --service-matrix tissue=/data/service.tsv
```

The source must retain its engine manifest and completed Nextflow search task. The command checks task completion, SHA256 equality between published and cached reports, matching FASTA/method/sample assignments, and complete Parquet decoding. Raw components receive size and SQLite checks; `--hash-raw` also hashes every raw byte and records a new baseline. It cannot establish historical raw identity without earlier checksums. A mismatched or missing search is refused; run the native pipeline to regenerate it.

Use `--independent-samples` only after confirming one independent biological sample per row. This uses sample IDs for model replicate labels and preserves the original labels as `source_replicate`; it does not establish biological independence. Models use all unordered condition pairs, no imputation and per-contrast BH FDR. The original searches are read without changes. `audit/` records the checks and code/input hashes, `logs/` holds execution logs, `analysis_details/` holds supporting evidence, and `final_results/` contains the deliverables. Historical phosphorylation evidence, when present, is retained separately without claiming newly fitted site-level comparisons.

## Import service quantification and benchmark detection

A complete synthetic example can be run without raw acquisitions or search engines:

```bash
python -m fun_proteomics quantify --matrix examples/quantification_matrix.csv \
  --samples examples/quantification_samples.csv --service-matrix examples/service_quantification.csv \
  --outdir results/example_quantification --differential
```

**`examples/quantification_matrix.csv`** — first rows from the bundled synthetic example:

```csv
Protein.Group,Genes,C1,C2,A1,A2,B1,B2
P00001,GENE_A,100,110,400,410,100,90
P00002,GENE_B,200,210,190,205,800,810
```

**`examples/quantification_samples.csv`** — first rows from the bundled synthetic example:

```csv
sample_id,condition,replicate
C1,Control,1
C2,Control,2
A1,Treatment A,1
A2,Treatment A,2
B1,Treatment B,1
B2,Treatment B,2
```

**`examples/service_quantification.csv`** — first rows from the bundled synthetic example:

```csv
Protein.Group,Genes,C1,C2,A1,A2,B1,B2
P00001,GENE_A,100,100,400,400,100,100
P00002,GENE_B,200,200,200,200,800,800
```

These invented measurements demonstrate file formats and testing, not biological findings.

Import CSV/TSV matrices from DIA-NN, FragPipe or a provider, including metadata in arbitrary positions.
Excel `.xlsx` input is supported when `openpyxl` is installed (included in the main runtime).
Provide a sample table with `sample_id`; differential testing also requires `condition` and `replicate`.
Use one cohort/sample type per import and independently numbered biological replicates within each condition.

```bash
python -m fun_proteomics quantify --matrix /data/protein_matrix.tsv \
  --samples /data/samples.csv --service-matrix /data/provider.tsv \
  --outdir /data/quantification_review --differential
```

**`protein_matrix.tsv / provider.tsv`** — plain abundance table example; matching sample columns:

```tsv
Protein.Group	Genes	C1	C2	T1	T2
P04637	TP53	120000	135000	260000	245000
P60709	ACTB	480000	505000	490000	510000
```

**`samples.csv`**:

```csv
sample_id,condition,replicate
C1,Control,1
C2,Control,2
T1,Treatment,1
T2,Treatment,2
```

**`targets.fasta`** — format illustration; supply a complete matching database:

```fasta
>sp|DEMO001|EXAMPLE_A GN=GENE_A
MPEPTIDEK
>sp|DEMO002|EXAMPLE_B GN=GENE_B
MPEPTIDER
```

With no `--control`, every condition pair is tested. Omit `--service-matrix` to analyze a provider table by itself.
Use `--fasta /data/targets.fasta` to supplement missing gene symbols from UniProt `GN=` headers.
Add `--normalization median` only when median centering fits the study; the default uses reported quantities unchanged.

For a direct-spectral-count detection benchmark (without intensity statistics):

```bash
python -m fun_proteomics quantify --matrix /data/pipeline_proteins.tsv \
  --samples /data/samples.csv --quantity spectral_count \
  --service-matrix /data/provider_combined_protein.tsv --service-quantity spectral_count \
  --outdir /data/detection_benchmark
```

**`pipeline_proteins.tsv / provider_combined_protein.tsv`** — spectral-count input example:

```tsv
Protein	C1 Spectral Count	C2 Spectral Count	T1 Spectral Count	T2 Spectral Count
P04637	12	15	28	30
P60709	45	50	47	49
```

FragPipe import selects MaxLFQ, intensity or spectral-count columns by their names; `--quantity` overrides selection.
For a generic matrix with bare sample columns, use `--quantity value` and explicitly set `--measure spectral_count`
if it contains counts. Count models use `reanalyze_counts.py` or the DDA workflow.
Protein-group identifiers use `Protein.Group`, `Protein`, `FeatureID` or `Protein ID`; `--id-column` overrides this.
Unknown samples, duplicate groups/columns and invalid measurements fail before being silently combined.

A sample table may include `service_sample` for explicit provider names. Use `--service-samples mapping.csv`
when that mapping should apply only to the benchmark table. Run paths and `.d` suffixes are normalized;
treatment labels are never inferred by fuzzy matching. Record acquisition evidence when correcting a provider label.
**`mapping.csv`** — service aliases for the C1/C2/T1/T2 example:

```csv
sample_id,service_sample
C1,Lab_C1
C2,Lab_C2
T1,Lab_T1
T2,Lab_T2
```

The benchmark exports per-sample and union detection counts, shared groups, pipeline-only groups,
service-only groups, signed differences, percentages and complete membership lists.
The percentage is undefined for a zero service baseline. These are detection comparisons, not evidence of matched FDR or higher accuracy.

## Launch the study interface

From the repository folder:

```bash
nextflow run main.nf --interactively
```

Nextflow and Java must already be available to start this command. On Linux x86_64, the launcher automatically prepares its local Python/Java/Nextflow runtime and the R reporting packages, then opens the study interface. The first launch needs internet and may take several minutes; later launches check and reuse installed dependencies. Setup does not start a spectrum search.

No PATH exports or separate setup commands are required. R is selected automatically (project-local R, system R, then PATH); if unavailable, a local R environment is installed with micromamba. Missing reporting packages are installed into `.runtime/R-library`, without sudo. An explicit `--rscript /path/to/Rscript` overrides detection. `PROTEOMICS_RUNTIME_DIR` or `--runtime_dir` changes the runtime root; the reporting library remains project-local.

If setup fails, the terminal shows the error. Resolve connectivity or installation problems and rerun the same command. For manual preparation, `setup_environment.sh` and `scripts/install_differential_tools.sh` remain available.

For DDA searches, obtain the registered MSFragger, IonQuant and diaTracer versions listed in [tool_versions.json](../assets/tool_versions.json). Enter their ZIP paths in the interface or study JSON. The workflow verifies archive checksums. FragPipe can be downloaded by the pipeline; the optional phosphorylation adapter detects or downloads its pinned FragPipe ZIP automatically. Binaries and registration credentials are not redistributed.

## Interactive setup

```bash
nextflow run main.nf --interactively
```

A fresh checkout opens an empty study. Choose input/output folders, detect acquisitions, assign conditions and cohorts, confirm the acquisition method, and select or download a FASTA. The input and output buttons open the operating system’s folder chooser. Use its New Folder action, or enter a name under **Create a subfolder** and click **Create & use folder**. The newly created directory becomes the output root for the run. DDA creates `final_results/` and `analysis_details/`; DIA uses the output folders described below. Zenity or KDialog and a desktop session are required; the window opens on the analysis computer. Folder paths can also be entered directly. Resources are detected automatically and remain editable. A cohort is a group of acquisitions that should be searched and quantified together.

Finish setup, then click **Run analysis**. Input checks occur during setup and Nextflow validates inputs again before the search. No search starts automatically. Browser-started analyses use plain scrolling logs automatically; no logging flag is needed. **Stop run** is enabled only for an active analysis. Repeated stop requests after completion are harmless, and the interactive Nextflow launcher is excluded from scientific-run detection; the browser and launcher terminal show progress. Runs write directly to their log files and retain their exit status if the interface or launcher terminal closes. Reopen the same interface command to reconnect to the latest run, including its Stop control. Closing the interface only detaches it; use **Stop run** to stop analysis. Use `--ui_port 8766` if another interface already occupies the default port.

The setup server runs independently of the launcher terminal. Ctrl+C or closing that terminal detaches the log display; the browser remains usable. Reopening the same command reuses the server and its current configuration. Server diagnostics are saved in `.runtime/ui_server_<port>.log`. If the server itself stops, every setup tab shows a connection notice. Keep the tab open and reopen the setup command; it reconnects without reloading your fields. Requests that might change files or start/stop a run are not automatically repeated after a lost response.

An interface disconnect is never reported as workflow completion. The terminal retries temporary disconnects and, when available, follows the saved run log and exit status until the analysis finishes. `[INTERFACE CLOSED]` describes only the interface. An analysis is complete only when its own run exits successfully. During quiet tasks, heartbeats show the actual time since the last task output and its last reported operation; they do not prove computation is advancing.

The executable adapters support native Bruker DDA-PASEF and DIA-PASEF `.d` directories containing `analysis.tdf` and `analysis.tdf_bin`. Confirming DIA-PASEF selects the global label-free DIA-NN workflow. DDA supports global label-free analysis, the optional phosphorylation search and dedicated HLA-enriched immunopeptidomics. Glycoproteomics, isobaric labeling, and targeted analysis adapters are not implemented. Acquisition metadata is checked again before launching.

**Save draft** keeps incomplete setup, sample names, assignments and current form fields in a new local version without starting an analysis. DIA preserves saved DDA modes, tools, differential and optional-analysis settings, but those settings are inactive during DIA runs.

## DIA-PASEF through main.nf

The DIA-NN adapter uses native Nextflow tasks in the same `main.nf` entry point.
It runs DIA-NN 2.7.0, builds one empirical library per sample type, quantifies
each type separately, and combines counts and reports. Existing DDA configurations
continue to use `--config` and FragPipe.

For browser setup, run `nextflow run main.nf --interactively`, detect the acquisitions,
and confirm **DIA-PASEF / global label-free**. Review the pinned method and supported
outputs before running. Each cohort becomes one DIA sample type with an independent
library. Sample names, conditions and replicates are retained. Unless explicit library
references were supplied, the first 20 acquisitions per cohort, sorted by raw path,
are selected (all samples when fewer than 20). Input checking previews the exact names.

The browser reuses a saved or downloaded FASTA only after verifying its SHA256.
It writes a separate target FASTA, removes only correctly paired `rev_` decoys,
and retains contaminant targets; inconsistent decoy metadata is rejected. An optional
DIA target-only FASTA can be selected instead. Generated DIA JSON, TSV and FASTA files
are stored under `ui_configs/dia/` with stable paths for identical inputs. The editable
browser study remains separate. Python 3.12 is detected automatically, and the default
installation cache is `.runtime/dia-cache`.

```bash
export NXF_VER=26.04.2
nextflow run main.nf --dia_config study.json --validate_only
nextflow run main.nf --dia_config study.json -work-dir work/dia -ansi-log false
nextflow run main.nf --dia_config study.json -work-dir work/dia -resume -ansi-log false
```

Start from [examples/diapasef.json](../examples/diapasef.json) and
[examples/diapasef_samples.tsv](../examples/diapasef_samples.tsv). Every path in the
DIA JSON is relative to that JSON file, except a CLI `--outdir`, which is relative
to the launch directory. Raw paths in the sample table are relative to `raw_root`.
The sample table requires `sample_id`, `sample_type`, `raw_path` and
`library_reference`. IDs must be unique within type; the same subject ID may
occur in multiple types. At least one library reference is required per type.
Condition/replicate and other additional columns are retained during preparation.

**`diapasef_samples.tsv`** — first sample type from the native DIA example:

```tsv
sample_id	sample_type	condition	replicate	raw_path	library_reference
001	serum	Control	1	serum/sample_001.d	true
002	serum	Control	2	serum/sample_002.d	false
003	serum	Treatment	1	serum/sample_003.d	true
004	serum	Treatment	2	serum/sample_004.d	false
```

**`study.json`** — native DIA excerpt; copy the complete template before editing:

```json
{
  "acquisition": "DIA-PASEF",
  "samplesheet": "diapasef_samples.tsv",
  "raw_root": "/data/raw",
  "canonical_fasta": "/data/targets.fasta",
  "outdir": "results/dia"
}
```

`--threads` (or `--cpus`), `--memory`, and `--outdir` override JSON resources and
output. Defaults are 32 threads and 96 GiB per engine task; the existing local
queue processes one task at a time. These memory declarations schedule resources,
not a Linux memory cap. A system Python 3.12 and Java/Nextflow are required.
The JSON `python` field selects the interpreter for the pinned private
environment. `PROTEOMICS_DIA_PYTHON` selects the launch-side config checker if
the command `python3` is unsuitable. This branch does not require the DDA tools
or R reporting runtime unless intensity differential analysis is enabled. That option requires R with limma and jsonlite; the configuration is checked before searching.

`install_cache` points to an existing directory that may be empty. An existing
DIA-NN archive and `wheels/` directory allow offline installation; otherwise the
environment task downloads pinned dependencies. No binaries are packaged with
the source. `canonical_fasta` supplies the target proteome and `contaminant_fasta`
optionally supplies contaminants. An optional `service_matrices` map benchmarks
each type against a supplied protein matrix; it is not used to select references.

The branch rejects DDA, mixed acquisitions, incomplete raws, blanks, duplicate
samples and missing service columns before searching. `--validate_only` also
installs/checks the pinned engine and prepares the FASTA, then stops without
searching. Use native `-resume` for a configured DIA run, or `fun_proteomics --continue /data/run` to discover a saved run. The bare DDA `--continue` flag, `--report_only`, and the DDA
`--dry_run` are not implemented for this branch.

Scientific settings are in [modules/diapasef/method.toml](../modules/diapasef/method.toml).
The inherited method uses fixed carbamidomethylation, no variable modifications,
one missed cleavage, 15 ppm MS1/MS2 tolerances, MS2 evidence and 1% run/global/library
confidence filters. It was developed on LFO; evaluate performance on other sample
types. DIA protein quantities are kept separate from DDA spectral-count statistics.
This adapter exports protein/gene quantification and optional intensity differential results. DDA HLA, phosphorylation and interaction reports remain separate.

DIA results use `final_results/START_HERE.html` and `final_results/MATRICES.html`.
`matrices/<type>/` contains annotated protein and unambiguous gene intensities;
`differential_expression/<type>/<protein|gene>/<comparison>/` contains full and significant tables;
`benchmark/<type>/` contains service detection comparisons. Gene assignments are retained from the report,
supplemented from the searched FASTA, and unresolved/ambiguous assignments are reported explicitly.
Gene rollups sum unambiguous protein-group quantities; they are not RNA expression, DIA-NN Gene.MaxLFQ,
or independently FDR-controlled gene identifications.
Native search output, libraries, prediction and intermediate group evidence are under `analysis_details/`;
input manifests are in `prepared/` and logs in `logs/`. Keep `work/` and `.nextflow/` for resume.

Native DIA JSON accepts `"differential": {"enabled": true, "comparisons": "all_pairwise", "fdr": 0.05,
"min_replicates": 2, "normalization": "none", "gene_level": true}`. The browser keeps these settings in
`dia.differential`. An optional `rscript` selects the R executable. Condition and replicate columns are required
when enabled. Each contrast uses a separate limma-trend fit and BH correction, with insufficient observations
retained as untested. Repeated-measures or batch-adjusted designs require a tailored model.
Native `service_matrices` (browser: `dia.service_matrices`) maps sample types to provider tables.
Older results remain readable; import an existing matrix with `analyze_quantification.py` to add comparisons
without searching again.

```bash
python3 -m unittest discover -s tests -p test_diapasef_config.py -v
python3 tests/run_diapasef_smoke.py
```

The smoke test uses synthetic acquisitions and does not produce scientific results.

## HLA-enriched immunopeptidomics

**Search modes in Study settings:** DDA and DDA+ are separate searches of the same DDA-PASEF acquisitions. Choose **DDA+ and DDA** to generate both sets of outputs; selecting one generates only that mode. Settings and Review explicitly state which results will be generated. Each mode adds a spectrum search; report-only processing cannot produce a mode that was never searched. Saved configurations retain their modes. The HLA preset requires DDA only; after switching back to ordinary proteomics, review the mode selection explicitly.

Drug-condition selections are retained per cohort when samples are rediscovered. The Statistics reference applies to peptide/isoform comparisons; drug interaction uses its own control, drug A, drug B and combination assignments. Missing or duplicate assignments open the Drugs subtab and focus the field to correct.

Step 6, **Optional analyses**, follows DEPs. In its **Immunopeptidomics** subtab, answer **HLA-enriched samples? Yes / No**. Choose Yes only for HLA-enriched, undigested peptide samples and select the enriched class (I or II). No is the default. It uses the bundled [HLA workflow](../assets/fragpipe_hla.workflow), derived from the pinned FragPipe 24.0 Nonspecific-HLA preset, with search settings adapted from [Chen et al., Benchmarking Software for DDA-PASEF Immunopeptidomics](https://doi.org/10.1016/j.mcpro.2025.101492).

The preset searches DDA-PASEF in DDA mode, with nonspecific digestion, 7–25 residues, charge 1–4, mass 200–5000 Da, precursor tolerance ±15 ppm, and fragment tolerance 0.03 Da. It has no fixed modifications; variable modifications are oxidation (M), protein N-terminal acetylation and cysteinylation (C), with at most three modifications per peptide. MSBooster/Percolator and IonQuant are enabled; match-between-runs is disabled. The ordinary phosphorylation selection is inactive in this mode.

Identification uses **1% PSM FDR**, with peptide/ion/protein thresholds at 100%, following the benchmark's main search comparison. This does not establish 1% peptide or protein FDR. Peptide spectral counts are built directly from filtered PSMs, independently of protein inference. The optional organism filter uses FASTA `OS=` annotations. Shared mappings remain visible and do not multiply PSM counts.

DEPs methods, comparisons and optional low-count filtering apply to **peptides only** for HLA searches. Protein/gene/modification matrices remain supporting source annotations, not measurements of source-protein expression. Open `final_results/START_HERE.html` for peptide counts, sample summaries, peptide-length distributions and sequence exports for external binding prediction. Class I exports 8–13 residue sequences; class II exports 9–22. No binding predictor, allele assignment, motif deconvolution or immunogenicity inference is included.

CLI configurations use `proteomics.analysis: "immunopeptidomics"`, confirmed acquisition `DDA-PASEF`, `modes: ["DDA"]`, and `immunopeptidomics: {"hla_class": "I"}` (or `"II"`). Omit or disable `differential.modifications`. The template is tested through configuration, command generation, synthetic PSM reporting and statistical regressions; a representative HLA-enriched acquisition has not been searched for validation. This adapts the paper's settings to current pinned tools rather than reproducing its software-version benchmark. Ordinary digested proteomics acquisitions must not be interpreted as HLA presentation data.

## Optional analyses

These analyses apply to the DDA branch. DIA exports identification and quantification
without the DDA phosphorylation, HLA, differential or drug-interaction reports.

Step 6 contains **Immunopeptidomics**, **Peptides & isoforms**, **Drugs**, and **Statistics** subtabs. Every analysis is optional. HLA enrichment is independent of drug treatment: the peptide, isoform and interaction options work with either global or HLA-enriched data, using the appropriate evidence type.

- **Peptide-level and isoform analysis** exports peptide evidence, detection summaries by condition, sequence coverage and unique-peptide spectral counts. Peptides are mapped against every target in the selected FASTA with I/L treated as indistinguishable. Shared mappings are retained. An isoform-discriminating peptide must match one target whose UniProt accession family includes alternatives. Select an isoform-inclusive FASTA in Database & tools before searching; canonical-only databases cannot establish isoform specificity. This is sequence evidence, not full proteoform reconstruction or isoform-level FDR. HLA results describe source sequences, not source-protein expression.
- **Drugs** first asks **Were drugs applied? Yes / No**. Yes reveals **Are you testing synergy or single drugs?** Single drugs leaves combination testing inactive; use DEPs for single-drug comparisons. Synergy reveals an optional **Drug-combination interaction** analysis. Choose a cohort and four different conditions: control, drug A alone, drug B alone and combination. Confirming matching single-drug doses and treatment times is an optional annotation and never blocks a run. When unchecked or omitted, the results are labeled an exploratory combination contrast because dose or timing differences can contribute to the contrast. At least two independently numbered biological replicates per condition are required. Select edgeR, limma-voom and/or Monte Carlo to test `combination − A − B + control` on a log2 scale. This runs independently of DEPs, testing proteins, genes and peptides for global proteomics, or only peptides for HLA enrichment. It describes proteomic responses at the observed doses/time, not cell-viability synergy or a Bliss/Loewe/ZIP score.

Drug-interaction methods share the four-condition design and TMM normalization but have separate result tables and BH FDR families for each method, cohort, search mode and feature type:

| Method | Test and additional output |
| --- | --- |
| `edgeR` (default for existing configurations) | Negative-binomial quasi-likelihood test of the interaction contrast. |
| `limma` | Voom precision weights and empirical Bayes variance moderation, followed by a moderated t-test. Exports 95% confidence intervals, voom log2 CPM and weights. This is not a fully Bayesian posterior; see the [limma user guide](https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf). |
| `monte_carlo` | Negative-binomial parametric bootstrap from a fitted model containing both single-drug main effects and zero interaction. Refit the full and null coefficient models in every simulation and compare likelihood-ratio deviances. Empirical p-value is `(exceedances + 1) / (simulations + 1)`, followed by BH correction. Both positive and negative departures are tested. |

Monte Carlo settings are **simulations per feature** (default 9,999; range 99–1,000,000), **random seed** (default 1) and **workers** (default 1; up to 8, capped by study threads). Reproducible RNG streams preserve results across worker counts on the same inputs and software. The method holds fitted dispersions and normalization offsets fixed: it is an approximate plug-in bootstrap, not a treatment-label permutation test or a fully Bayesian analysis. Simulation standard errors and exact binomial intervals describe Monte Carlo uncertainty in the conditional tail probability, not biological effect uncertainty. The minimum p-value is `1/(simulations+1)`; 9,999 simulations gives 0.0001, which can limit BH discoveries in large families. Increasing simulations improves resolution and costs time; it does not increase biological replication. Do not treat the smallest p-value across methods as a combined test.

**Statistics** adds independent comparisons for selected peptide spectral counts and/or isoform-discriminating peptide counts. Enable comparisons, select edgeR and/or limma-voom, choose a reference condition or all pairs, and set BH FDR and optional count filtering. These controls work with DEPs off. Each feature type, method, contrast, cohort and search mode has a separate test family. The engine uses TMMwsp and the edgeR small-count quasi-likelihood model, with an overall condition test and named pairwise tests; limma-voom is an optional comparison method. At least two conditions and two independent replicates per condition are required. There is no inferred batch or subject adjustment. Drug interaction retains its separate four-condition model and FDR.

Isoform comparisons include only counts supported by uniquely mapping peptides whose accession family includes alternatives in the searched FASTA; other accessions are excluded with reasons. These do not establish isoform switching. HLA comparisons describe peptide evidence, not source-protein expression. Differential-test FDR is separate from identification FDR.

**Low coverage does not fail the statistical analyses.** Nonzero features remain eligible unless the user explicitly enables the corresponding low-count filter. Full tables retain nonsignificant results. Optional-statistics tables also retain all-zero or filtered features with explicit status and blank p-values. When a model cannot be estimated, the output retains the features with an insufficient-evidence `test_status`, blank p-values/FDR and `significant=false`; this is not evidence of no effect. Empty isoform-evidence tables are valid outputs. Invalid inputs and genuine search execution failures remain errors.

Open `final_results/START_HERE.html` for an offline dashboard with separate **Main results**, **Optional analyses**, **Counts & study QC**, and interpretation views. Filter test tables by analysis, search mode, method and result status. Every test links to its full table and significant subset; zero discoveries and insufficient evidence are distinguished. `results_index.csv` lists all generated tests without pooling discoveries. `--report_only` regenerates statistical analyses from completed searches.

Open-search discovery has been removed. Older configurations automatically drop its setting and mass-shift comparisons while preserving selected peptide and isoform statistics. To apply the new layout to completed results without recomputing statistics, run `.runtime/env/bin/python pipeline_lib/result_layout.py <outdir>/final_results --preserve-original`. Original result tables are archived under `analysis_details/previous_result_layout_*`; the delivered files are copied with verified SHA-256 checksums. This same layout is generated automatically for future runs.

Native DDA configuration excerpt; merge this optional block into your complete configuration:

```json
{
  "optional_analyses": {
    "hla_enriched": false,
    "peptide_isoforms": true,
    "statistics": {
      "enabled": true,
      "features": [
        "peptide",
        "isoform"
      ],
      "methods": [
        "edgeR",
        "limma"
      ],
      "comparisons": "control",
      "control": "Control",
      "fdr": 0.05,
      "filter_low_counts": false
    },
    "drugs": {
      "applied": true,
      "study_type": "synergy",
      "interaction": {
        "enabled": true,
        "cohort": "study",
        "control": "Control",
        "drug_a": "Drug A",
        "drug_b": "Drug B",
        "combination": "Combination",
        "matched_design": true,
        "fdr": 0.05,
        "methods": [
          "edgeR",
          "limma",
          "monte_carlo"
        ],
        "monte_carlo": {
          "simulations": 9999,
          "seed": 1,
          "workers": 2
        }
      }
    }
  }
}
```

The condition/cohort names must match the sample assignments. For single drugs, use `{"applied": true, "study_type": "single"}`; for untreated studies, `{"applied": false}`. Inactive drug choices are omitted when saving the interface. Older configurations without the optional block continue to work.

## Command-line setup

```bash
cp examples/run.json config.local.json
# Edit sample paths, conditions, tools, resources, and FASTA settings.
sha256sum /data/reference/search.fasta
nextflow run main.nf --config config.local.json --validate_only
nextflow run main.nf --config config.local.json --verbosity 3
```

The example contains synthetic names and placeholder paths/checksum. Replace them before running. The organism filter defaults to `null` (all target organisms); choose a scientific species name to restrict reporting. Use a FASTA appropriate to the experiment regardless of this reporting filter. Local FASTA input requires a SHA256; the browser download option builds the database and fills this automatically. For a target-only FASTA, provide the contaminant archive and checksum as described by the schema.

For terminal runs, use `nextflow run main.nf --config config.local.json -ansi-log false`. The workflow prepares its runtime automatically. Advanced logging overrides are documented in `nextflow run main.nf --help`.

**`config.local.json`** — native DDA excerpt; replace paths, samples and database checksum in the complete template:

```json
{
  "outdir": "results/my_study",
  "modes": [
    "DDA"
  ],
  "resources": {
    "cpus": 8,
    "memory_gib": 32
  },
  "database": {
    "fasta": "/data/reference/search.fasta",
    "contains_decoys": true
  }
}
```

**`search.fasta`** — format illustration only; use the complete validated target/contaminant/decoy database:

```fasta
>sp|DEMO001|EXAMPLE_A GN=GENE_A
MPEPTIDEK
>rev_sp|DEMO001|EXAMPLE_A GN=GENE_A
KEDITPEPM
>contam_DEMO
MPEPTIDER
>rev_contam_DEMO
REDITPEPM
```

JSON paths are relative to the JSON file. CLI `--outdir` is relative to the launch directory. Every supported field is described in [config.schema.json](../assets/config.schema.json); `nextflow run main.nf --help` lists pipeline flags. GPU execution is not exposed.

```bash
# Override resources and output location.
nextflow run main.nf --config config.local.json --outdir results/study \
  --threads 16 --memory 64 --verbosity 4
# Continue this study. Preserve work/ and .nextflow/.
nextflow run main.nf --config config.local.json --continue
# Rebuild tables and differential results without another spectrum search.
nextflow run main.nf --config config.local.json --report_only
```

`--continue` selects the latest real analysis for the same output directory, including a completed run, and ignores setup, validation and preview sessions. Matching Nextflow tasks are cached. Completed FragPipe searches also survive reporting-only changes when their acquisition audits, database, effective workflow settings and tool versions match. New search checkpoints verify result checksums. Legacy reuse requires retained work-directory plans and provenance. Missing or changed steps run again; an interrupted search may restart. Native `-resume [RUN_NAME]` remains supported. The browser shows **Continue completed work** in both Review and Run & progress. Its checked state is retained between steps and applies to the Run button and saved terminal command; uncheck it for a new analysis.

The output subfolder entered in setup is created and selected when you continue; clicking **Create & use folder** separately is optional.

## Results and differential analysis

Start by opening `final_results/START_HERE.html` in your browser. Protein analyses are primary. Main statistical results are in `main_results/differential_expression/`, with methods together under each cohort, search mode, comparison and feature type. Raw and normalized protein/peptide matrices are in `main_results/counts/`. Optional gene counts, comparisons, normalized matrices and z-scores are in `optional_analyses/gene_counts/`; sequence evidence and optional statistical results are in `optional_analyses/`. Detailed model inputs, settings and diagnostics remain under `analysis_details/reporting/`. DDA and DDA+ remain separate searches, never biological replicates. [OUTPUTS.md](../OUTPUTS.md) describes all files and count definitions.

Differential analysis is optional and disabled for a new study. Uncheck **Filter low-count features before testing** to disable count prefiltering: inactive thresholds are omitted from the saved configuration, and legacy numeric zeros are ignored. Re-enabling the filter restores usable defaults for invalid fields. All-zero features are still excluded because no count model can test them. In step 5 (**DEPs**), leave **Detect differential proteins (primary analysis) and supporting peptides** unchecked to generate counts and normalized matrices without DEPs, or enable it and select edgeR, edgeRun and/or limma. The normal run includes the selected differential methods. Use CLI `--report_only` to rebuild results from completed count tables. At least two conditions and two independent replicates per condition are required. No paired-subject or batch adjustment is inferred.

- edgeR fits a joint condition GLM for every differential analysis, including two conditions, and produces an overall condition test.
- limma uses voom and a joint condition linear model.
- edgeRun adds pairwise unconditional exact tests. Its `edgeRun_upper` default/minimum is **100000**, in increments of 100: a deterministic summation bound, not simulation iterations.
- Choose any condition as the control/reference in DEPs, or select **No control — compare all pairs**. Filenames name both groups, such as `edgeR_CTRL_vs_Treatment.csv`; optional short labels simplify long condition names. Positive logFC means the second group exceeds the first.
- Enable **Also test modification counts** to test distinct modified peptide counts per gene and per protein, separately for each modification type.
- BH-adjusted p-values are calculated separately per method, contrast, cohort and search mode. Default FDR is 0.05.

Protein-group spectral counts are used for differential testing. Gene counts sum peptide spectral counts with unambiguous single-gene assignments, using both Gene and Mapped Genes. Protein, peptide and protein-modification matrices are delivered under `final_results/main_results/counts/<cohort>/<mode>/`. Select **Optional analyses → Gene counts** (`optional_analyses.gene_counts=true`) to include gene and gene-modification matrices and comparisons under `final_results/optional_analyses/gene_counts/`. Their comparisons use the DEPs methods, comparisons and filtering settings with separate model fitting and FDR correction; selecting them does not change protein counts or protein statistics. With DEPs disabled, the selected gene matrices are delivered without comparisons. The same Gene counts option includes gene-level drug-interaction tests when that analysis is enabled; those tests use its own methods and FDR. Existing saved studies retain their gene analyses; new studies default this option off. Supporting model matrices and audit files are retained under `analysis_details/reporting/`. Modification counts do not measure site occupancy and are not filtered by localization probability. Shared/ambiguous gene assignments are reported separately.

Run logs, traces and execution reports are stored under `<outdir>/logs/`. The interface displays the live log in that folder. Direct CLI runs keep the controller log at `.nextflow.log` unless `nextflow -log PATH run ...` is supplied.

Results are grouped by cohort, search mode, comparison and feature type under `final_results/main_results/differential_expression/`, with each method’s full and significant tables together. Sparse modification analyses use TMMwsp normalization and the edgeR v4 small-count quasi-likelihood GLM; edgeRun and limma remain comparison methods. Low-count filtering is **off by default** for all feature types; only all-zero rows are excluded. In step 5, enable **Filter low-count features before testing** to apply configurable thresholds (proteins/genes/peptides: adjusted count 10, total 15; modifications: adjusted count 1, total 3). Saved configurations without `differential.filter_low_counts` also use no low-count filter. Each new run records this choice in its results README and filtering reports. Shared gene/protein assignments are excluded. Detailed definitions and limitations are in [OUTPUTS.md](../OUTPUTS.md).

## Validation and packaging

```bash
python scripts/validate_repo.py
bash tests/run_all.sh
bash tests/run_smoke.sh
python scripts/package_source.py
```

The smoke test uses synthetic acquisitions and stub search outputs; it does not validate scientific search performance. Statistical tests require the R packages above. Packaging uses an explicit source allowlist and produces a local ZIP and checksum manifest. It does not upload anything.

See [LICENSE_NOTES.md](../LICENSE_NOTES.md) for redistribution status. A repository-wide license has not yet been selected; no third-party software license is granted here.

## Tool discovery and repair

Steps 3 and 4 share a tool panel. Detection checks configured paths, the local runtime, home/Downloads, the study directory, `/opt`, and `/usr/local` with bounded recursion. Provide an additional tools folder to search elsewhere. Only archives matching the pinned SHA256 are adopted; installed binaries without a verified ZIP are reported separately.

**Install missing tools** downloads public FragPipe and installs available registered archives into `.runtime/tools/`. **Force reinstall tools** replaces only managed installations, verifies archives again, and redownloads public FragPipe. MSFragger, IonQuant and diaTracer still require their official registration/download; the panel links to those providers. Their registration is never bypassed. The workflow continues to stage tools from the verified ZIPs for each Nextflow search.

Step 3 has a separate **Phosphorylation (optional)** section, including module installation and modification count definitions. Phosphorylation can be selected there before installation. Continue to Database & tools or install the module directly; its FragPipe dependency is resolved automatically. Input validation still requires the module to be installed before an analysis starts.

Differential reports separate `pairwise/` group comparisons from `omnibus/` overall condition tests within each method. Each cohort, search mode and feature type has its own folder, with a `significant/` subfolder. `FDR` is the sole adjusted-p-value column. All tables include mean raw counts, mean normalized CPM and detected/total sample counts for each named condition. Omnibus tables omit ambiguous fold-change columns; use named pairwise tests for direction and significance.

## Preparing a GitHub repository

See [GITHUB_UPLOAD.md](../GITHUB_UPLOAD.md) for source validation, packaging and upload commands.
The repository name is `pasef-proteomics`. Experimental data, results, installed tools
and local study configurations are excluded from the source package.

Intensity output conventions follow the [DIA-NN documentation](https://github.com/vdemichev/DiaNN).
Intensity comparisons use [limma](https://bioconductor.org/packages/release/bioc/html/limma.html);
FragPipe imports follow its [output guide](https://fragpipe.nesvilab.org/docs/tutorial_fragpipe_outputs.html).
