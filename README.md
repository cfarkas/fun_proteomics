# fun_proteomics

Analyze raw Bruker DDA-PASEF, DIA-PASEF and conventional DDA BAF data, or import a protein quantification table. Get protein and gene matrices, all-versus-all comparisons, quality control, and service detection benchmarks. Every matrix is available as **CSV and TSV**.

## Interactive example: drag .d samples into conditions

After the one-time [installation](#install), start in your working directory. Replace `/path/to/dir` with the directory you want to work in:

```bash
cd /path/to/dir
```

Invoke the pipeline interactively:

```bash
fun_proteomics setup
```

The setup page opens in your browser. Keep the terminal open, then:

1. In **Discover samples**, choose your data type, then click **Choose folder** (raw `.d` data) or **Choose file** (a protein matrix). Samples appear automatically. The top **Discover samples** button also opens the picker when no input is selected.
2. **Assign samples to conditions.** Drag cards into **Control**, **Treatment**, or groups you add. Each card also has a condition dropdown. Leave excluded samples unassigned; raw-data discovery skips blanks.
3. Click **Next: review & save**, then **Choose folder** for your results. Browse to a parent folder and use **New folder** to create a project by name. Choose the required FASTA and confirm the experimental design. Optional settings stay collapsed.
4. Click **Save and continue**, then **Run / continue analysis**. Open **Results** or **Matrices · CSV & TSV** when the run finishes.

No full paths need to be typed. The picker includes folder shortcuts, breadcrumbs and search. Use **Open saved project** at the top to browse to an earlier project.

The example uses independent biological replicates. Conditions come from your assignments; filenames do not determine them. Keep the raw-data and output folders separate. The input guide below explains the FASTA and registered-tool requirements for each acquisition method.

For DIA and imported intensity comparisons, **BH remains the default**. Optional `fdrtool_pvalue` and `fdrtool_normal` methods are available in advanced setup and with `--fdr-method`. The [FDR-method guide](docs/fdr_methods.md) explains the separate original and recalibrated P-value columns.

## Install

Linux, Python 3.11 or newer:

```bash
pip install --upgrade fun_proteomics
fun_proteomics install --conda
fun_proteomics doctor
```

[Installation and existing runtimes](docs/installation.md) · [PyPI](https://pypi.org/project/fun-proteomics/) · [Walkthrough videos](https://cfarkas.github.io/fun_proteomics/#videos)

## Quick start with real data

Download a small, published ubiquitin-affinity LFQ dataset: **12 samples, four conditions, three experimental replicates per condition**. The archive includes protein quantities, sample assignments, original source files, checksums and data attribution.

```bash
wget https://github.com/cfarkas/fun_proteomics/releases/download/v0.9.0/ubiquitin-lfq-example.tar.gz
wget https://github.com/cfarkas/fun_proteomics/releases/download/v0.9.0/ubiquitin-lfq-example.tar.gz.sha256
sha256sum -c ubiquitin-lfq-example.tar.gz.sha256
tar -xzf ubiquitin-lfq-example.tar.gz
```

Watch the interactive walkthrough below. Start the browser once with `fun_proteomics setup`, then choose the extracted files and a new **output/project folder**, such as `ubiquitin-study`. The video shows file selection, condition assignment, all-pairs analysis, CSV/TSV downloads and continuation.

https://github.com/user-attachments/assets/bcd925bb-c021-4ccb-92fd-8ac3f4603049

This recording uses the downloaded real data and runs the complete matrix analysis. It needs no raw spectrum search. [Step-by-step quick start, terminal and scripted versions](docs/quick_start.md).

## Set up your analysis

Choose **one** of the following launch methods. **`--project` is your output folder**; `setup --outdir` is an equivalent option. Use a new project for a new analysis.

### A) Interactive settings

Follow the embedded video above. Open the browser with `fun_proteomics setup`, then choose your data. Keep the terminal open. In the local browser:

1. Choose a matrix, raw-data folder or previous run with its file/folder picker.
2. Samples are discovered automatically when you select an input. Drag sample cards into condition columns. Add or rename conditions; move excluded samples to **Unassigned**. A sample table can prefill assignments.
3. Click **Next: review & save**, choose or create the output folder, then confirm independent samples or paired matrix/BAF comparisons with donor IDs. Keep all-versus-all comparisons; four conditions produce six pairs at each feature level.
4. Click **Save and continue**, then **Run / continue analysis**.
5. Open **Results** or **Matrices · CSV & TSV** in the same page.

Use **Choose file** and **Choose folder** to browse the computer running the analysis; no full paths need to be typed. [Input formats, technical injections and SSH browser access](docs/setup.md).

### B) Terminal questions

```bash
fun_proteomics setup --terminal --project /data/my-study --run
```

Answer the questions about input files, conditions, replicate design and output folder. This uses the same analysis as the browser and works over SSH. Omit `--run` to save and review the configuration first.

### C) Full command-line options

Supply every input explicitly to run without questions:

```bash
fun_proteomics setup --non-interactive --project /data/my-study \
  --mode matrix --matrix /data/proteins.csv --samples /data/samples.csv \
  --independent-samples --normalization none \
  --cpus 8 --memory-gib 32 --run
```

Example input files (illustrative intensity values):

**`samples.csv`** — one row per sample:

```csv
sample_id,condition
C1,Control
C2,Control
T1,Treatment
T2,Treatment
```

**`proteins.csv`** — one row per protein group:

```csv
Protein.Group,Genes,C1,C2,T1,T2
P04637,TP53,120000,135000,260000,245000
P60709,ACTB,480000,505000,490000,510000
P02768,ALB,90000,85000,,100000
```

Each `sample_id` must exactly match a quantity column in `proteins.csv`. Use non-log-transformed intensities; leave missing measurements blank. `Genes` is optional: include gene symbols here or supply a matching FASTA with `GN=` annotations for gene-level results.

For paired matrix/BAF data, add a `donor` column to the sample table and replace `--independent-samples` with `--subject-column donor`. Add `--service-matrix /data/provider.tsv` to benchmark a matching service table. [All command-line options and raw DIA/DDA/BAF examples](docs/command_line.md).

**`provider.tsv`** — optional service quantities using the same sample names:

```tsv
Protein.Group	Genes	C1	C2	T1	T2
P04637	TP53	120000	135000	260000	245000
P60709	ACTB	480000	505000	490000	510000
```

To check and run an existing saved configuration:

```bash
fun_proteomics check --config /data/my-study/config/run.json
fun_proteomics run --config /data/my-study/config/run.json
```

All three methods save `config/run.json`. Browser sample assignments are also saved as CSV and TSV.

<details>
<summary>What the saved configuration looks like</summary>

**`config/run.json`** — generated excerpt; setup fills the complete file:

```json
{
  "project": "/data/my-study",
  "mode": "matrix",
  "matrix": "/data/proteins.csv",
  "samples": "/data/samples.csv"
}
```

Create and review this file through setup. [Complete configuration example](docs/input_files.md#saved-configuration).

</details>

## Results and continuation

Open **`my-study/START_HERE.html`** for results or **`my-study/MATRICES.html`** for direct CSV/TSV downloads. Detailed outputs live in `my-study/results/final_results/`.

```bash
fun_proteomics --continue /data/my-study
```

Completed native results are verified and reused. Native BAF checks converted spectra and complete searches. Raw PASEF workflows resume their Nextflow session. [Integrity checks and recovery](docs/continuation.md).

## Guides

| Task | Guide |
| --- | --- |
| Try the real LFQ dataset | [Quick start](docs/quick_start.md) |
| Watch interactive launches | [Videos](https://cfarkas.github.io/fun_proteomics/#videos) · [transcripts](docs/videos.md) |
| Import service Excel workbooks | [Provider workbooks](docs/provider_workbooks.md) |
| Re-search raw BAF data | [BAF guide](docs/raw_baf.md) |
| Compare pipeline and service detections | [Benchmarking](docs/benchmarking.md) |
| See file contents before preparing inputs | [Input file previews](docs/input_files.md) |
| Analyze your own data | [Input and sample setup](docs/setup.md) |
| Find matrices and comparisons | [Outputs](docs/outputs.md) |
| Use advanced native Python or Nextflow settings | [Advanced reference](docs/advanced_reference.md) |
| Fix an error | [Troubleshooting](docs/troubleshooting.md) |

MIT licensed. Maintained by [cfarkas](https://github.com/cfarkas). The public example retains its upstream data license and citation.
