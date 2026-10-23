# fun_proteomics

Analyze raw Bruker DDA-PASEF, DIA-PASEF and conventional DDA BAF data, or import a protein quantification table. Get protein and gene matrices, all-versus-all comparisons, quality control, and service detection benchmarks. Every matrix is available as **CSV and TSV**.

## Install

Linux, Python 3.11 or newer:

```bash
pip install --upgrade fun_proteomics
fun_proteomics install --conda
fun_proteomics doctor
```

[Installation and existing runtimes](https://github.com/cfarkas/fun_proteomics/blob/main/docs/installation.md) · [PyPI](https://pypi.org/project/fun-proteomics/) · [Walkthrough videos](https://cfarkas.github.io/fun_proteomics/#videos)

## Quick start with real data

Download a small, published ubiquitin-affinity LFQ dataset: **12 samples, four conditions, three experimental replicates per condition**. The archive includes protein quantities, sample assignments, original source files, checksums and data attribution.

```bash
wget https://github.com/cfarkas/fun_proteomics/releases/download/v0.9.0/ubiquitin-lfq-example.tar.gz
wget https://github.com/cfarkas/fun_proteomics/releases/download/v0.9.0/ubiquitin-lfq-example.tar.gz.sha256
sha256sum -c ubiquitin-lfq-example.tar.gz.sha256
tar -xzf ubiquitin-lfq-example.tar.gz
fun_proteomics setup --project "$PWD/ubiquitin-study" \
  --mode matrix --matrix "$PWD/ubiquitin-lfq-example/proteins.csv" \
  --samples "$PWD/ubiquitin-lfq-example/samples.csv"
```

The browser opens with the dataset selected. Discover the samples, review the four condition columns, confirm the replicate design, then save and run. This analyzes real quantified data; it does not repeat a raw spectrum search. [Step-by-step quick start, terminal and scripted versions](https://github.com/cfarkas/fun_proteomics/blob/main/docs/quick_start.md).

## Set up your analysis

Choose **one** of the following launch methods. **`--project` is your output folder**; `setup --outdir` is an equivalent option. Use a new project for a new analysis.

### A) Interactive settings

```bash
fun_proteomics setup --project /data/my-study
```

Keep the terminal open. In the local browser:

1. Choose a matrix, raw-data folder or previous run, and review the output folder.
2. Click **Discover samples**. Drag sample cards into condition columns. Add or rename conditions; move excluded samples to **Unassigned**. A sample table can prefill assignments.
3. Choose independent samples or paired matrix/BAF comparisons with donor IDs. Keep all-versus-all comparisons; four conditions produce six pairs at each feature level.
4. Click **Save configuration and check**, then **Run / continue analysis**.
5. Open **Results** or **Matrices · CSV & TSV** in the same page.

Use the Browse buttons to select files on the computer running the analysis. [Input formats, technical injections and SSH browser access](https://github.com/cfarkas/fun_proteomics/blob/main/docs/setup.md).

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

For paired matrix/BAF data, add a `donor` column to the sample table and replace `--independent-samples` with `--subject-column donor`. Add `--service-matrix /data/provider.tsv` to benchmark a matching service table. [All command-line options and raw DIA/DDA/BAF examples](https://github.com/cfarkas/fun_proteomics/blob/main/docs/command_line.md).

To check and run an existing saved configuration:

```bash
fun_proteomics check --config /data/my-study/config/run.json
fun_proteomics run --config /data/my-study/config/run.json
```

All three methods save `config/run.json`. Browser sample assignments are also saved as CSV and TSV.

## Results and continuation

Open **`my-study/START_HERE.html`** for results or **`my-study/MATRICES.html`** for direct CSV/TSV downloads. Detailed outputs live in `my-study/results/final_results/`.

```bash
fun_proteomics --continue /data/my-study
```

Completed native results are verified and reused. Native BAF checks converted spectra and complete searches. Raw PASEF workflows resume their Nextflow session. [Integrity checks and recovery](https://github.com/cfarkas/fun_proteomics/blob/main/docs/continuation.md).

## Guides

| Task | Guide |
| --- | --- |
| Try the real LFQ dataset | [Quick start](https://github.com/cfarkas/fun_proteomics/blob/main/docs/quick_start.md) |
| Watch interactive launches | [Videos](https://cfarkas.github.io/fun_proteomics/#videos) · [transcripts](https://github.com/cfarkas/fun_proteomics/blob/main/docs/videos.md) |
| Import service Excel workbooks | [Provider workbooks](https://github.com/cfarkas/fun_proteomics/blob/main/docs/provider_workbooks.md) |
| Re-search raw BAF data | [BAF guide](https://github.com/cfarkas/fun_proteomics/blob/main/docs/raw_baf.md) |
| Compare pipeline and service detections | [Benchmarking](https://github.com/cfarkas/fun_proteomics/blob/main/docs/benchmarking.md) |
| Review the full EV_PPC execution and service comparisons | [Execution report](https://github.com/cfarkas/fun_proteomics/blob/main/docs/validation.md) |
| Analyze your own data | [Input and sample setup](https://github.com/cfarkas/fun_proteomics/blob/main/docs/setup.md) |
| Find matrices and comparisons | [Outputs](https://github.com/cfarkas/fun_proteomics/blob/main/docs/outputs.md) |
| Use advanced native Python or Nextflow settings | [Advanced reference](https://github.com/cfarkas/fun_proteomics/blob/main/docs/advanced_reference.md) |
| Fix an error | [Troubleshooting](https://github.com/cfarkas/fun_proteomics/blob/main/docs/troubleshooting.md) |

MIT licensed. Maintained by [cfarkas](https://github.com/cfarkas). The public example retains its upstream data license and citation.
