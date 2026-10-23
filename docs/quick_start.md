# Quick start: real ubiquitin-affinity proteomics

This small download contains **12 real LFQ samples**, with three experimental replicates in each of four conditions: `Ctrl`, `Ubi1`, `Ubi4` and `Ubi6`. These are ubiquitin-affinity enrichment measurements from [Zhang et al., Molecular Cell (2017)](https://doi.org/10.1016/j.molcel.2017.01.004), distributed as `UbiLength` by the authors' [DEP package](https://github.com/arnesmits/DEP). They measure enrichment of interactors, not whole-cell gene expression.

The example contains existing quantification, so it needs no raw acquisitions, FASTA or registered search engines. It runs the real matrix, QC and statistical pipeline.

## Install and download

On Linux with Python 3.11 or newer:

```bash
pip install --upgrade fun_proteomics
fun_proteomics install --conda
fun_proteomics doctor

wget https://github.com/cfarkas/fun_proteomics/releases/download/v0.9.0/ubiquitin-lfq-example.tar.gz
wget https://github.com/cfarkas/fun_proteomics/releases/download/v0.9.0/ubiquitin-lfq-example.tar.gz.sha256
sha256sum -c ubiquitin-lfq-example.tar.gz.sha256
tar -xzf ubiquitin-lfq-example.tar.gz
```

The archive is smaller than 1 MB. It contains `proteins.csv` / `.tsv`, `samples.csv` / `.tsv`, original upstream data, a reproducible export script, checksums, citation and the upstream Artistic-2.0 license. The measured LFQ values are unchanged. Reverse, contaminant and site-only flags are retained and filtered by the importer.

<details>
<summary>Preview the downloaded input files</summary>

**`samples.csv`** — four selected rows; the full file contains all 12 samples:

```csv
sample_id,condition,replicate
Ctrl_1,Ctrl,1
Ubi1_1,Ubi1,1
Ubi4_1,Ubi4,1
Ubi6_1,Ubi6,1
```

**`proteins.csv`** — three actual rows and selected columns, reordered for this preview:

```csv
Protein.Group,Genes,Ctrl_1,Ubi1_1,Ubi4_1,Ubi6_1
Q99798;A2A274,ACO2,25905000.0,26919000.0,26040000.0,16447000.0
A6NDG6;H3BV17,PGP,186950000.0,188320000.0,214240000.0,212640000.0
B4E241;P84103,SFRS3;SRSF3,261820000.0,245840000.0,367490000.0,244180000.0
```

Use the full downloaded files for analysis. They include three replicates per condition, more protein rows, annotations and filtering flags. The matching TSV files contain the same tables with tabs instead of commas.

</details>

Choose **A, B or C** below. Each method uses the same inputs and writes to `ubiquitin-study/`; run only one method for a fresh project. Use a new project folder to try a different design.

## A) Interactive settings

The video shows the complete browser workflow with this download. Start it once with `fun_proteomics setup`, then make the choices in the page.

https://github.com/user-attachments/assets/bcd925bb-c021-4ccb-92fd-8ac3f4603049

Choose a new **output/project folder**, such as `ubiquitin-study`. Browse to the extracted `proteins.csv` and `samples.csv`.

1. Keep the terminal open and click **Discover samples**.
2. The sample table prefills four condition columns. Review the three cards in each column. Drag cards to change assignments, or to **Unassigned** to exclude them.
3. Keep all-versus-all comparisons and confirm independent experimental replicates for this tutorial. No donor pairing is supplied by the source design.
4. Keep additional normalization set to **none**; these are already LFQ quantities. Missing values are not imputed.
5. Click **Save configuration and check**, then **Run / continue analysis**.
6. Open **Matrices · CSV & TSV** and **Results**.

To practice assigning all cards yourself, leave the sample-table field empty, discover the columns, and create `Ctrl`, `Ubi1`, `Ubi4` and `Ubi6` conditions. Assign the corresponding `_1`, `_2` and `_3` sample cards to each condition.

## B) Terminal questions

```bash
fun_proteomics setup --terminal --project "$PWD/ubiquitin-study" --run
```

Choose **Quantification matrix**. Enter the paths to `ubiquitin-lfq-example/proteins.csv` and `ubiquitin-lfq-example/samples.csv`. Leave FASTA and service benchmark inputs empty. Keep differential analysis enabled, confirm the independent replicate design, and use no extra normalization. Review the output folder before launch.

## C) Full command-line options

```bash
fun_proteomics setup --non-interactive --project "$PWD/ubiquitin-study" \
  --mode matrix --matrix "$PWD/ubiquitin-lfq-example/proteins.csv" \
  --samples "$PWD/ubiquitin-lfq-example/samples.csv" \
  --independent-samples --normalization none \
  --cpus 4 --memory-gib 8 --run
```

This uses all six condition pairs, a minimum of two observed replicates per group, and BH FDR 0.05 separately for each contrast and feature level. This tutorial uses the source's condition-only replicate design. It demonstrates this pipeline's analysis; its missing-value handling and statistical method differ from the original publication.

For all input, pairing, service and raw-data flags, see [command-line options](command_line.md).

## Find results and continue

Open these files:

```text
ubiquitin-study/START_HERE.html
ubiquitin-study/MATRICES.html
```

The checked example exports **2,902 protein-group rows and 2,696 unambiguous gene rows**, with **six protein and six gene comparison tables**. Rows with insufficient observations remain untested. These are matrix row counts, not a claim that every row is positively quantified in every sample.

Detailed files are under `ubiquitin-study/results/final_results/`: `matrices/`, `quality_control/` and `differential_expression/`. Every matrix and result table has matching CSV and TSV versions. The example has no separate service measurement and therefore does not create a service benchmark.

```bash
fun_proteomics --continue "$PWD/ubiquitin-study"
```

A completed unchanged delivery is verified and reused. [Watch the interactive walkthroughs](https://cfarkas.github.io/fun_proteomics/#videos) · [Understand the outputs](outputs.md).

For a tiny offline, synthetic format demonstration, `fun_proteomics demo --project "$PWD/synthetic-study"` remains available. Its measurements are invented and are separate from this public-data example.
