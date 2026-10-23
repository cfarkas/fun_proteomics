# Input files at a glance

Start with the preview that matches your input. These values illustrate formats; the [downloadable quick start](quick_start.md) supplies real measurements. CSV uses commas; TSV uses tabs. Keep sample IDs identical across the files you match. Binary acquisitions, Excel workbooks and search-tool ZIPs are selected as files; Excel previews below show the worksheet cells.

## Protein matrix and sample table

**`samples.csv`**:

```csv
sample_id,condition,replicate
C1,Control,1
C2,Control,2
T1,Treatment,1
T2,Treatment,2
```

**`proteins.csv`**:

```csv
Protein.Group,Genes,C1,C2,T1,T2
P04637,TP53,120000,135000,260000,245000
P60709,ACTB,480000,505000,490000,510000
P02768,ALB,90000,85000,,100000
```

Use linear nonnegative intensities. Blank/NA cells remain missing; zero is not a positive detection. `Genes` can be omitted when a matching annotated FASTA is supplied. These are protein quantities; spectral-count tables are a separate measurement.

## Paired donors

**`samples.csv`** — paired version:

```csv
sample_id,condition,donor
D1_CTL,Control,D1
D1_Treatment,Treatment,D1
D2_CTL,Control,D2
D2_Treatment,Treatment,D2
```

Select `--subject-column donor` for this design. One donor can contribute one biological sample in each condition. Aggregate technical injections before matrix import; [BAF grouping](raw_baf.md) happens during raw quantification.

## Service quantification and worksheet names

**`service.xlsx`** — preview of its Identified proteins worksheet:

```csv
Protein IDs,Gene names,LFQ intensity Lab_C1,LFQ intensity Lab_C2,LFQ intensity Lab_T1,LFQ intensity Lab_T2
P04637,TP53,118000,132000,255000,240000
P60709,ACTB,475000,500000,485000,505000
```

`--service-sheet 'Identified proteins'` selects that worksheet. The sheet name is not a path. CSV/TSV inputs need no worksheet selection.

**`service_samples.csv`** — only when service and pipeline sample names differ:

```csv
sample_id,service_sample
C1,Lab_C1
C2,Lab_C2
T1,Lab_T1
T2,Lab_T2
```

The measurement prefix `LFQ intensity` is recognized automatically; aliases identify the remaining sample name. See the [complete service comparison](benchmarking.md).

## Raw PASEF sample table

**`samples.csv`**:

```csv
sample_id,condition,replicate,raw_path
C1,Control,1,C1.d
C2,Control,2,C2.d
T1,Treatment,1,T1.d
T2,Treatment,2,T2.d
```

These paths are relative to `--raw-root`. Each `.d` folder must contain its original binary acquisition components. For BAF `biological_sample` and donor columns, use the [single BAF chapter](raw_baf.md#1-install-and-gather-inputs).

## Annotated FASTA

**`targets.fasta`** — format illustration, not a search database:

```fasta
>sp|DEMO001|EXAMPLE_A GN=GENE_A
MPEPTIDEK
>sp|DEMO002|EXAMPLE_B GN=GENE_B
MPEPTIDER
```

Use the complete organism database for analysis. `GN=` supplies gene names. DIA uses target-only FASTA; DDA/BAF needs a prepared target/contaminant/decoy search FASTA as illustrated in their guides.

## Saved configuration

**`config/run.json`** — minimal matrix configuration normally written by setup:

```json
{
  "fun_proteomics_project": 1,
  "project": "/data/my-study",
  "mode": "matrix",
  "matrix": "/data/proteins.csv",
  "samples": "/data/samples.csv",
  "independent_samples": true,
  "normalization": "none"
}
```

Setup records additional settings and assignment audits. Review saved choices with `fun_proteomics check`; use a new project after launch if the study design changes. Native raw-workflow JSON has a different schema; the [advanced reference](advanced_reference.md) shows those previews and complete templates.
