# Full command-line options

All launch methods use the same saved configuration. `--project` selects the **output folder**, containing `config/`, `results/`, logs and the result/matrix shortcuts. Use a new project for changed inputs or design.

## A) Interactive settings

```bash
fun_proteomics setup --project /data/my-study
```

Discover samples, drag cards into conditions, review the design, save and run. [Browser instructions](setup.md).

## B) Terminal questions

```bash
fun_proteomics setup --terminal --project /data/my-study --run
```

Omit `--run` to save before launching. [Try the real dataset](quick_start.md).

## C) Full command-line options

### Independent matrix samples

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

### Paired matrix samples and a service benchmark

```bash
fun_proteomics setup --non-interactive --project /data/paired-study \
  --mode matrix --matrix /data/proteins.csv --samples /data/samples.csv \
  --subject-column donor --normalization none \
  --service-matrix /data/provider.xlsx --run
```

Add a `donor` column to `samples.csv`, with the same donor ID for samples from the same person. Use one biological quantity per donor and condition. A sample table can include `service_sample` aliases, or add `--service-samples /data/service_samples.csv` for a separate benchmark mapping. [Workbook formats](provider_workbooks.md).

**`samples.csv`** — paired version; one row per donor and condition:

```csv
sample_id,condition,donor
C1,Control,D1
T1,Treatment,D1
C2,Control,D2
T2,Treatment,D2
```

**`provider.xlsx`** — preview of an Excel quantity sheet:

```csv
Protein IDs,Gene names,LFQ intensity C1,LFQ intensity T1,LFQ intensity C2,LFQ intensity T2
P04637,TP53,120000,260000,135000,245000
P60709,ACTB,480000,490000,505000,510000
```

**`service_samples.csv`** — optional mapping when a provider instead uses Lab_C1/Lab_C2/Lab_T1/Lab_T2:

```csv
sample_id,service_sample
C1,Lab_C1
C2,Lab_C2
T1,Lab_T1
T2,Lab_T2
```

A mapping uses the sample IDs from your own sample table. In the paired example above, service names already match and no mapping is needed.

### Raw BAF

Use the [complete BAF chapter](raw_baf.md) for its embedded video, acquisition table, search FASTA, service table, terminal/scripted commands and continuation.

### Raw DIA-PASEF

```bash
fun_proteomics setup --non-interactive --project /data/dia-study \
  --mode dia --raw-root /data/raw --samples /data/samples.csv \
  --fasta /data/human_targets.fasta --independent-samples \
  --cpus 16 --memory-gib 64 --run
```

**`samples.csv`** — raw PASEF version; paths are relative to the selected raw folder:

```csv
sample_id,condition,replicate,raw_path
C1,Control,1,C1.d
C2,Control,2,C2.d
T1,Treatment,1,T1.d
T2,Treatment,2,T2.d
```

**`human_targets.fasta`** — header/sequence format illustration only:

```fasta
>sp|DEMO001|EXAMPLE_A GN=GENE_A
MPEPTIDEK
>sp|DEMO002|EXAMPLE_B GN=GENE_B
MPEPTIDER
```

Use a complete target-only FASTA with `GN=` annotations; the tiny sequences above only show the format. Native DIA execution uses Nextflow. [Raw sample tables and inputs](setup.md).

### Raw DDA-PASEF

```bash
fun_proteomics setup --non-interactive --project /data/dda-study \
  --mode dda --raw-root /data/raw --samples /data/samples.csv \
  --fasta /data/search_targets_contaminants_decoys.fasta \
  --msfragger-zip /data/tools/MSFragger-4.4.1.zip \
  --ionquant-zip /data/tools/IonQuant-1.11.20.zip \
  --diatracer-zip /data/tools/diatracer-2.2.1.zip \
  --independent-samples --cpus 16 --memory-gib 64 --run
```

**`search_targets_contaminants_decoys.fasta`** — format illustration only; use a complete validated search database:

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

The raw sample table has the same shape shown in the DIA example above. Native DDA execution uses Nextflow. Advanced DDA+, phosphoproteomics and HLA settings use the [advanced reference](advanced_reference.md).

## Setup flags

| Flags | Purpose |
| --- | --- |
| `--project FOLDER`, `--outdir FOLDER` | Output/project folder. These are equivalent in `setup`. |
| `--mode matrix\|dia\|dda\|baf\|continue` | Input and workflow type. |
| `--terminal`, `--non-interactive`, `--run` | Answer questions, supply flags, and optionally launch after saving. |
| `--matrix FILE`, `--samples FILE`, `--sheet NAME` | Matrix, sample table and optional Excel worksheet. |
| `--raw-root FOLDER`, `--fasta FILE` | Raw acquisitions and matching search/annotation database. |
| `--independent-samples` | Confirm independent replicate samples. |
| `--subject-column COLUMN` | Paired donor/subject analysis for matrix or BAF inputs. Do not combine with `--independent-samples`. |
| `--service-matrix FILE`, `--service-samples FILE`, `--service-sheet NAME` | Optional matched service benchmark and mapping. |
| `--quantity auto\|value\|maxlfq\|intensity\|spectral_count` | Matrix quantity-column selection. |
| `--normalization none\|median` | Additional intensity normalization; default `none`. |
| `--no-differential` | Export matrices and QC without statistical comparisons. |
| `--cpus N`, `--memory-gib N` | Requested analysis resources. |
| `--msfragger-zip FILE`, `--ionquant-zip FILE`, `--diatracer-zip FILE`, `--fragpipe-zip FILE` | Registered DDA/BAF archives and optional local FragPipe archive. |
| `--python FILE` | Existing BAF converter Python runtime. |
| `--study-config FILE` | Existing native DDA/DIA JSON configuration. |
| `--source-run FOLDER` | Previous analysis when using `--mode continue`. |
| `--port N`, `--no-browser` | Browser port or print the URL without opening it. Useful over SSH. |

`fun_proteomics setup --help` lists the installed version's options. Guided intensity analyses default to all pairs, BH FDR 0.05 and at least two observed replicates/pairs. Native JSON exposes specialist settings.

## Run saved settings and continue

```bash
fun_proteomics check --config /data/my-study/config/run.json
fun_proteomics run --config /data/my-study/config/run.json
fun_proteomics status --project /data/my-study
fun_proteomics --continue /data/my-study
```

**`config/run.json`** — minimal matrix configuration; normally created by setup:

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

For an existing R installation, global runtime options come before the subcommand:

```bash
fun_proteomics --rscript /usr/bin/Rscript --r-library /data/R-library \
  run --config /data/my-study/config/run.json
```

`--workflow-dir` selects a source checkout; `--cache-dir` selects the launcher cache. `--print-command` previews a command. `python -m fun_proteomics` accepts the same arguments. [Installation](installation.md) · [Continuation and integrity](continuation.md).

For direct native matrix analysis with custom FDR or minimum replicate thresholds, use `fun_proteomics quantify --help`; for the complete Nextflow configuration schema and specialized analyses, see the [advanced reference](advanced_reference.md).
