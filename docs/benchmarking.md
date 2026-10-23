# Compare a run with service quantification

Use this guide for matrix, DDA-PASEF, DIA-PASEF or [raw BAF](raw_baf.md) results. Benchmarking compares positive quantities over the same biological samples. It reports shared groups, extra groups and missing groups; it does not use a significance cutoff.

## Prepare the files

The numbers below illustrate file formats. Use your own measured quantities.

**`proteins.csv`** — your pipeline protein quantities:

```csv
Protein.Group,Genes,C1,C2,T1,T2
P04637,TP53,120000,135000,260000,245000
P60709,ACTB,480000,505000,490000,510000
P02768,ALB,90000,85000,,100000
```

**`samples.csv`** — the biological samples to match:

```csv
sample_id,condition,replicate
C1,Control,1
C2,Control,2
T1,Treatment,1
T2,Treatment,2
```

**`service.xlsx`** — CSV-style preview of the Excel worksheet named Identified proteins:

```csv
Protein IDs,Gene names,LFQ intensity Lab_C1,LFQ intensity Lab_C2,LFQ intensity Lab_T1,LFQ intensity Lab_T2
P04637,TP53,118000,132000,255000,240000
P60709,ACTB,475000,500000,485000,505000
```

**`--service-sheet` is a worksheet name inside the workbook**, not another file. A service CSV/TSV has the same table structure and needs no sheet option. Automatic worksheet detection is used when the sheet name is omitted.

Use the following mapping only when the service names differ. Here `C1` matches `LFQ intensity Lab_C1`; the measurement prefix is recognized automatically.

**`service_samples.csv`** — exact pipeline-to-service sample aliases:

```csv
sample_id,service_sample
C1,Lab_C1
C2,Lab_C2
T1,Lab_T1
T2,Lab_T2
```

Omit the mapping file when both tables already use the same sample names. Each table must contain one quantity per biological sample; combine technical injections before import.

## Run the comparison

```bash
fun_proteomics benchmark \
  --matrix /data/proteins.csv --samples /data/samples.csv \
  --service-matrix /data/service.xlsx --service-sheet 'Identified proteins' \
  --service-samples /data/service_samples.csv \
  --outdir /data/service-benchmark
```

For a new analysis, choose the same service inputs in browser setup. To calculate differential abundance as well, add `--differential`; for paired data, use `--subject-column donor` and the [paired sample table](provider_workbooks.md#repeated-donors-and-technical-injections).

## Read the results

Open `service-benchmark/final_results/START_HERE.html`. All tables have CSV and TSV versions. This is an illustrative excerpt of `detection_summary.csv`:

```csv
sample_id,pipeline_detected,service_detected,shared,pipeline_only,service_only,net_difference,percent_change
C1,100,90,85,15,5,10,11.111111
C2,95,100,90,5,10,-5,-5
```

| Table | Meaning |
| --- | --- |
| `reported_detection_summary` | All positive groups reported by each analysis, including incomplete groups. |
| `detection_summary`, `protein_membership` | Exact complete-group counts and shared/additional/missing identifiers. |
| `group_completeness_audit` | Incomplete groups excluded from exact identity matching. |
| `accession_overlap_summary`, `accession_membership` | Supporting accession overlap; group members are not independent protein detections. |
| `sample_mapping`, `service_quantification` | Matched source columns and the service quantities used. |

Positive quantities count as detections; zero, blank and NA do not. Decoys and contaminants are excluded. The union counts a group once across matched samples. Percentage change is undefined if the service has zero detections. Incomplete provider headers stay in reported counts but cannot establish exact group identity.

Compare compatible measurements. Databases, grouping, FDR and match-between-runs can change coverage; more detections alone do not establish higher accuracy. Keep cohorts and search modes separate.

<details>
<summary>Native configuration preview</summary>

Guided setup writes these settings automatically. In a native DDA JSON, `service_quantification` is keyed by cohort; this is an excerpt:

```json
{
  "service_quantification": {
    "study": {"matrix": "/data/service.xlsx", "sheet": "Identified proteins", "samples": [{"sample_id": "C1", "service_sample": "Lab_C1"}]}
  }
}
```

Native DIA uses `service_matrices`, keyed by sample type, for example:

```json
{"service_matrices": {"serum": "/data/service.tsv"}}
```

The service TSV has the same quantities as the workbook preview above, separated by tabs. Use the complete [native configuration examples](advanced_reference.md) when launching a raw workflow.

</details>
