# Benchmark any run against service quantification

Service benchmarking is available for matrix imports, raw DDA-PASEF, DIA-PASEF and raw BAF projects. Select the service matrix during setup. For an existing run, use the independent benchmark command; previous search results remain available.

```bash
fun_proteomics benchmark \
  --matrix /data/search/combined_protein.tsv \
  --samples /data/samples.csv \
  --service-matrix /data/service.xlsx \
  --service-sheet 'Identified proteins' \
  --service-samples /data/service_samples.csv \
  --outdir /data/service-benchmark
```

The pipeline matrix may instead be a DIA-NN matrix or an exported protein-intensity CSV/TSV. Omit `--service-sheet` for automatic Excel detection or text inputs. Omit `--service-samples` when provider columns already match sample IDs. `service_samples.csv` contains `sample_id,service_sample`, with exact provider aliases. `samples.csv` contains unique pipeline `sample_id` values; include condition and replicate for differential analysis. To also calculate every condition contrast, add `--differential`; for paired data, add `--subject-column donor` and the donor column. Technical injections must already be aggregated when importing matrices.

Open `OUTPUT/final_results/START_HERE.html`. All matrices and benchmark tables have CSV and TSV versions.

| Table | Meaning |
| --- | --- |
| `reported_detection_summary` | All positive groups as reported by each analysis, including incomplete groups; no shared/gained/lost identity claim is made for incomplete groups. |
| `detection_summary` | Exact complete protein groups: pipeline detected, service detected, shared, pipeline-only, service-only, net difference and percentage change, per matched sample and in the union. |
| `protein_membership` | The actual shared, additional and missing group identifiers. |
| `group_completeness_audit` | All reported/positive groups on each side, and the incomplete groups excluded from exact matching. |
| `accession_overlap_summary`, `accession_membership` | Supporting overlap of accessions from complete groups. These are group members, not independently identified proteins. |
| `sample_mapping` | Exact source columns used in the service comparison. |
| `service_quantification` | The matched service quantities used to count detection. |

Positive abundance counts as detected; blank/NA and zero do not. No imputation or significance threshold is used to inflate detections. Decoys and contaminants are excluded. The union counts a group once across matched samples. Percentage change is undefined when the service has zero detections. Clipped provider FASTA headers cannot establish complete protein-group membership and are excluded from exact matching, with their number explicitly reported.

Use matching measurements and samples. Comparing spectral counts to LFQ positivity, different databases or different FDR thresholds answers different questions. Match-between-runs can quantify a protein without direct MS2 evidence in that sample. A positive net difference is additional reported coverage, not proof of higher accuracy. Separate search modes and cohorts; never count a DDA and DDA+ search of one sample as biological replicates.

Native DDA JSON accepts `service_quantification` keyed by cohort, with `matrix`, optional `sheet`, and optional `samples` containing exact service aliases. Native DIA retains its `service_matrices` configuration. Guided setup writes these settings automatically.
