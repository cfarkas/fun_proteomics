# Quick start: a complete small analysis

[Install once](installation.md). This example uses six invented samples, three conditions and ten protein groups. It runs the real matrix and statistical reporting code; the measurements have no biological interpretation.

```bash
fun_proteomics demo --project "$PWD/proteomics-demo"
fun_proteomics setup --project "$PWD/proteomics-demo"
```

The first command creates `proteomics-demo-inputs/` beside the project. It does not launch an analysis. In the browser:

1. Choose **Quantification matrix**. Browse to `proteomics-demo-inputs/proteins.csv`, then click **Discover samples**. Drag C1/C2 into Control, A1/A2 into Treatment A and B1/B2 into Treatment B (rename/add condition columns). Alternatively, select `samples.csv` before discovery to prefill these assignments.
2. Choose **Review analysis**. Keep all-versus-all comparisons. Confirm independent biological replicates for this synthetic example. Select `service.csv` for the optional detection benchmark.
3. Choose **Save configuration and check**, then **Run / continue analysis**.
4. Open **Matrices · CSV & TSV**. There are three condition pairs at protein and gene levels, six differential result tables in total.

For the same interactive experience in your terminal:

```bash
fun_proteomics setup --terminal --project "$PWD/proteomics-demo" --run
```

For a reproducible scripted example:

```bash
fun_proteomics setup --non-interactive --project "$PWD/proteomics-demo" \
  --mode matrix --matrix "$PWD/proteomics-demo-inputs/proteins.csv" \
  --samples "$PWD/proteomics-demo-inputs/samples.csv" \
  --service-matrix "$PWD/proteomics-demo-inputs/service.csv" \
  --independent-samples --run
```

Use a new project folder for different inputs. Repeat the same saved project to continue. Try the integrity check:

```bash
fun_proteomics --continue "$PWD/proteomics-demo"
```

A completed unchanged delivery prints `Completed outputs verified. No analysis needs repeating.`

[Watch the browser and terminal walkthroughs](https://cfarkas.github.io/fun_proteomics/) · [Understand the results](outputs.md).
