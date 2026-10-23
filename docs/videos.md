# Interactive walkthrough videos

[Open the video guide](https://cfarkas.github.io/fun_proteomics/). The five captioned recordings show actual application sessions. The quick start analyzes the published dataset; the other four use explicitly synthetic examples.

| Video | What it shows | Transcript |
| --- | --- | --- |
| [Real-data quick start](https://github.com/cfarkas/fun_proteomics/releases/download/v0.9.1/ubiquitin-quickstart.mp4) | Choose files and an output folder, assign four conditions, run 12 protein/gene comparisons, open CSV/TSV matrices, and continue verified results. | [Captions](ubiquitin-quickstart.vtt) |
| [Browser: matrix analysis](https://github.com/cfarkas/fun_proteomics/releases/download/v0.7.0/browser-analysis.mp4) | Discover matrix samples, drag them into conditions, confirm an independent design, save, run real all-pairs statistics, and open CSV/TSV matrices. | [Captions](browser-analysis.vtt) |
| [Terminal: launch and continue](https://github.com/cfarkas/fun_proteomics/releases/download/v0.7.0/terminal-analysis.mp4) | Answer terminal questions, run the synthetic analysis, and verify completed outputs with `--continue`. | [Captions](terminal-analysis.vtt) · [session transcript](terminal-analysis.txt) |
| [Raw DIA: interactive setup](https://github.com/cfarkas/fun_proteomics/releases/download/v0.9.0/raw-dia-setup.mp4) | Scan `.d` folders, drag acquisition cards into conditions, choose FASTA, review all pairs, and save a validated raw project. This fixture has synthetic acquisition metadata; no spectrum search is performed. | [Captions](raw-dia-setup.vtt) |
| [Raw BAF: paired setup](https://github.com/cfarkas/fun_proteomics/releases/download/v0.9.0/raw-baf-setup.mp4) | Assign acquisitions, combine technical injections under biological sample IDs, enter donors, review six paired comparisons, and save. Synthetic placeholders are never searched. | [Captions](raw-baf-setup.vtt) |

To follow along, use the [quick start](quick_start.md). The browser and terminal videos execute the same analysis code available after `pip install fun_proteomics`. The terminal video replays an actual captured terminal session at a readable pace. The raw videos demonstrate configuration only because their fixtures contain no real spectra.

The public-data recording uses `scripts/record_public_quickstart.py`; the four synthetic recordings use `scripts/record_walkthroughs.py` after installing the package and the recording-only Playwright/ffmpeg tools. Media files are release assets, keeping them outside the Python wheel. Source, version, commands and media checksums are recorded with the release. The BAF recorder can validate user-supplied registered archives and reader paths; those dependencies are not included in the recordings. Its command-line help lists recording options.
