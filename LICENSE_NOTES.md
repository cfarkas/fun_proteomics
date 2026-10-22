# Licensing and redistribution

The original workflow code and documentation are distributed under the MIT license in LICENSE. Copyright (c) 2026 cfarkas and fun_proteomics contributors. Third-party tools retain their own licenses; the MIT license does not relicense them.

FragPipe, MSFragger, IonQuant, diaTracer, DIA-NN, Nextflow and the Python/R dependencies retain their own licenses and citation requirements. Obtain registered software through the official distribution links in assets/tool_versions.json. Their binaries, vendor SDKs, private URLs, registration details and keys are not included.

This source package contains no experimental acquisitions or sample-level study data. Aggregate execution and benchmark results are documented separately. The FragPipe settings file is a configuration, not a bundled copy of the application.

Native BAF conversion optionally installs [pyBaf2Sql](https://github.com/gtluu/pyBaf2Sql) (Apache-2.0 wrapper, separately supplied Bruker Baf2Sql 2.9.0 library and third-party notices) and [psims](https://github.com/mobiusklein/psims) (Apache-2.0). These dependencies are installed in a separate user runtime, not redistributed in the MIT fun_proteomics wheel or source archive. Their original licenses remain applicable.
