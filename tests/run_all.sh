#!/usr/bin/env bash
# Select the same installed runtime used by the workflow.
set -euo pipefail
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
cd "$repo"
python_bin=${PROTEOMICS_TEST_PYTHON:-"$repo/.runtime/env/bin/python"}
if [[ ! -x "$python_bin" ]]; then python_bin=$(command -v python3); fi
if [[ -z ${PROTEOMICS_RSCRIPT:-} ]]; then
  if [[ -x "$repo/.runtime/r-env/bin/Rscript" ]]; then
    export PROTEOMICS_RSCRIPT="$repo/.runtime/r-env/bin/Rscript"
  elif [[ -x /usr/bin/Rscript ]]; then
    export PROTEOMICS_RSCRIPT=/usr/bin/Rscript
  else
    export PROTEOMICS_RSCRIPT=$(command -v Rscript)
  fi
fi
export R_LIBS="$repo/.runtime/R-library${R_LIBS:+:$R_LIBS}"
"$python_bin" scripts/validate_repo.py
"$python_bin" -m unittest discover -s tests -v
