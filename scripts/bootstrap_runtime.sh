#!/usr/bin/env bash
# Prepare project-local dependencies for a Nextflow launch. Never uses sudo.
set -euo pipefail
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
runtime=$1
requested_r=$2
settings=$3
need_r=${4:-true}
mkdir -p "$runtime"
# Two launchers must not modify the same environment simultaneously.
exec 9>"$runtime/setup.lock"
flock 9
export PROTEOMICS_RUNTIME_DIR="$runtime"
echo '[SETUP] Checking Python, Java and Nextflow runtime…'
bash "$repo/setup_environment.sh"
export PATH="$runtime/env/bin:$PATH"
selected_r=''
if [[ "$need_r" == true ]]; then
    if [[ "$requested_r" != auto ]]; then
        selected_r=$(command -v "$requested_r") || { echo "Rscript not found: $requested_r" >&2; exit 1; }
    elif [[ -x "$runtime/r-env/bin/Rscript" ]]; then
        selected_r="$runtime/r-env/bin/Rscript"
    elif [[ -x /usr/bin/Rscript ]]; then
        selected_r=/usr/bin/Rscript
    elif command -v Rscript >/dev/null 2>&1; then
        selected_r=$(command -v Rscript)
    else
        echo '[SETUP] Installing a local R environment…'
        "$runtime/bin/micromamba" --root-prefix "$runtime/mamba" create -y --prefix "$runtime/r-env" --file "$repo/envs/r-environment.yml" --strict-channel-priority
        selected_r="$runtime/r-env/bin/Rscript"
    fi
    export PATH="$(dirname -- "$selected_r"):$PATH"
    export R_LIBS="$repo/.runtime/R-library${R_LIBS:+:$R_LIBS}"
    echo '[SETUP] Checking R reporting packages…'
    if ! "$selected_r" -e 'p<-c("edgeR","edgeRun","limma","jsonlite","statmod","data.table");quit(status=if(all(vapply(p,requireNamespace,logical(1),quietly=TRUE)) && packageVersion("edgeR") >= "4.4.2" && packageVersion("limma") >= "3.62.2")0 else 1)'; then
        bash "$repo/scripts/install_differential_tools.sh" "$selected_r"
    fi
    "$selected_r" -e 'p<-c("edgeR","edgeRun","limma","jsonlite","statmod","data.table");stopifnot(all(vapply(p,requireNamespace,logical(1),quietly=TRUE)));if(packageVersion("edgeR") < "4.4.2" || packageVersion("limma") < "3.62.2")stop("edgeR >= 4.4.2 and limma >= 3.62.2 are required")'
fi
python - "$settings" "$selected_r" <<'PY'
import json,sys
from pathlib import Path
Path(sys.argv[1]).write_text(json.dumps({'rscript':sys.argv[2]}))
PY
echo '[SETUP] Ready.'
