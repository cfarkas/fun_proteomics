#!/usr/bin/env bash
set -euo pipefail
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
rscript="${1:-Rscript}"
lib="$repo/.runtime/R-library"
mkdir -p "$lib"
export R_LIBS="$lib${R_LIBS:+:$R_LIBS}"
"$rscript" -e 'lib<-strsplit(Sys.getenv("R_LIBS"),":",fixed=TRUE)[[1]][1];options(repos=c(CRAN="https://cloud.r-project.org"));if(!requireNamespace("BiocManager",quietly=TRUE))install.packages("BiocManager",lib=lib);missing<-c("edgeR","limma","jsonlite","statmod","data.table","fdrtool");missing<-missing[!vapply(missing,requireNamespace,logical(1),quietly=TRUE)];if(length(missing))BiocManager::install(missing,lib=lib,ask=FALSE,update=FALSE)'
# Early edgeR 4.0 small-count code can fail on sparse proteomics matrices.
# These matching Bioconductor releases support R >= 3.6 and are installed locally.
for spec in 'limma 3.62.2 6abb60905719281fa3ccd5f2a6c201e2878a4c0a0a8e30cec662a541e67fb475' 'edgeR 4.4.2 8620592ce1f3ca573c8a59852a50f6b50c2c942891f92c1142f002eb23607e3f'; do
    read -r package version checksum <<< "$spec"
    if ! "$rscript" -e 'a<-commandArgs(TRUE);quit(status=if(requireNamespace(a[1],quietly=TRUE)&&packageVersion(a[1])>=a[2])0 else 1)' "$package" "$version"; then
        source_archive="$repo/.runtime/${package}_${version}.tar.gz"
        curl -L --fail -o "$source_archive" "https://bioconductor.org/packages/3.20/bioc/src/contrib/${package}_${version}.tar.gz"
        python - "$source_archive" "$checksum" <<'VERIFY'
import hashlib,sys
from pathlib import Path
if hashlib.sha256(Path(sys.argv[1]).read_bytes()).hexdigest()!=sys.argv[2]:
    raise SystemExit('R package archive checksum mismatch')
VERIFY
        "$rscript" -e 'a<-commandArgs(TRUE);install.packages(a[1],repos=NULL,type="source",lib=a[2])' "$source_archive" "$lib"
    fi
done
"$rscript" -e 'stopifnot(packageVersion("edgeR") >= "4.4.2",packageVersion("limma") >= "3.62.2")'
archive="$repo/.runtime/edgeRun_1.0.9.tar.gz"
curl -L --fail -o "$archive" https://cran.r-project.org/src/contrib/Archive/edgeRun/edgeRun_1.0.9.tar.gz
python - "$archive" <<'CHECK'
import hashlib,sys
from pathlib import Path
assert hashlib.sha256(Path(sys.argv[1]).read_bytes()).hexdigest()=='21b6f144a6d1598bd9ad810fdddc233252fc83fc9a3700e26ac154ca1d62b634','edgeRun archive checksum mismatch'
CHECK
"$rscript" -e 'a<-commandArgs(TRUE);install.packages(a[1],repos=NULL,type="source",lib=a[2])' "$archive" "$lib"
