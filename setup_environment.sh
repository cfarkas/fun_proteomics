#!/usr/bin/env bash
set -euo pipefail
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
runtime=${PROTEOMICS_RUNTIME_DIR:-"$repo/.runtime"}
[[ $(uname -s) == Linux && $(uname -m) == x86_64 ]] || { echo 'Linux x86_64 is required.' >&2; exit 1; }
mkdir -p "$runtime/bin"
manager="$runtime/bin/micromamba"
if [[ ! -x "$manager" ]]; then
    curl --fail --location --retry 3 'https://github.com/mamba-org/micromamba-releases/releases/download/2.3.3-0/micromamba-linux-64' -o "$manager.download"
    echo "9496f94a8b78c536573c93d946ec9bba74bd9ff79ee55aaa4b546e30db8f511b  $manager.download" | sha256sum --check
    mv "$manager.download" "$manager"
    chmod +x "$manager"
fi
spec="$repo/envs/environment.yml"
[[ ! -f "$repo/envs/linux-64.lock.txt" ]] || spec="$repo/envs/linux-64.lock.txt"
spec_hash=$(sha256sum "$spec" | cut -d ' ' -f 1)
if [[ -f "$runtime/environment.sha256" ]]; then
    [[ $(cat "$runtime/environment.sha256") == "$spec_hash" ]] || { echo 'Environment specification changed; choose a new PROTEOMICS_RUNTIME_DIR and rerun setup.' >&2; exit 1; }
else
    "$manager" --root-prefix "$runtime/mamba" create -y --prefix "$runtime/env" --file "$spec" --strict-channel-priority
    "$manager" --root-prefix "$runtime/mamba" list --prefix "$runtime/env" --explicit > "$runtime/installed-linux-64.lock.txt"
    printf '%s\n' "$spec_hash" > "$runtime/environment.sha256"
fi
"$runtime/env/bin/python" -c 'import openpyxl, jsonschema; print("Python reporting dependencies: OK")'
"$runtime/env/bin/java" -version
echo "Environment ready: $runtime/env"
