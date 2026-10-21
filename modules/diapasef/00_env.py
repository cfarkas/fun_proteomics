#!/usr/bin/env python3
"""Create and verify a project-local Python/DIA-NN environment.

python3 00_env.py
python3 00_env.py --root environment --requirements requirements-lock.txt --cache-dir software

A cache is optional: without one the official archive and pinned Python packages
are downloaded. This step does not modify system Python or system installations.
"""

import argparse
from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess
import sys
import urllib.request
import venv
import zipfile


ROOT = Path(__file__).resolve().parent
VERSION = "2.7.0"
ARCHIVE = f"DIA-NN-{VERSION}-Academia-Linux.zip"
URL = f"https://github.com/vdemichev/DiaNN/releases/download/2.0/{ARCHIVE}"
SHA256 = "246ca8baf2e9f308854976f94b8761e452275a97b3632070e275a2b4401cb223"
BINARY_SHA256 = "8108ea214d96043837f2881f0ca405f049d937b50fe6152b608b7015bf6cd49c"
PACKAGES = ["numpy==2.2.6", "pandas==2.2.3", "pyarrow==19.0.1", "matplotlib==3.10.1"]


def checksum(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--requirements", type=Path)
    parser.add_argument("--cache-dir", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    if sys.version_info[:2] != (3, 12):
        raise RuntimeError("Reproduction requires Python 3.12; use python3.12 to run this script")
    software, logs = root / "software", root / "logs"
    software.mkdir(parents=True, exist_ok=True)
    logs.mkdir(exist_ok=True)
    python = root / ".venv/bin/python"
    if not python.exists():
        print("Creating the project Python environment", flush=True)
        venv.EnvBuilder(with_pip=True).create(root / ".venv")
    requirements = args.requirements.read_text() if args.requirements else "\n".join(PACKAGES) + "\n"
    pins = [line.strip() for line in requirements.splitlines() if line.strip() and not line.startswith("#")]
    if any("==" not in pin for pin in pins):
        raise ValueError("Every requirement must have an exact version pin")
    expected = dict(pin.split("==", 1) for pin in pins)
    probe = "import importlib.metadata as m,json; print(json.dumps({d.metadata['Name'].lower():d.version for d in m.distributions()}))"
    actual = json.loads(subprocess.check_output([str(python), "-c", probe], text=True))
    (root / "requirements.txt").write_text(requirements)
    if any(actual.get(name.lower()) != version for name, version in expected.items()):
        command = [str(python), "-m", "pip", "install", "--disable-pip-version-check"]
        wheels = args.cache_dir / "wheels" if args.cache_dir else None
        if wheels and wheels.is_dir() and any(wheels.glob("*.whl")):
            command += ["--no-index", "--find-links", str(wheels.resolve())]
        subprocess.run([*command, "-r", str(root / "requirements.txt")], check=True)
    subprocess.run([str(python), "-m", "pip", "check"], check=True)
    actual = json.loads(subprocess.check_output([str(python), "-c", probe], text=True))
    if any(actual.get(name.lower()) != version for name, version in expected.items()):
        raise RuntimeError("Installed Python versions do not match the requirements lock")

    archive = software / ARCHIVE
    if not archive.is_file() or checksum(archive) != SHA256:
        cached = args.cache_dir / ARCHIVE if args.cache_dir else None
        if cached and cached.is_file():
            if checksum(cached) != SHA256:
                raise RuntimeError("Cached DIA-NN archive failed checksum verification")
            shutil.copy2(cached, archive)
        else:
            print(f"Downloading official DIA-NN {VERSION}", flush=True)
            partial = archive.with_suffix(".zip.part")
            request = urllib.request.Request(URL, headers={"User-Agent": "PASEF-Proteomics/0.3.0"})
            with urllib.request.urlopen(request, timeout=120) as response, partial.open("wb") as output:
                while chunk := response.read(8 * 1024 * 1024):
                    output.write(chunk)
            if checksum(partial) != SHA256:
                raise RuntimeError("DIA-NN archive failed official checksum verification")
            partial.replace(archive)
    engine_root = software / f"diann-{VERSION}"
    candidates = list(engine_root.rglob("diann-linux")) if engine_root.exists() else []
    if not candidates:
        engine_root.mkdir(exist_ok=True)
        with zipfile.ZipFile(archive) as bundle:
            for member in bundle.infolist():
                destination = (engine_root / member.filename).resolve()
                if engine_root.resolve() not in destination.parents:
                    raise RuntimeError(f"Unsafe archive path: {member.filename}")
            bundle.extractall(engine_root)
        candidates = list(engine_root.rglob("diann-linux"))
    if len(candidates) != 1:
        raise RuntimeError(f"Expected one DIA-NN executable, found {candidates}")
    binary = candidates[0]
    if checksum(binary) != BINARY_SHA256:
        raise RuntimeError("DIA-NN executable does not match the pinned official binary")
    binary.chmod(binary.stat().st_mode | 0o111)
    license_file = args.cache_dir / "LICENSE.DIA-NN.txt" if args.cache_dir else None
    if (license_file and license_file.is_file()
            and license_file.resolve() != (software / license_file.name).resolve()):
        shutil.copy2(license_file, software / license_file.name)
    environment = dict(os.environ)
    environment["LD_LIBRARY_PATH"] = str(binary.parent) + ":" + environment.get("LD_LIBRARY_PATH", "")
    result = subprocess.run([str(binary)], cwd=binary.parent, env=environment, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=60)
    (logs / "diann_help.txt").write_text(result.stdout)
    if "DIA-NN" not in result.stdout or "error while loading shared libraries" in result.stdout:
        raise RuntimeError(f"DIA-NN did not start: {result.stdout[:2000]}")
    print("Python and DIA-NN environment ready", flush=True)



if __name__ == "__main__":
    main()
