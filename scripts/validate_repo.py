#!/usr/bin/env python3
"""Check the explicit publication source set without traversing private data."""
import ast
import json
import re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ROOT_FILES={'README.md','GITHUB_UPLOAD.md','OUTPUTS.md','LICENSE_NOTES.md','SOURCE_SHA256SUMS','main.nf','nextflow.config','setup_environment.sh','.gitignore','.gitattributes','pyproject.toml','MANIFEST.in','PACKAGE_README.md','LICENSE'}
SOURCE_DIRS={'assets','pipeline_lib','scripts','tests','envs','examples','.github','modules','src','docs'}
EXTENSIONS={'.html','.py','.R','.sh','.json','.txt','.md','.nf','.config','.workflow','.yml','.yaml','.tsv','.csv','.toml','.vtt','.cast'}

def source_files(directory):
    for path in sorted(directory.iterdir()):
        if path.is_symlink(): continue
        if directory==ROOT and path.name not in ROOT_FILES|SOURCE_DIRS: continue
        if directory == ROOT/'src' and path.name != 'fun_proteomics': continue
        if path.name == '_workflow': continue
        if directory == ROOT/'modules' and path.name != 'diapasef': continue
        if path.name in {'__pycache__','.git'} or path.name.endswith(('.pyc','.local.json','.local.config')): continue
        if path.is_dir(): yield from source_files(path)
        else: yield path

def main():
    files=list(source_files(ROOT)); failures=[]
    for path in files:
        if path.name not in ROOT_FILES and path.suffix not in EXTENSIONS:
            failures.append(f'Unexpected source file: {path.relative_to(ROOT)}');continue
        text=path.read_text()
        if path.suffix=='.py': ast.parse(text,filename=str(path))
        if path.suffix=='.json': json.loads(text)
        prohibited=['/media/'+'server/','/home/'+'server/','cfarkas'+'@ucsc.cl','carlosfarkas'+'@gmail.com','download.php?'+'token=', 'TH'+'P1','th'+'p1']
        if any(token in text for token in prohibited): failures.append(f'Private or study-specific value: {path.relative_to(ROOT)}')
        if re.search(r'-----BEGIN (?:RSA |OPENSSH )?PRIVATE KEY-----',text): failures.append(f'Credential: {path.relative_to(ROOT)}')
        if path.suffix=='.md':
            for href in re.findall(r'\]\(([^)]+)\)',text):
                if '://' in href or href.startswith('#'):continue
                dest=href.split('#')[0]
                if dest and not (path.parent/dest).exists():failures.append(f'Broken link: {path.relative_to(ROOT)}: {href}')
    if failures:raise SystemExit('\n'.join(failures))
    print(f'PASS: {len(files)} publication files; syntax, JSON, links and privacy checks')
if __name__=='__main__':main()
