"""Discover pinned archives and install tools locally without bypassing registration."""
import hashlib,json,os,shutil,tempfile,threading,urllib.request,zipfile
from pathlib import Path
from functools import lru_cache
ROOT=Path(__file__).resolve().parents[1]
SPECS=json.loads((ROOT/'assets/tool_versions.json').read_text())
KEYS={'FragPipe':'fragpipe_zip','MSFragger':'msfragger_zip','IonQuant':'ionquant_zip','diaTracer':'diatracer_zip'}
ENTRIES={'FragPipe':'fragpipe-24.0/bin/fragpipe','MSFragger':'MSFragger-4.4.1/MSFragger-4.4.1.jar','IonQuant':'IonQuant-1.11.20/IonQuant-1.11.20.jar','diaTracer':'diatracer-2.2.1/diaTracer-2.2.1.jar'}
LOCK=threading.RLock()
@lru_cache(maxsize=128)
def _hash(path,size,mtime):
 with open(path,'rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def verified(path,name):
 try:
  p=Path(path);s=p.stat()
  return p.is_file() and _hash(str(p.resolve()),s.st_size,s.st_mtime_ns)==SPECS[name]['archive_sha256']
 except (OSError,TypeError,ValueError):return False
def roots(extra=''):
 values=[ROOT/'.runtime',ROOT.parent,Path.home(),Path.home()/'Downloads',Path.home()/'.local/share','/opt','/usr/local']
 if extra:values.insert(0,Path(extra).expanduser())
 return list(dict.fromkeys(Path(p).resolve() for p in values if Path(p).is_dir()))
def discover(tools=None,search_root=''):
 tools=tools or {};candidates={n:[] for n in KEYS};binaries={n:[] for n in KEYS};scanned=[]
 for n,k in KEYS.items():
  if tools.get(k):candidates[n].append(Path(tools[k]).expanduser())
 seen=set();budget=12000
 for root in roots(search_root):
  scanned.append(str(root))
  for base,dirs,files in os.walk(root,followlinks=False):
   path=Path(base)
   dirs[:]=[d for d in dirs if not d.startswith('.') and d not in {'node_modules','mamba','env','r-env','raw','RawDataDDA','RawDataDIA'} and not d.endswith('.d') and not d.startswith(('work','results'))]
   if len(path.relative_to(root).parts)>=4:dirs[:]=[]
   if str(path) in seen:continue
   seen.add(str(path));budget-=1
   for f in files:
    for n in KEYS:
     if n.lower() in f.lower():
      if f.lower().endswith('.zip'):candidates[n].append(path/f)
      elif f.endswith('.jar') or f=='fragpipe':binaries[n].append(str(path/f))
   if budget<=0:break
  if budget<=0:break
 result={};found={}
 for n,k in KEYS.items():
  archive=next((str(p.resolve()) for p in dict.fromkeys(candidates[n]) if verified(p,n)),None)
  if archive:found[k]=archive
  installed=ROOT/'.runtime/tools'/n
  entry=installed/ENTRIES[n]
  result[n]={'key':k,'version':SPECS[n]['version'],'archive':archive,'installed':entry.is_file(),'binaries':binaries[n][:8], 'downloadable':'url' in SPECS[n], 'registration':SPECS[n].get('registration'), 'state':'installed' if entry.is_file() and archive else 'archive verified' if archive else 'installation found; verified ZIP required' if binaries[n] else 'missing'}
 return {'tools':found,'status':result,'searched':scanned,'limited':budget<=0}
def ensure_archive(name,tools=None,force=False):
 if name not in KEYS:raise ValueError('Unknown tool')
 tools=tools or {};spec=SPECS[name];path=tools.get(KEYS[name])
 if not force or 'url' not in spec:
  if path and verified(path,name):return Path(path).resolve()
  found=discover(tools)['tools'].get(KEYS[name])
  if found:return Path(found)
 if 'url' not in spec:raise ValueError(f'{name} {spec["version"]}: register at {spec["registration"]}, download the ZIP, then Detect tools again. Force reinstall requires the registered ZIP too.')
 cache=ROOT/'.runtime/downloads';cache.mkdir(parents=True,exist_ok=True)
 destination=cache/Path(spec['url']).name
 with tempfile.NamedTemporaryFile(dir=cache,suffix='.download',delete=False) as f:tmp=Path(f.name)
 try:
  with urllib.request.urlopen(spec['url'],timeout=90) as response,tmp.open('wb') as out:shutil.copyfileobj(response,out)
  if not verified(tmp,name):raise ValueError(f'{name} download checksum mismatch; no installation changed.')
  os.replace(tmp,destination)
 finally:tmp.unlink(missing_ok=True)
 return destination

def install(name,tools=None,force=False):
 with LOCK:
  archive=ensure_archive(name,tools,force)
  folder=ROOT/'.runtime/tools';folder.mkdir(parents=True,exist_ok=True);target=folder/name
  if force or not (target/ENTRIES[name]).is_file():
   with tempfile.TemporaryDirectory(dir=folder,prefix='install-') as temp:
    stage=Path(temp)/'new';stage.mkdir()
    with zipfile.ZipFile(archive) as z:
     for info in z.infolist():
      if not (stage/info.filename).resolve().is_relative_to(stage.resolve()) or ((info.external_attr>>16)&0o170000)==0o120000:raise ValueError('Unsafe tool archive entry')
     z.extractall(stage)
     for info in z.infolist():
      path=stage/info.filename
      if path.is_file():path.chmod(0o755 if (info.external_attr>>16)&0o111 else 0o644)
    if not (stage/ENTRIES[name]).is_file():raise ValueError('Unexpected tool archive layout')
    if name=='FragPipe':(stage/ENTRIES[name]).chmod(0o755)
    previous=Path(temp)/'previous'
    if target.exists():target.rename(previous)
    try:stage.rename(target)
    except Exception:
     if previous.exists():previous.rename(target)
     raise
  return {'name':name,'path':str(target),'tools':{KEYS[name]:str(archive)},'installed':True}
