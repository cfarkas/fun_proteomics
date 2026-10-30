"""Launch and reuse a session-bound setup server independently of the terminal."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import threading
from urllib.error import URLError
from urllib.request import Request, urlopen
import webbrowser

from . import __version__
from .continuation import read_json, write_json


def reachable(record):
    from .cli import process_identity
    if not record.get('identity') or process_identity(record.get('pid')) != record['identity']:
        return False
    try:
        request = Request(record['url'].split('/#')[0] + '/api/initial',
                          headers={'X-Session-Token': record['token']})
        with urlopen(request, timeout=2) as response:
            return response.status == 200
    except (KeyError, ValueError, OSError, URLError):
        return False


def start(args, workflow, timeout=30):
    from .setup_web import worker_environment
    settings = vars(args).copy()
    settings['no_browser'] = True
    request = {'args': settings, 'workflow': str(Path(workflow).resolve()),
               'cwd': str(Path.cwd()), 'version': __version__,
               'desktop': {name: os.environ.get(name, '') for name in ('DISPLAY', 'WAYLAND_DISPLAY')}}
    key = hashlib.sha256(json.dumps(request, sort_keys=True).encode()).hexdigest()[:20]
    directory = Path(workflow) / '.runtime' / 'setup-sessions'
    directory.mkdir(parents=True, exist_ok=True, mode=0o700)
    state = directory / (key + '.json')
    logfile = directory / (key + '.log')
    with (directory / (key + '.lock')).open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if state.is_file():
            record = read_json(state)
            if reachable(record):
                return record
        request_path = directory / (key + '.request.json')
        write_json(request_path, request)
        request_path.chmod(0o600)
        command = [sys.executable, '-u', '-m', 'fun_proteomics.setup_session', str(request_path), str(state)]
        with logfile.open('a') as log:
            logfile.chmod(0o600)
            process = subprocess.Popen(command, cwd=request['cwd'], env=worker_environment(workflow, args),
                                       stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT,
                                       start_new_session=True)
        try:
            deadline = time.monotonic() + timeout
            while time.monotonic() < deadline:
                if process.poll() is not None:
                    raise ValueError(f'Setup server exited with code {process.returncode}. See {logfile}')
                if state.is_file():
                    record = read_json(state)
                    if record.get('pid') == process.pid and reachable(record):
                        threading.Thread(target=process.wait, daemon=True).start()
                        return record
                time.sleep(.1)
            raise ValueError(f'Setup server did not start within {timeout}s. See {logfile}')
        except BaseException:
            process.terminate()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
            raise


def launch(args, workflow):
    record = start(args, workflow)
    print('fun_proteomics | Guided setup\n' + record['url'], flush=True)
    print('Setup stays available after this terminal closes. Run the same command to reopen it.', flush=True)
    webbrowser.open(record['url'])
    return 0


def main():
    from .cli import process_identity
    from .setup_web import serve_foreground
    request_path, state_path = map(Path, sys.argv[1:])
    request = read_json(request_path)
    args = argparse.Namespace(**request['args'])
    os.chdir(request['cwd'])
    def ready(url, app, server):
        write_json(state_path, {'pid': os.getpid(), 'identity': process_identity(os.getpid()),
                               'url': url, 'token': app.token, 'version': __version__})
        state_path.chmod(0o600)
    return serve_foreground(args, Path(request['workflow']), ready=ready)


if __name__ == '__main__':
    raise SystemExit(main())
