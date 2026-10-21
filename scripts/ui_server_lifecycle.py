"""Keep the local setup server independent of its launcher terminal."""
import fcntl
import json
from pathlib import Path
import subprocess
import sys
import threading
import time
import urllib.error
import urllib.request


def server_info(url, root):
    try:
        with urllib.request.urlopen(url + '/api/config', timeout=1) as response:
            info = json.load(response)
    except urllib.error.HTTPError as error:
        raise ValueError(f'Cannot use {url}: HTTP {error.code}. Choose another --ui_port.') from error
    except (OSError, urllib.error.URLError):
        return None
    if not isinstance(info.get('root'), str) or not info['root'] or Path(info['root']).resolve() != Path(root).resolve():
        raise ValueError(f'{url} belongs to another project. Choose another --ui_port.')
    return info


def start_server(root, config, port, overrides=None, startup_timeout=20):
    """Start once, verify readiness, and return (URL, log path, new process or None)."""
    root = Path(root).resolve()
    url = f'http://127.0.0.1:{port}'
    state_dir = root / '.runtime'
    state_dir.mkdir(parents=True, exist_ok=True)
    log_path = state_dir / f'ui_server_{port}.log'
    with (state_dir / f'ui_server_{port}.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if server_info(url, root) is not None:
            return url, log_path, None
        command = [sys.executable, '-u', str(root / 'scripts/config_server.py'),
                   '--config', str(Path(config).resolve()), '--port', str(port),
                   '--overrides', json.dumps(overrides or {}), '--no-log-mirror']
        with log_path.open('a') as log:
            process = subprocess.Popen(command, cwd=root, stdin=subprocess.DEVNULL,
                                       stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        try:
            deadline = time.monotonic() + startup_timeout
            while time.monotonic() < deadline:
                if process.poll() is not None:
                    raise RuntimeError(f'Setup server exited with code {process.returncode}. See {log_path}.')
                if server_info(url, root) is not None:
                    threading.Thread(target=process.wait, daemon=True).start()
                    return url, log_path, process
                time.sleep(0.1)
            raise RuntimeError(f'Setup server did not respond within {startup_timeout}s. See {log_path}.')
        except Exception:
            process.terminate()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
            raise
