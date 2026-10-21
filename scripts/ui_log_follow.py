"""Mirror the same on-disk run log used by the browser, including new runs."""
import json
from pathlib import Path
import sys
import time
import urllib.request
from ui_run_state import ACTIVE, read_state

class LogMirror:
    def __init__(self, output=None):
        self.output=output or sys.stdout
        self.path=None
        self.offset=0
    def update(self, state):
        path=state.get('log_path')
        if not path:return
        if path!=self.path:self.path=path;self.offset=0
        target=Path(path)
        if not target.is_file():return
        if target.stat().st_size<self.offset:self.offset=0
        with target.open('rb') as handle:
            handle.seek(self.offset)
            chunk=handle.read()
            self.offset=handle.tell()
        if chunk:self.output.write(chunk.decode('utf-8',errors='replace'));self.output.flush()

def follow(url, retry_seconds=10):
    mirror=LogMirror()
    state={};disconnected=None
    print('Attached to the existing interface. This terminal follows the browser run log.',flush=True)
    print('Ctrl+C detaches this terminal only. Use Stop run in the browser to stop analysis.',flush=True)
    try:
        while True:
            try:
                with urllib.request.urlopen(url+'/api/status',timeout=5) as response:state=json.load(response)
                if disconnected is not None:
                    print('[INTERFACE RECONNECTED] Following the current run log.',flush=True)
                disconnected=None
                mirror.update(state)
            except (OSError,ValueError) as error:
                if disconnected is None:
                    disconnected=time.monotonic()
                    print('[INTERFACE DISCONNECTED] '+str(error)+'. Retrying; analysis completion is not yet verified.',flush=True)
                # Keep following the saved run even if its web server is gone.
                if state.get('state_path'):
                    try:
                        state=read_state(state['state_path'])
                        mirror.update(state)
                    except (OSError,ValueError,KeyError) as saved_error:
                        if time.monotonic()-disconnected>=retry_seconds:
                            print('Cannot verify saved run status: '+str(saved_error),flush=True)
                            return 1
                    else:
                        if state['state'] not in ACTIVE:
                            print(f"[ANALYSIS STATUS] {state['state']} | exit {state.get('exit_code')} | log: {state['log_path']}",flush=True)
                            return 0 if state['state']=='completed' and state.get('exit_code')==0 else 1
                else:
                    mirror.update(state)
                    if time.monotonic()-disconnected>=retry_seconds:
                        print('Interface unavailable. Analysis completion could not be verified. Reopen the interface and inspect the output logs.',flush=True)
                        return 1
            time.sleep(1)
    except KeyboardInterrupt:
        print('Detached. The existing interface and its analysis were not stopped.',flush=True)
        return 0


def mirror_runner(runner, stopped):
    """Terminal output is optional; a closed terminal must not affect a run."""
    mirror=LogMirror()
    try:
        while not stopped.is_set():
            mirror.update(runner.status())
            stopped.wait(0.5)
        mirror.update(runner.status())
    except (OSError,ValueError):
        # The supervisor still writes the complete log and final status to disk.
        return
