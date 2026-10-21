"""Durable run status and a supervisor independent of the interface process."""
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time

ACTIVE = ('running', 'stopping')


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode='w', dir=path.parent, delete=False) as handle:
        temporary = Path(handle.name)
        json.dump(data, handle)
        handle.write('\n')
    try:
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def process_identity(pid):
    """Reject dead processes and PID reuse, including after a reboot."""
    if not pid:
        return None
    try:
        fields = Path(f'/proc/{int(pid)}/stat').read_text().rsplit(')', 1)[1].split()
        if fields[0] == 'Z':
            return None
        boot = Path('/proc/sys/kernel/random/boot_id').read_text().strip()
        return boot + ':' + fields[19]
    except (OSError, ValueError, IndexError):
        return None


def process_matches(pid, identity):
    return identity is not None and process_identity(pid) == identity


def read_state(path):
    data = json.loads(Path(path).read_text())
    if data['state'] in ACTIVE:
        supervisor = process_matches(data.get('supervisor_pid'), data.get('supervisor_identity'))
        child = process_matches(data.get('pid'), data.get('pid_identity'))
        starting = not data.get('supervisor_pid') and time.time() - data['started'] < 10
        if not (supervisor or child or starting):
            # Losing a controller is never evidence of successful completion.
            data.update(state='interrupted', exit_code=None,
                        message='Run processes disappeared before an exit status was saved. Inspect the run log before resuming.')
        elif Path(data['stop_path']).exists():
            data['state'] = 'stopping'
    return data


def supervise(path):
    path = Path(path)
    data = json.loads(path.read_text())
    data.update(supervisor_pid=os.getpid(), supervisor_identity=process_identity(os.getpid()))
    write_json(path, data)
    stopping = False

    def request_stop(signum, frame):
        nonlocal stopping
        stopping = True

    signal.signal(signal.SIGTERM, request_stop)
    signal.signal(signal.SIGINT, request_stop)
    with Path(data['log_path']).open('a') as log:
        try:
            # Nextflow writes directly to disk, never to a pipe owned by the UI.
            proc = subprocess.Popen(data['command'], cwd=data['root'], stdin=subprocess.DEVNULL,
                                    stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        except OSError as error:
            data.update(state='failed', ended=time.time(), exit_code=None, message=str(error))
            log.write(f'[RUN FAILED] Could not start analysis: {error}\n')
        else:
            data.update(pid=proc.pid, pid_identity=process_identity(proc.pid))
            write_json(path, data)
            while True:
                if stopping or Path(data['stop_path']).exists():
                    data['state'] = 'stopping'
                    write_json(path, data)
                    try:
                        os.killpg(proc.pid, signal.SIGTERM)
                    except ProcessLookupError:
                        pass
                    # Do not reap the leader before cleaning up its group: its
                    # children may ignore TERM and its PID must not be reused.
                    time.sleep(5)
                    try:
                        os.killpg(proc.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    code = proc.wait()
                    data['state'] = 'stopped'
                    break
                try:
                    code = proc.wait(timeout=0.2)
                    data['state'] = 'completed' if code == 0 else 'failed'
                    break
                except subprocess.TimeoutExpired:
                    continue
            data.update(exit_code=code, ended=time.time())
            log.write(f"[RUN {data['state'].upper()}] exit {code}\n")
        log.flush()
    write_json(path, data)


if __name__ == '__main__':
    supervise(sys.argv[1])
