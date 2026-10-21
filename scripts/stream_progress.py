#!/usr/bin/env python3
"""Stream local task logs without changing scientific scripts or cache inputs."""
import argparse
from datetime import datetime
import os
from pathlib import Path
import re
import signal
import time


def alive(pid):
    try:
        os.kill(pid, 0)
        status = Path(f'/proc/{pid}/stat').read_text().rsplit(')', 1)[1].split()[0]
        return status != 'Z'
    except (ProcessLookupError, FileNotFoundError):
        return False


class Monitor:
    def __init__(self, work, since, verbosity=4):
        self.verbosity = verbosity
        self.work = work
        self.since = since
        self.tasks = {}

    def poll(self):
        for marker in self.work.glob('*/*/.command.begin'):
            task = marker.parent
            if task not in self.tasks and marker.stat().st_mtime >= self.since:
                script = task / '.command.run'
                header = script.read_text(errors='replace') if script.exists() else ''
                match = re.search(r"^### name: '(.+)'$", header, re.M)
                name = match[1] if match else task.name
                now = time.monotonic()
                self.tasks[task] = {'name': name, 'offsets': {}, 'last': now, 'heartbeat': now, 'activity': '', 'done': False, 'pending': {}}
                if self.verbosity >= 2:
                    print(f'[START] {name} | work directory: {task}', flush=True)
                if self.verbosity >= 4:
                    command = task / '.command.sh'
                    if command.exists():
                        print(f'[COMMAND] {name}\n{command.read_text(errors="replace")}', flush=True)
        for task, state in self.tasks.items():
            if state['done']:
                continue
            for rel in ['.command.out', '.command.err', 'run/pipeline.log']:
                path = task / rel
                try:
                    offset = state['offsets'].get(rel, 0)
                    if path.stat().st_size < offset:
                        offset = 0
                    with path.open('rb') as handle:
                        handle.seek(offset)
                        chunk = handle.read()
                        state['offsets'][rel] = handle.tell()
                    if chunk:
                        if self.verbosity >= 4:
                            print(f"[{state['name']} | {rel}]\n" + chunk.decode('utf-8', errors='replace'), end='', flush=True)
                        else:
                            data = state['pending'].get(rel, b'') + chunk
                            lines = data.split(b'\n')
                            state['pending'][rel] = lines.pop()
                            for line in lines:
                                self.emit_line(state['name'], line.decode('utf-8', errors='replace'))
                        state['last'] = time.monotonic()
                        activity = chunk.decode('utf-8', errors='replace').strip().splitlines()
                        if activity:
                            state['activity'] = activity[-1].strip()[:200]
                except FileNotFoundError:
                    continue
            exitfile = task / '.exitcode'
            if exitfile.exists() and exitfile.read_text().strip():
                state['done'] = True
                for remainder in state['pending'].values():
                    if remainder:
                        self.emit_line(state['name'], remainder.decode('utf-8', errors='replace'))
                code = exitfile.read_text().strip()
                if self.verbosity >= 2 or code != '0':
                    print(f"[END] {state['name']} | exit={code}", flush=True)
                if code != '0' and self.verbosity < 4:
                    error_log = task / '.command.err'
                    if error_log.exists():
                        print('\n'.join(error_log.read_text(errors='replace').splitlines()[-50:]), flush=True)
            elif self.verbosity >= 2 and time.monotonic() - max(state['last'], state['heartbeat']) >= 30:
                quiet = int(time.monotonic() - state['last'])
                duration = f'{quiet // 60}m {quiet % 60}s' if quiet >= 60 else f'{quiet}s'
                activity = f" Last output: {state['activity']}" if state['activity'] else ''
                print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] [HEARTBEAT] {state['name']} has not exited; no task log output for {duration}.{activity}", flush=True)
                state['heartbeat'] = time.monotonic()

    def emit_line(self, name, line):
        warning = re.search(r'\b(?:WARN(?:ING)?|ERROR|FATAL|Exception)\b', line, re.I)
        stage = re.search(r'Starting|Running|Executing|Finished|Completed|Processing (?:file|run|sample)|Matching-between-runs|Transferring|MSFragger|MSBooster|Percolator|ProteinProphet|Philosopher|IonQuant|progress|\d+(?:\.\d+)?%', line, re.I)
        if warning or (self.verbosity >= 3 and stage):
            print(f'[{name}] {line}', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent', type=int, required=True)
    parser.add_argument('--work-dir', type=Path, required=True)
    parser.add_argument('--since', type=float, required=True)
    parser.add_argument('--verbosity', type=int, choices=range(1, 5), default=3)
    args = parser.parse_args()
    monitor = Monitor(args.work_dir, args.since, args.verbosity)
    stopped = False

    def stop(signum, frame):
        nonlocal stopped
        stopped = True

    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGINT, stop)
    while not stopped and alive(args.parent):
        monitor.poll()
        time.sleep(0.5)
    monitor.poll()


if __name__ == '__main__':
    try:
        main()
    except BrokenPipeError:
        pass
