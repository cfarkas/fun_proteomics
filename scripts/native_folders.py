"""Open system file and folder choosers on the analysis computer."""
import os
from pathlib import Path
import shutil
import subprocess
import threading

_DIALOG_LOCK = threading.Lock()


class FolderChooserUnavailable(ValueError):
    """The current session cannot provide a desktop dialog."""


def choose_folder(initial='', target='input'):
    titles = {'output': 'Choose or create output folder', 'input': 'Choose input folder', 'tools': 'Choose tools folder'}
    if target not in titles:
        raise ValueError('Choose an input, output or tools folder.')
    return _choose(initial, titles[target], directory=True)


def choose_file(initial='', title='Choose a file'):
    return _choose(initial, title, directory=False)


def _choose(initial, title, *, directory):
    kind = 'folder' if directory else 'file'
    if not (os.environ.get('DISPLAY') or os.environ.get('WAYLAND_DISPLAY')):
        raise FolderChooserUnavailable('No desktop session. Open the interface on the analysis computer, or enter the path directly.')
    path = Path(initial).expanduser().resolve() if initial else Path.home()
    selected_file = not directory and path.is_file()
    while not path.is_dir() and not selected_file and path != path.parent:
        path = path.parent
    initial_path = str(path) + (os.sep if path.is_dir() else '')
    zenity = shutil.which('zenity')
    kdialog = shutil.which('kdialog')
    if zenity:
        command = [zenity, '--file-selection', *(['--directory'] if directory else []), '--title=' + title, '--filename=' + initial_path]
    elif kdialog:
        command = [kdialog, '--getexistingdirectory' if directory else '--getopenfilename', str(path), '--title', title]
    else:
        raise FolderChooserUnavailable('System chooser unavailable. Install Zenity or KDialog, or enter the path directly.')
    if not _DIALOG_LOCK.acquire(blocking=False):
        raise ValueError('A file or folder chooser is already open. Select a path or cancel that window first.')
    try:
        try:
            result = subprocess.run(command, capture_output=True, text=True, timeout=900)
        except subprocess.TimeoutExpired:
            raise ValueError('System chooser timed out. Click Choose ' + kind + ' to try again.') from None
        if result.returncode == 1:
            return {'cancelled': True}
        if result.returncode != 0:
            raise ValueError('Could not open the system ' + kind + ' chooser. Check the desktop session or enter the path directly.')
        selected = result.stdout.rstrip('\r\n')
        if not selected:
            return {'cancelled': True}
        path = Path(selected).resolve()
        if not (path.is_dir() if directory else path.is_file()):
            raise ValueError('The selected ' + kind + ' no longer exists or has the wrong type.')
        return {'cancelled': False, 'path': str(path)}
    finally:
        _DIALOG_LOCK.release()
