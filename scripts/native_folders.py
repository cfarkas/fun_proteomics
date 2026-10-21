"""Open a folder chooser on the computer running the local interface."""
import os
from pathlib import Path
import shutil
import subprocess
import threading

_DIALOG_LOCK = threading.Lock()

def choose_folder(initial='', target='input'):
    if target not in ('input','output','tools'):
        raise ValueError('Choose an input, output or tools folder.')
    if not (os.environ.get('DISPLAY') or os.environ.get('WAYLAND_DISPLAY')):
        raise ValueError('No desktop session. Open the interface on the analysis computer, or enter the folder path directly.')
    folder=Path(initial).expanduser().resolve() if initial else Path.home()
    while not folder.is_dir() and folder!=folder.parent:
        folder=folder.parent
    title={'output':'Choose or create output folder','input':'Choose input folder','tools':'Choose tools folder'}[target]
    zenity=shutil.which('zenity')
    kdialog=shutil.which('kdialog')
    if zenity:
        command=[zenity,'--file-selection','--directory','--title='+title,'--filename='+str(folder)+os.sep]
    elif kdialog:
        command=[kdialog,'--getexistingdirectory',str(folder),'--title',title]
    else:
        raise ValueError('System folder chooser unavailable. Install Zenity or KDialog, or enter the folder path directly.')
    if not _DIALOG_LOCK.acquire(blocking=False):
        raise ValueError('A folder chooser is already open. Select a folder or cancel that window first.')
    try:
        try:
            result=subprocess.run(command,capture_output=True,text=True,timeout=900)
        except subprocess.TimeoutExpired:
            raise ValueError('Folder chooser timed out. Click Choose folder to try again.') from None
        if result.returncode==1:
            return {'cancelled':True}
        if result.returncode!=0:
            raise ValueError('Could not open the system folder chooser. Check the desktop session or enter the folder path directly.')
        selected=result.stdout.rstrip('\r\n')
        if not selected:
            return {'cancelled':True}
        path=Path(selected).resolve()
        if not path.is_dir():
            raise ValueError('The selected folder no longer exists.')
        return {'cancelled':False,'path':str(path)}
    finally:
        _DIALOG_LOCK.release()
