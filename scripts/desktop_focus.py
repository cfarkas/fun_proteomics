"""Focus the local study browser only after an explicit user input gesture."""
import ctypes as c
import os
import re
import subprocess

def focus_browser():
    if not os.environ.get('DISPLAY'):return {'focused':False,'reason':'No local X11 desktop'}
    try:
        windows=subprocess.check_output(['xprop','-root','_NET_CLIENT_LIST_STACKING'],text=True,timeout=2)
        target=None
        for ident in reversed(re.findall(r'0x[0-9a-fA-F]+',windows)):
            props=subprocess.check_output(['xprop','-id',ident,'WM_CLASS','_NET_WM_NAME'],text=True,timeout=2)
            if 'Proteomics · Study setup' in props and any(name in props.lower() for name in ('firefox','chromium','google-chrome')):
                target=int(ident,16);break
        if target is None:return {'focused':False,'reason':'Study browser window not found'}
        x=c.CDLL('libX11.so.6')
        x.XOpenDisplay.restype=c.c_void_p
        x.XDefaultRootWindow.argtypes=[c.c_void_p];x.XDefaultRootWindow.restype=c.c_ulong
        x.XInternAtom.argtypes=[c.c_void_p,c.c_char_p,c.c_int];x.XInternAtom.restype=c.c_ulong
        class Data(c.Union):_fields_=[('l',c.c_long*5)]
        class Client(c.Structure):_fields_=[('type',c.c_int),('serial',c.c_ulong),('send_event',c.c_int),('display',c.c_void_p),('window',c.c_ulong),('message_type',c.c_ulong),('format',c.c_int),('data',Data)]
        class Event(c.Union):_fields_=[('client',Client),('pad',c.c_long*24)]
        x.XSendEvent.argtypes=[c.c_void_p,c.c_ulong,c.c_int,c.c_long,c.POINTER(Event)]
        x.XRaiseWindow.argtypes=[c.c_void_p,c.c_ulong]
        x.XSetInputFocus.argtypes=[c.c_void_p,c.c_ulong,c.c_int,c.c_ulong]
        x.XGetInputFocus.argtypes=[c.c_void_p,c.POINTER(c.c_ulong),c.POINTER(c.c_int)]
        x.XSync.argtypes=[c.c_void_p,c.c_int];x.XCloseDisplay.argtypes=[c.c_void_p]
        d=x.XOpenDisplay(None)
        if not d:return {'focused':False,'reason':'Cannot access local desktop'}
        try:
            e=Event();e.client.type=33;e.client.display=d;e.client.window=target;e.client.message_type=x.XInternAtom(d,b'_NET_ACTIVE_WINDOW',0);e.client.format=32;e.client.data.l[0]=2
            x.XSendEvent(d,x.XDefaultRootWindow(d),0,(1<<20)|(1<<19),c.byref(e))
            x.XRaiseWindow(d,target)
            x.XSetInputFocus(d,target,2,0)
            x.XSync(d,0)
            focused=c.c_ulong();revert=c.c_int();x.XGetInputFocus(d,c.byref(focused),c.byref(revert))
            return {'focused':focused.value==target,'window':hex(target)}
        finally:x.XCloseDisplay(d)
    except (OSError,subprocess.SubprocessError) as error:
        return {'focused':False,'reason':str(error)}

if __name__=='__main__':
    import json
    print(json.dumps(focus_browser()))
