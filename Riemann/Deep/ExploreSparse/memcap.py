# explore31: macOS physical footprint (what `top -stats mem` shows) and resident size of this process, in MB.
import ctypes, os
_lib = ctypes.CDLL('/usr/lib/libproc.dylib')
class _RI(ctypes.Structure):
    _fields_ = [('uuid', ctypes.c_uint8 * 16)] + [(n, ctypes.c_uint64) for n in (
        'user', 'sys', 'idle', 'intr', 'pageins', 'wired', 'resident', 'footprint', 'start', 'exit',
        'cuser', 'csys', 'cidle', 'cintr', 'cpageins', 'celapsed', 'dread', 'dwrite')]
def mem():
    ri = _RI()
    _lib.proc_pid_rusage(os.getpid(), 2, ctypes.byref(ri))
    return ri.footprint / 1e6, ri.resident / 1e6
def check(cap=280, where=''):
    f, r = mem()
    if f > cap:
        raise SystemExit(f'memory cap {cap} MB exceeded ({f:.0f} MB footprint) at {where}')
    return f
