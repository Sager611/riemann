#!/usr/bin/env python3
"""Sample RSS and physical footprint (top's MEM) of one process every 20 ms until it exits.
usage: python3 fpsamp.py <pgrep -f pattern>   (waits up to 10 s for the process)"""
import ctypes, ctypes.util, subprocess, sys, time
libc = ctypes.CDLL(ctypes.util.find_library('c'))
class RU(ctypes.Structure):  # struct rusage_info_v0 (sys/resource.h)
    _fields_ = [('uuid', ctypes.c_uint8 * 16), ('user', ctypes.c_uint64), ('system', ctypes.c_uint64),
                ('pkg_idle', ctypes.c_uint64), ('intr', ctypes.c_uint64), ('pageins', ctypes.c_uint64),
                ('wired', ctypes.c_uint64), ('resident', ctypes.c_uint64), ('footprint', ctypes.c_uint64),
                ('start', ctypes.c_uint64), ('exit', ctypes.c_uint64)]
def sample(pid):
    ru = RU()
    if libc.proc_pid_rusage(pid, 0, ctypes.byref(ru)) != 0:
        return None
    return ru.resident >> 20, ru.footprint >> 20
pid, t0 = None, time.time()
while pid is None and time.time() - t0 < 10:
    out = subprocess.run(['pgrep', '-f', sys.argv[1]], capture_output=True, text=True).stdout.split()
    pid = int(out[0]) if out else None
    time.sleep(0.005)
log, t1 = [], time.time()
while pid is not None:
    s = sample(pid)
    if s is None:
        break
    log.append((round(time.time() - t1, 2),) + s)
    time.sleep(0.02)
if log:
    print('max rss=%d MB max fp=%d MB after %.2f s (last sample: rss %d, fp %d)' % (
        max(e[1] for e in log), max(e[2] for e in log), log[-1][0], log[-1][1], log[-1][2]))
else:
    print('no samples')

