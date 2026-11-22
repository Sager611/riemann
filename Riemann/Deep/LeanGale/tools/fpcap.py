#!/usr/bin/env python3
"""fpcap: run a command; kill its whole process tree when any process's physical memory footprint
(the metric of `top -l 1 -pid <pid> -stats mem`, "Physical memory footprint of the process") exceeds CAP_MB,
or when the wall time exceeds LIMIT_S.  Uses proc_pid_rusage(RUSAGE_INFO_V4).ri_lifetime_max_phys_footprint,
so short spikes between samples are caught too.  Usage: python3 fpcap.py CAP_MB LIMIT_S cmd [args...]
Prints 'fpcap: peak N MB (lifetime max footprint) ...' on stderr; exit 137 when killed."""
import ctypes, os, signal, subprocess, sys, time

lib = ctypes.CDLL('/usr/lib/libproc.dylib')
lib.proc_pid_rusage.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_void_p]
lib.proc_listchildpids.argtypes = [ctypes.c_int, ctypes.c_void_p, ctypes.c_int]

def maxfp(pid):
    buf = (ctypes.c_uint64 * 48)()
    if lib.proc_pid_rusage(pid, 4, ctypes.byref(buf)) != 0:
        return 0
    # layout: uuid (16 bytes = 2 words), then 64-bit fields; ri_lifetime_max_phys_footprint is field 28
    return max(buf[2 + 28], buf[2 + 7])

def children(pid):
    arr = (ctypes.c_int * 4096)()
    n = lib.proc_listchildpids(pid, ctypes.byref(arr), ctypes.sizeof(arr))
    if n <= 0:
        return []
    n = min(n, 4096)
    return [p for p in arr[:n] if p > 0]

def tree(pid):
    out, todo = [], [pid]
    while todo:
        p = todo.pop()
        out.append(p)
        todo.extend(children(p))
    return out

def main():
    cap = int(sys.argv[1]) * 1024 * 1024
    lim = float(sys.argv[2])
    cmd = sys.argv[3:]
    t0 = time.time()
    proc = subprocess.Popen(cmd, start_new_session=True)
    peak = 0
    status = ''
    seen = {}
    while True:
        rc = proc.poll()
        for p in tree(proc.pid):
            m = maxfp(p)
            if m:
                seen[p] = max(seen.get(p, 0), m)
        cur = max(seen.values(), default=0)
        peak = max(peak, cur)
        if rc is not None:
            break
        if cur > cap:
            status = 'fpcap: killed at %d MB > %d MB' % (cur >> 20, cap >> 20)
            break
        if time.time() - t0 > lim:
            status = 'fpcap: killed after %d s' % lim
            break
        time.sleep(0.2)
    if status:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.wait()
        print('%s (peak %d MB)' % (status, peak >> 20), file=sys.stderr)
        sys.exit(137)
    print('fpcap: peak %d MB (lifetime max footprint), rc %d, %d s' % (peak >> 20, rc, time.time() - t0), file=sys.stderr)
    sys.exit(rc)

main()
