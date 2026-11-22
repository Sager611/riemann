# Samples the physical footprint and RSS of the process whose full command line matches argv[1] (ERE, pgrep -f).
import ctypes, os, subprocess, sys, time
lib = ctypes.CDLL('/usr/lib/libproc.dylib')
F = 'user_time system_time pkg_idle_wkups interrupt_wkups pageins wired_size resident_size phys_footprint proc_start_abstime proc_exit_abstime child_user_time child_system_time child_pkg_idle_wkups child_interrupt_wkups child_pageins child_elapsed_abstime diskio_bytesread diskio_byteswritten qos_default qos_maintenance qos_background qos_utility qos_legacy qos_user_initiated qos_user_interactive billed_system_time serviced_system_time logical_writes lifetime_max_phys_footprint instructions cycles billed_energy serviced_energy interval_max_phys_footprint runnable_time'.split()
class R(ctypes.Structure):
    _fields_ = [('uuid', ctypes.c_uint8 * 16)] + [(f, ctypes.c_uint64) for f in F]
t0 = time.time(); pid = None
while pid is None and time.time() - t0 < 30:
    ps = [int(p) for p in subprocess.run(['pgrep', '-f', sys.argv[1]], capture_output=True, text=True).stdout.split()]
    ps = [p for p in ps if p != os.getpid()]
    if ps: pid = max(ps)
    else: time.sleep(0.01)
if pid is None: sys.exit('fp: no process matched')
rss = fp = n = 0; t1 = time.time()
while True:
    r = R()
    if lib.proc_pid_rusage(pid, 4, ctypes.byref(r)) != 0: break
    n += 1; rss = max(rss, r.resident_size); fp = max(fp, r.lifetime_max_phys_footprint)
    time.sleep(0.02)
print(f'fp: {n} samples over {time.time() - t1:.1f} s: max RSS {rss >> 20} MB, max physical footprint {fp >> 20} MB')
