# t3b_bisect.py -- refine the DH threshold h_A between 0.29 and 0.30 (imports t3_lagarias helpers)
import importlib.util, sys, math, io, contextlib
spec = importlib.util.spec_from_file_location('t3', 't3_lagarias.py')
src = open('t3_lagarias.py').read().split("t0 = time.time()")[0]
g = {}
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(src, 't3_lagarias.py', 'exec'), g)
lo, hi = 0.29, 0.30
for it in range(7):
    mid = (lo + hi) / 2
    c, mx, mn = g['count'](g['logPhi'], mid, 0.5005, 1.25, 80.0, 92.0, step=0.001)
    print("h=%.6f zeros right of line: %.3f (max step %.3f)" % (mid, c, mx)); sys.stdout.flush()
    if round(c) >= 1: lo = mid
    else: hi = mid
print("DH local threshold h_A in (%.5f, %.5f); rho0 - 1/2 = 0.30852" % (lo, hi))
