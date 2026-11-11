# threshold lambda*(d): smallest lam with max_s max_f J > 1 (float SDP); usage: python3 scan.py d1,d2,...
import sys, time, numpy as np
from sdp import solve, best_s
ds = [int(x) for x in sys.argv[1].split(',')]
lo0, hi0 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.775, float(sys.argv[3]) if len(sys.argv) > 3 else 0.7892
for d in ds:
    t = time.time(); lo, hi = lo0, hi0; sb = 1.0
    while hi - lo > 2e-6:
        mid = (lo+hi)/2
        s, v = best_s(d, mid, lo=0.75, hi=1.35, it=10)
        if v > 1: hi, sb = mid, s
        else: lo = mid
    print(f"d={d}: lambda* in ({lo:.6f}, {hi:.6f}], s={sb:.4f}  ({time.time()-t:.0f}s)", flush=True)
