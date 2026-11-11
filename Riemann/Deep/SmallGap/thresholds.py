# thresholds lambda*(d) for odd d, max over a few s (float SDP); usage: python3 thresholds.py 21,31,41,51,61
import sys, time
from sdp import solve
for d in [int(x) for x in sys.argv[1].split(',')]:
    t = time.time(); best = (9, None)
    for s in (1.0, 1.1, 1.2, 1.3):
        lo, hi = 0.7840, 0.7870
        while hi - lo > 4e-6:
            mid = (lo+hi)/2
            v = solve(d, s, mid)[0]
            if v == v and v > 1: hi = mid
            else: lo = mid
        if hi < best[0]: best = (hi, s)
    print(f'd={d}: lambda* <= {best[0]:.6f} (s={best[1]})  [{time.time()-t:.0f}s]', flush=True)
