# sanity: theta = 1 (RH-type input, BGMM Thm 3 / CGdL setting): threshold lambda*(d)
import sys
from sdp import solve
th = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
for d in [int(x) for x in sys.argv[1].split(',')]:
    lo, hi = 0.55, 0.65; best = None
    while hi - lo > 1e-5:
        mid = (lo+hi)/2
        v = max((solve(d, s, mid, theta=th)[0] for s in (0.9, 1.1, 1.3)), key=lambda x: x if x == x else -1)
        if v > 1: hi = mid
        else: lo = mid
    print(f'theta={th} d={d}: lambda* in ({lo:.5f}, {hi:.5f}]', flush=True)
