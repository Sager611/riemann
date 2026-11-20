# ceilaudit9 E30c: anchored LP thresholds (block ends at S = S_1 = 0, explore30's 1-unit margins + rm units for integer rounding).
import sys, math
sys.path.insert(0, 'Riemann/Deep/CeilAudit9')
from e30 import single_anchored, bounds
def thr(f, lo=0.6, hi=2.2):
    for _ in range(14):
        m = (lo + hi)/2
        if f(m) < 0: hi = m
        else: lo = m
    return hi
D, Sm, S1m = bounds(3e12 + 1e4)
for rm in (0.0, 1.0):
    a = thr(lambda a: single_anchored(a, 0.03, D, Sm - rm, S1m - rm)); print(f"T=3e12, w=0.03a, anchored, rounding margin {rm}: a* = {a:.4f}, c* = {math.exp(2*a):.2f}", flush=True)
a = thr(lambda a: single_anchored(a, 0.03, 1.0, 0.0946*2*math.pi, 0.059*2*math.pi)); print(f"t->inf, w=0.03a, anchored: a* = {a:.4f}, c* = {math.exp(2*a):.2f}", flush=True)
a = thr(lambda a: single_anchored(a, 0.0, 1.0, 0.0946*2*math.pi, 0.059*2*math.pi)); print(f"t->inf, w->0, anchored: a* = {a:.4f}, c* = {math.exp(2*a):.2f}", flush=True)
