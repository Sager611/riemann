# ceilaudit9 E30b: LP thresholds a* (W_norm = 0) by bisection, single harmonic, two bumps (optimal by Cauchy-Schwarz), explore30's constraints.
import sys, math
sys.path.insert(0, 'Riemann/Deep/CeilAudit9')
from e30 import single, bounds
def thr(f, lo=0.5, hi=2.5):
    for _ in range(30):
        m = (lo + hi)/2
        if f(m) < 0: hi = m
        else: lo = m
    return hi
D, Sm, S1m = bounds(3e12 + 1e4)
for wf in (0.03, 0.0):
    for rm in (0.0, 1.0):
        a = thr(lambda a: single(a, wf, D, Sm - rm, S1m - rm)); print(f"T=3e12, w={wf}a, extra rounding margin {rm}: a* = {a:.4f}, c* = {math.exp(2*a):.2f}", flush=True)
    a = thr(lambda a: single(a, wf, 1.0, 0.0946*2*math.pi, 0.059*2*math.pi)); print(f"t->inf, w={wf}a: a* = {a:.4f}, c* = {math.exp(2*a):.2f}", flush=True)
for t in (1e13, 1e15, 1e20, 1e30, 1e50, 1e100, 1e300):
    D, Sm, S1m = bounds(t); a = thr(lambda a: single(a, 0.03, D, Sm - 1, S1m - 1)); print(f"t={t:.0e}, w=0.03a, margin 1+1: a* = {a:.4f}, c* = {math.exp(2*a):.2f}", flush=True)
