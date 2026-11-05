import sys, mpmath as mp
from flint import arb, ctx
from semi import build_L
ctx.prec = 300; mp.mp.dps = 80
Ss = {'inf': set(), 'inf2': {2}, 'inf23': {2, 3}, 'inf235': {2, 3, 5}, 'all': None}
def lm(a, S, N=50):
    M = build_L(arb(2*a), S, N)
    A = mp.matrix(N+1, N+1)
    for i in range(N+1):
        for j in range(N+1): A[i, j] = mp.mpf(M[i, j].mid().str(90, radius=False))
    E = mp.eigsy(A, eigvals_only=True)
    return min(E)
def thr(S, lo, hi, steps=11):
    for _ in range(steps):
        mid = (lo+hi)/2
        if lm(mid, S) > 0: lo = mid
        else: hi = mid
    return lo, hi
for name, lo, hi in [('inf', 0.40, 0.45), ('inf2', 0.55, 0.60), ('inf23', 0.80, 0.90), ('inf235', 0.90, 1.00)]:
    a0, a1 = lm(lo, Ss[name]), lm(hi, Ss[name])
    print(name, 'lm(lo)=%s lm(hi)=%s' % (mp.nstr(a0, 4), mp.nstr(a1, 4)), 'threshold in', thr(Ss[name], lo, hi) if a0 > 0 > a1 else 'n/a', flush=True)
