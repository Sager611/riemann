# hdx_growth.py (explore64): growth of the ceiling with the window. Lam*(c, sqrt c) lower bounds (arb, hdx_ceiling long-range part,
# optimized envelope, N = 4000 cells) against the constant-envelope sum S(c) = sum_{sqrt c <= n < c} 2 w_n (1 - ln n/ln c) (exact) and the
# PNT asymptotic 8 sqrt(c)/ln c of S(c).
import sys, math
sys.path.insert(0, 'Riemann/Deep/ExploreHDX')
from flint import arb, fmpq, ctx
from hdx_ceiling import optimize, certify, atoms
ctx.prec = 128
print('# c, c_link = sqrt c, Lam_long (opt. envelope) >= arb, S(c) constant envelope (arb), 8 sqrt(c)/ln c')
for c in [10, 47, 100, 200, 500, 1000]:
    cq = fmpq(c); cl = math.sqrt(c)
    clq = fmpq(int(math.ceil(cl*10**6)), 10**6)   # c_link >= sqrt c: by monotonicity the bound also holds for links of length ln sqrt(c)
    u, v, t, F = optimize(cq, clq, 4000, longonly=True)
    cv = certify(cq, clq, u, [x for x in t if x[2] == 1], 4000)
    S = arb(0); L = arb(c).log()
    for n, p in atoms(cq):
        if fmpq(n) >= clq: S += 2*arb(p).log()/arb(n).sqrt()*(1 - arb(n).log()/L)
    print(f'{c:5d} {float(clq):9.4f}  {v:8.4f} >= {float(cv.lower()):8.4f}   S {float(S.lower()):8.4f}   {8*math.sqrt(c)/math.log(c):8.4f}', flush=True)
