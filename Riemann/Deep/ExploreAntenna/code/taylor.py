# taylor.py (explore19): Taylor n-bar / Schelkunoff null placement in the Woodward-Lawson (CCM Fourier) basis.
# For n-bar = 1..Nmax+1 the unique (up to scale) even v in span{phi_0..phi_{nbar-1}} whose transform vanishes at gamma_1..gamma_{nbar-1}
# (Taylor's line source with its first nbar-1 nulls moved onto the zeta zeros; the remaining nulls stay on the Nyquist lattice 2 pi j/L).
# Prints RQ(nbar) = Q(v)/|v|^2 against the Galerkin minimum.  Usage: taylor.py c
import sys, math, json, time
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, arb_mat, ctx, fmpq
import mpmath as mp
import wmat

c = float(sys.argv[1])
bits = int(4*math.pi*c*1.4427*2) + 160      # null systems are Cauchy-like: double precision budget
ctx.prec = bits; mp.mp.prec = bits
Nmax = max(24, int(math.ceil(2.7*c*math.log(c))) + 8)
cq = fmpq(int(round(c*64)), 64)
F = wmat.Form(cq, 'zeta')
T, idx = F.block(Nmax, 'even')
M = wmat.to_mp(T, mp.mp)
E = mp.mp.eigsy(M, eigvals_only=True)
lmin = min(E)
L = mp.mp.log(mp.mp.mpf(float(cq)))
w = [2*mp.mp.pi*j/L for j in range(Nmax+1)]
Nj = [mp.mp.sqrt(1/L)] + [mp.mp.sqrt(2/L)]*Nmax
t0 = time.time()
zs = [mp.mp.zetazero(k).imag for k in range(1, Nmax+1)]
res = []
for nb in range(1, Nmax+2):
    if nb == 1:
        v = mp.mp.matrix([[1]])
    else:
        A = mp.mp.matrix(nb-1, nb)
        for k in range(nb-1):
            g = zs[k]
            for j in range(nb):
                A[k, j] = Nj[j]*(-1)**j/(g*g - w[j]*w[j])
        # null vector: fix v_0 = 1, solve the rest
        B = A[:, 1:nb]; r = -A[:, 0]
        x = mp.mp.lu_solve(B, r)
        v = mp.mp.matrix(nb, 1); v[0] = 1
        for j in range(1, nb): v[j] = x[j-1]
    num = mp.mp.mpf(0)
    for i in range(nb):
        for j in range(nb): num += v[i]*M[i, j]*v[j]
    den = sum(v[i]**2 for i in range(nb))
    res.append((nb, num/den))
best = min(res, key=lambda r: r[1])
Tstar = 2*math.pi*c
NT = (Tstar/(2*math.pi))*math.log(Tstar/(2*math.pi*math.e)) + 7/8
out = {'c': c, 'Nmax': Nmax, 'bits': bits, 'lambda_min_galerkin': mp.mp.nstr(lmin, 6), 'best_nbar': best[0], 'best_RQ': mp.mp.nstr(best[1], 6),
       'ratio_best_over_lmin': mp.mp.nstr(best[1]/lmin, 4), 'N(2 pi c)': round(NT, 2),
       'RQ_by_nbar': {r[0]: mp.mp.nstr(r[1], 3) for r in res if r[0] % max(1, Nmax//24) == 0 or r[0] == best[0]}, 's': round(time.time()-t0, 1)}
print(json.dumps(out), flush=True)
