# run5.py (explore17): deflation test. V_K = span{Phi^{(j)} 1_W} (even j for the even block, odd j for the odd block),
# j exact via the recursion (2 z d/dz) on z^a e^{-b z}, z = e^{2u}. Report lambda_min of Q restricted to the L2-complement of V_K
# inside E_N (Q from WeilMat Arb closed forms, rounded to float64 after assembly).
import os, sys, time, math, numpy as np, mpmath as mp
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'WeilSimple'))
from flint import arb, arb_mat, ctx, fmpq
from weil_matrix import WeilMat
t0 = time.time(); mp.mp.dps = 50
NT = 10   # terms n of the Phi series
def phi_derivs(u, J):   # [Phi^{(j)}(u) for j = 0..J]
    z = mp.e**(2*mp.mpf(u)); out = [mp.mpf(0)]*(J + 1)
    for n in range(1, NT + 1):
        b = mp.pi*n*n
        for (cf, a) in ((2*mp.pi**2*n**4, mp.mpf(9)/4), (-3*mp.pi*n*n, mp.mpf(5)/4)):
            poly = [cf]   # coefficients of z^{a+m} e^{-bz}
            for j in range(J + 1):
                out[j] += mp.e**(-b*z)*z**a*mp.fsum(poly[m]*z**m for m in range(len(poly)))
                new = [mp.mpf(0)]*(len(poly) + 1)
                for m, cm in enumerate(poly):
                    new[m] += 2*(a + m)*cm; new[m + 1] += -2*b*cm
                poly = new
    return out
def basis_vals(L, sector, N, x):
    idx = np.arange(0, N + 1) if sector == 'even' else np.arange(1, N + 1)
    w = 2*math.pi*idx/L
    if sector == 'even':
        Nj = np.where(idx == 0, 1/math.sqrt(L), math.sqrt(2/L)); return np.array([Nj[a]*np.cos(w[a]*x) for a in range(len(idx))])
    return np.array([math.sqrt(2/L)*np.sin(w[a]*x) for a in range(len(idx))])
