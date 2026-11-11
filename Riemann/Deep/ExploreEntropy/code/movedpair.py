# movedpair.py (explore19): xi-like function with the zeros gamma_k, gamma_{k+1} merged into a pair m +- i delta (off the line).
# Q_mod(f) = Q_zeta(f) + sum_{+-} [2 Re F(+-m + i delta)^2 - F(+-gamma_k)^2 - F(+-gamma_{k+1})^2]  (real even or odd f).
# For each (k, delta): smallest eigenvalue of the Galerkin matrix of Q_mod on window c (even and odd blocks), bisection for the first
# c where it is negative (lambda_min is non-increasing in c: windows nest).  Usage: movedpair.py k delta cmin cmax [steps]
import sys, time, math, json
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, arb_mat, ctx, fmpq
import mpmath as mp
import wmat

import os
NMUL = float(os.environ.get('NMUL', '1'))
def Nof(c): return max(int(24*NMUL), int(math.ceil(NMUL*2.7*c*math.log(c))) + 8)

def lam_min(cq, k, delta, zeros, sector):
    c = float(cq); bits = int(4*math.pi*c*1.4427) + 96
    ctx.prec = bits; mp.mp.prec = bits
    F = wmat.Form(fmpq(cq), 'zeta')
    N = Nof(c)
    T, idx = F.block(N, sector)
    M = wmat.to_mp(T, mp.mp)
    g1, g2 = zeros[k-1], zeros[k]
    m = mp.mp.nstr((mp.mp.mpf(g1) + mp.mp.mpf(g2))/2, 60)
    ft = F.ft_even if sector == 'even' else F.ft_odd
    def vec(z):
        return [ft(j, z) for j in idx]
    a = vec(acb(arb(str(m)), arb(str(delta))))
    b1 = vec(acb(arb(str(g1)))); b2 = vec(acb(arb(str(g2))))
    n = len(idx)
    D = mp.mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            # Re(a_i a_j): for even sector a real-symmetric product; odd: F_j(t) imaginary on reals -> use F = i*G with G real:
            # generic formula Re[F(z)^2] with F = sum v_j F_j: v^T Re(a a^T) v; F(gamma)^2 -> |F(gamma)|^2 = v^T Re(b conj(b)^T) v
            aa = (a[i]*a[j]); bb1 = (b1[i]*b1[j].conjugate()); bb2 = (b2[i]*b2[j].conjugate())
            # for odd f, F is purely imaginary on R: F(conj z) = -conj F(z); ghat(z) = F(z) conj F(conj z) = -F(z)^2; ghat(gamma) = |F|^2
            sgn = 1 if sector == 'even' else -1
            val = 2*(2*sgn*aa.real) - 2*bb1.real - 2*bb2.real
            D[i, j] = mp.mp.mpf(val.mid().str(40 + int(bits*0.302), radius=False))
    E = mp.mp.eigsy(M + D, eigvals_only=True)
    E0 = mp.mp.eigsy(M, eigvals_only=True)
    return float(min(E)), float(min(E0)), N

if __name__ == '__main__':
    k = int(sys.argv[1]); delta = float(sys.argv[2]); cmin = float(sys.argv[3]); cmax = float(sys.argv[4])
    steps = int(sys.argv[5]) if len(sys.argv) > 5 else 7
    mp.mp.prec = 200
    zeros = [str(mp.mp.zetazero(i).imag) for i in range(1, k+2)]
    zs = [mp.mp.mpf(z) for z in zeros]
    m = (zs[k-1] + zs[k])/2
    print(json.dumps({'k': k, 'gk': str(zs[k-1])[:12], 'gk1': str(zs[k])[:12], 'm': float(m), 'pred_c': float(m/(2*mp.mp.pi)), 'delta': delta}), flush=True)
    zl = zeros
    res = {}
    def run(c):
        t0 = time.time()
        cq = fmpq(int(round(c*64)), 64)
        le, l0e, N = lam_min(cq, k, delta, zl, 'even')
        lo, l0o, _ = lam_min(cq, k, delta, zl, 'odd')
        print(json.dumps({'c': float(cq), 'N': N, 'lmin_even_mod': le, 'lmin_odd_mod': lo, 'lmin_even_zeta': l0e, 'lmin_odd_zeta': l0o, 's': round(time.time()-t0, 1)}), flush=True)
        return min(le, lo)
    lo_c, hi_c = cmin, cmax
    vlo = run(lo_c); vhi = run(hi_c)
    if not (vlo > 0 and vhi < 0):
        print(json.dumps({'note': 'no bracket', 'vlo': vlo, 'vhi': vhi})); sys.exit()
    for _ in range(steps):
        mid = (lo_c + hi_c)/2
        v = run(mid)
        if v > 0: lo_c = mid
        else: hi_c = mid
    print(json.dumps({'threshold_c': [lo_c, hi_c], 'm_over_2pi': float(m/(2*mp.mp.pi))}), flush=True)
