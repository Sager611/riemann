# certdel.py (morirate): certified finite-mu instance of the lost-sample law. For the explicit even vector v (ground state of the
# Galerkin matrix of W - 2|F(gamma_1)|^2), Arb gives an upper bound X for the Rayleigh quotient of the form W - 2|u^(gamma_1)|^2
# (W = Weil form via the closed-form Galerkin matrix of wmat.py, entries as Arb balls; F_j(gamma) as balls for gamma in a certified
# interval around the first zero). Then lambda_min(W, even) / lambda_min(W - 2|u^(gamma_1)|^2) >= eps_c / X with eps_c certified
# (Chekmachev; re-certified in ../ChekAudit). Usage: certdel.py c k
import sys, os, json, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'ExploreAntenna', 'code'))
from flint import arb, acb, arb_mat, ctx, fmpq
import mpmath as mp
import wmat

def Zfun(t):
    # Hardy Z(t) = exp(i theta(t)) zeta(1/2 + i t), theta(t) = Im lgamma(1/4 + i t/2) - (t/2) log pi
    s = acb(arb(1)/2, t)
    th = acb(arb(1)/4, t/2).lgamma().imag - t/2*arb.pi().log()
    z = s.zeta()*acb(th.cos(), th.sin())
    return z.real, z.imag

def main():
    c = float(sys.argv[1]); k = int(sys.argv[2])
    bits = int(4*math.pi*c*1.4427) + 260
    ctx.prec = bits; mp.mp.prec = bits; dps = int(bits*0.302) + 10
    mp.mp.dps = dps
    g = mp.zetazero(k).imag
    gs = mp.nstr(g, dps)
    # certify a sign change of Z on [g - r, g + r]
    r = arb(10)**(-60)
    zl, _ = Zfun(arb(gs) - r); zr, _ = Zfun(arb(gs) + r)
    sign_change = bool((zl < 0 and zr > 0) or (zl > 0 and zr < 0))
    gball = arb(gs).union(arb(gs) - r).union(arb(gs) + r)
    mp.mp.prec = bits
    F = wmat.Form(fmpq(int(round(c*64)), 64), 'zeta')
    N = max(24, int(math.ceil(2.7*c*math.log(c))) + 8)
    T, idx = F.block(N, 'even'); n = len(idx)
    b = [F.ft_even(j, acb(gball)).real for j in idx]
    M = wmat.to_mp(T, mp.mp)
    bm = [mp.mp.mpf(x.mid().str(dps, radius=False)) for x in b]
    A = M.copy()
    for i in range(n):
        for j in range(n): A[i, j] -= 2*bm[i]*bm[j]
    E, V = mp.mp.eigsy(A)
    i0 = min(range(n), key=lambda i: E[i])
    v = [arb(mp.nstr(V[rr, i0], dps)) for rr in range(n)]
    vTv = sum((x*x for x in v), arb(0))
    vTv_W = arb(0)
    for i in range(n):
        row = arb(0)
        for j in range(n): row += T[i, j]*v[j]
        vTv_W += v[i]*row
    Fv = sum((v[j]*b[j] for j in range(n)), arb(0))
    Q = vTv_W - 2*Fv*Fv
    X = Q/vTv
    print(json.dumps({'c': c, 'k': k, 'gamma': float(g), 'y': float(g)/(2*math.pi*c), 'sign_change_certified': sign_change, 'N': N,
                      'rayleigh_W_minus_zero': X.str(8), 'upper': X.upper().str(6), 'galerkin_min': mp.nstr(E[i0], 8),
                      'pred_4g2': 4*math.log((2 + math.sqrt(4 - (float(g)/(2*math.pi*c))**2))/(float(g)/(2*math.pi*c)))}), flush=True)

if __name__ == '__main__':
    main()
