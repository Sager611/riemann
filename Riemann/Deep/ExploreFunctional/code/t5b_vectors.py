# t5b_vectors.py (explore39): even eigenvectors of Q (CCM basis, wmat, mpmath) -> upper bound on the Poincare margin
# P - R <= (lam1 - lam0)/W(phi1) (t8 notation), and the relative margin (P - R)/R, at c where float64 cannot resolve the cluster.
import sys, math, json, time, mpmath as mp, numpy as np
sys.path.insert(0, '.scratch/explore39'); import wmat
from flint import ctx, fmpq
for c in [float(x) for x in sys.argv[1:]]:
    bits = int(4*math.pi*c*1.4427) + 128; ctx.prec = bits; mp.mp.prec = bits
    N = int(math.ceil(0.6*c*math.log(c))) + 20
    F = wmat.Form(fmpq(int(round(c*64)), 64), 'zeta'); L = math.log(c); a = L/2
    T, idx = F.block(N, 'even'); E, Qv = mp.mp.eigsy(wmat.to_mp(T, mp.mp))
    order = sorted(range(len(idx)), key=lambda i: E[i])
    lam = [E[i] for i in order[:3]]
    vec = lambda k: np.array([float(Qv[r, order[k]]) for r in range(len(idx))])
    Cv = np.array([float(F.pole_vec(j, 'even').mid()) for j in idx])
    xs, ws = np.polynomial.legendre.leggauss(2000); xs = a*xs; ws = a*ws
    Bm = np.array([np.full_like(xs, 1/math.sqrt(L))] + [math.sqrt(2/L)*np.cos(2*math.pi*j*xs/L) for j in idx[1:]]).T
    v0, v1 = vec(0), vec(1); psi = Bm @ v0
    if np.median(psi) < 0: v0 = -v0; psi = -psi
    phi1 = Bm @ v1; Cpsi = float(Cv @ v0); R = 2*Cpsi**2
    m = psi > 1e-10*psi.max()   # dropping the unresolved edge zone lowers W(phi1): a valid lower bound, so P-R bound stays an upper bound
    W1 = (np.sum((ws*np.cosh(xs/2)/np.where(m, psi, 1.0)*phi1**2)[m]) - float(Cv @ v1)**2/Cpsi)/Cpsi
    gap = lam[1] - lam[0]
    print(json.dumps({'c': c, 'N': N, 'lam0': mp.nstr(lam[0], 6), 'lam1': mp.nstr(lam[1], 6), 'psi>0 (|x|<=0.98a)': bool(np.all(psi[np.abs(xs) <= 0.98*a] > 0)), 'mask_frac': float('%.3g' % (1 - m.mean())),
                      'psi_min/max': float('%.3g' % (psi.min()/psi.max())), 'R': round(R, 6), 'W(phi1)': float('%.4g' % W1),
                      'P-R<=': mp.nstr(gap/W1, 5), '(P-R)/R<=': mp.nstr(gap/W1/R, 5)}), flush=True)
