# env.py: high-frequency envelope constants of the lattice symbol.
#  2sigma  = sum_{n<c} 2 w_n                  (no extension; Chekmachev's scan needs log(r/2pi) > 2 sigma)
#  beta_sep = sum_p beta_p, per-prime Fejer-Riesz optimum with extension atoms at k log p (p^k >= c, k <= K_p)
#  normC   = ||C_p|| on L^2(-L/2, L/2) (Perron eigenvalue, discretized) = lower bound for ANY extension's envelope
import sys, math, numpy as np
from scipy.optimize import linprog
from weil_lat import atoms
def beta_p(p, c, extra, ngrid=2048):
    L = math.log(c)
    ks_fixed = [k for k in range(1, 200) if p**k < c]
    k0 = (ks_fixed[-1] + 1) if ks_fixed else 1
    ks_free = list(range(k0, k0 + extra))
    th = np.linspace(0, math.pi, ngrid)
    fixed = -2 * sum((math.log(p) / p**(k / 2)) * np.cos(k * th) for k in ks_fixed)
    if not ks_free:
        return float(-fixed.min()), k0, k0 - 1
    # variables: beta, e_k ; minimize beta s.t. beta + fixed + 2 sum e_k cos(k th) >= 0
    nf = len(ks_free)
    A = np.hstack([-np.ones((ngrid, 1)), -2 * np.cos(np.outer(th, ks_free))])
    b = fixed
    res = linprog(np.concatenate([[1.0], np.zeros(nf)]), A_ub=A, b_ub=b, bounds=[(None, None)] * (nf + 1), method='highs')
    return float(res.x[0]) * (1 + 1e-6) + 1e-9, k0, ks_free[-1]
def normC(c, ncell=6000):
    L = math.log(c); at = atoms(c); h = L / ncell
    v = np.ones(ncell)
    # shift operator on cells: (S_s f)(x) = f(x - s); with s = q h + r h (0<=r<1): average of two neighbouring shifts
    def apply(f):
        out = np.zeros(ncell)
        for ln, w in at:
            for s in (ln, -ln):
                q = math.floor(s / h); fr = s / h - q
                for qq, ww in ((q, 1 - fr), (q + 1, fr)):
                    if ww == 0: continue
                    if qq >= 0:
                        if qq < ncell: out[qq:] += w * ww * f[:ncell - qq]
                    else:
                        if -qq < ncell: out[:ncell + qq] += w * ww * f[-qq:]
        return out
    lam = 0
    for it in range(300):
        u = apply(v); lam_new = np.linalg.norm(u) / np.linalg.norm(v); v = u / np.linalg.norm(u)
        if abs(lam_new - lam) < 1e-10: break
        lam = lam_new
    return lam_new
print('c | 2sigma | beta_sep (extra=2,4,8 free powers) | ||C_p||_window | log10 scan nodes: P=2L,no ext | sep ext(8) | norm-envelope')
for c in [10, 16, math.exp(3.2), 47, 71, 200]:
    at = atoms(c); L = math.log(c)
    s2 = 2 * sum(w for _, w in at)
    primes = sorted(set(round(math.exp(ln)) for ln, _ in at if all(round(math.exp(ln)) % q for q in range(2, int(math.exp(ln)**0.5) + 1))))
    bs = []
    for extra in (2, 4, 8):
        tot = 0; Pneed = 2 * L
        for p in primes:
            b, k0, kmax = beta_p(p, c, extra)
            tot += b; Pneed = max(Pneed, 2 * kmax * math.log(p))
        bs.append((tot, Pneed))
    nc = normC(c)
    lg = lambda P, beta: math.log10(P / (2 * math.pi) * 2 * math.pi * math.exp(beta + 1))
    print(f'{c:.3g} | {s2:.3f} | ' + ', '.join(f'{b:.3f} (P>={P:.1f})' for b, P in bs) + f' | {nc:.3f} | {lg(2*L, s2):.1f} | {lg(bs[2][1], bs[2][0]):.1f} | {lg(2*L, nc):.1f}', flush=True)
