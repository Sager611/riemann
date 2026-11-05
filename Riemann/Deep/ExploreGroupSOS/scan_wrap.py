# scan_wrap.py: like scan_ext but P fixed (= fac*2L) and extension atoms at k log p mod P (wrapped into the gap arc).
import sys, math, numpy as np
from scipy.optimize import linprog
from weil_lat import atoms, m_sym, hplus
def beta_p(p, c, ks_free, ngrid=4096):
    ks_fixed = [k for k in range(1, 200) if p**k < c]
    th = np.linspace(0, math.pi, ngrid)
    fixed = -2 * sum((math.log(p) / p**(k / 2)) * np.cos(k * th) for k in ks_fixed)
    if not ks_free: return float(-fixed.min()), []
    nf = len(ks_free)
    A = np.hstack([-np.ones((ngrid, 1)), -2 * np.cos(np.outer(th, ks_free))])
    res = linprog(np.concatenate([[1.0], np.zeros(nf)]), A_ub=A, b_ub=fixed, bounds=[(None, None)] * (nf + 1), method='highs')
    return float(res.x[0]), [(k, float(e)) for k, e in zip(ks_free, res.x[1:])]
def count(c, P, E, sig, L, at):
    r = 20.0
    while True:
        env = hplus(r) - 2 * sum(math.exp(-(2*n+0.5)*L) * ((2*n+0.5) + r) / ((2*n+0.5)**2 + r*r) for n in range(40)) - 4 * math.exp(L/2) / (r - 0.5) - sig
        if env > 0: break
        r *= 1.02
    K = int(r * P / (2 * math.pi)) + 2
    neg = 0; mx = 0.0
    for k0 in range(0, K + 1, 200000):
        k = np.arange(k0, min(K + 1, k0 + 200000)); rr = 2 * math.pi * k / P
        w = m_sym(rr, L, at)
        for x, e in E: w = w + 2 * e * np.cos(rr * x)
        neg += int(np.sum(w < 0))
        if np.any(w < 0): mx = max(mx, float(rr[w < 0].max()))
    return K, neg, mx
c = float(sys.argv[1]); extra = int(sys.argv[2]); kmax = int(sys.argv[3])
L = math.log(c); at = atoms(c)
primes = [p for p in range(2, int(c) + 1) if all(p % q for q in range(2, int(p**0.5) + 1)) and p < c]
for fac in [1.0, 1.25, 1.5, 2.0]:
    P = 2 * L * fac
    K0, n0, m0 = count(c, P, [], 2 * sum(w for _, w in at), L, at)
    if fac == 1.0:
        print(f'c={c:g} P=2L (Chekmachev lattice): nodes {K0}, negative {n0}, largest negative r {m0:.0f}', flush=True); continue
    E = []; btot = 0
    for p in primes:
        k0 = max([k for k in range(1, 200) if p**k < c]) + 1
        ks = [k for k in range(k0, kmax + 1) if L <= (k * math.log(p)) % P <= P - L][:extra]
        b, co = beta_p(p, c, ks); btot += b
        E += [(k * math.log(p), e) for k, e in co]
    K1, n1, m1 = count(c, P, E, btot, L, at)
    print(f'c={c:g} P={fac}*2L: no ext nodes {K0} neg {n0} (max r {m0:.0f}) | wrapped sep ext ({len(E)} atoms) envelope {btot:.3f} vs 2sigma {2*sum(w for _, w in at):.3f}: nodes {K1}, neg {n1}, max neg r {m1:.0f}', flush=True)
