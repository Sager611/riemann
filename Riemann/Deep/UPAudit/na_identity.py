# N-A: exact check of Lemma Z (Z2)+(Z3) bookkeeping and the zero-side formula, real zeta zeros.
# M(s) = sum over squarefree n (primes in P) mu(n) f(log n/L) prod rho(log q/L) n^{-s}, s = 1/2 +- y + it.
# A: Newton power sums (collision packets q^h, coefficient rho^h q^{-hs}).
# B: free primes -> (nu - nu_pp)^{(x)m}; prime-power packets q^h with rho(h log q/L)/h (q down to T^{eta/h}); G_r direct.
# C: as B with G_r from the zeros: G_r(c) = int f((c+v)/L) g^{*r}(v) dv, g = -phi(u) sum_rho e^{(rho-s)u}, |gamma-t| <= Gam.
import numpy as np, math, sys, json, itertools
from math import log
ZEROS = np.load('/Users/adriansager/playground/riemann/Riemann/Deep/ExploreHelio/results/zeta_zeros.npy')
def S(x):
    x = np.clip(np.asarray(x, float), 0.0, 1.0)
    g = lambda y: np.where(y > 0, np.exp(-1/np.maximum(y, 1e-300)), 0.0)
    return g(x)/(g(x)+g(1-x))
def primes_upto(n):
    s = np.ones(n+1, bool); s[:2] = False
    for p in range(2, int(n**0.5)+1):
        if s[p]: s[p*p::p] = False
    return np.nonzero(s)[0]
def partitions(j, maxpart=None):
    if maxpart is None: maxpart = j
    if j == 0: yield (); return
    for p in range(min(j, maxpart), 0, -1):
        for rest in partitions(j-p, p): yield (p,)+rest
def zlam(lam):
    z = 1
    for h in set(lam):
        m = lam.count(h); z *= h**m*math.factorial(m)
    return z
class Case:
    def __init__(c, T, th, eta, b):
        c.T, c.th, c.eta, c.b, c.L = T, th, eta, b, log(T); c.lX = th*c.L
        c.P = [int(q) for q in primes_upto(int(T**b)+2) if c.rho(log(q)/c.L) > 0]
        pp = []
        for q in primes_upto(int(T**b)+2):
            h, n = 1, int(q)
            while log(n) < b*c.L:
                w = c.rho(log(n)/c.L)/h
                if w > 0: pp.append((n, h, float(w)))
                h += 1; n *= int(q)
        pp.sort(key=lambda x: x[0]); c.PP = pp; c.PPh = [x for x in pp if x[1] >= 2]  # sorted: tuples() prunes by break
        c.jmax = int(c.lX/log(c.P[0]))
    def rho(c, v): return float((1-S((v-c.th)/(c.b-c.th)))*S(v/c.eta-1))
    def f(c, x):
        x = np.asarray(x, float); ins = np.abs(x) < c.th
        z = np.where(ins, 1-(x/c.th)**2, 1.0)
        return np.where(ins, np.exp(1-1/np.maximum(z, 1e-300)), 0.0)*(1+0.5j*x)
    def phi(c, u): u = np.asarray(u, float); return np.array([c.rho(v) for v in u/c.L])/np.maximum(u, 1e-300)
def tuples(items, keyf, r, budget):  # ordered r-tuples with sum of keyf < budget
    if r == 0: yield (); return
    for it in items:
        k = keyf(it)
        if k >= budget: break
        for rest in tuples(items, keyf, r-1, budget-k): yield (it,)+rest
def direct_M(c, s):
    tot = 0j
    def rec(start, lg, coef, mu):
        nonlocal tot
        tot += mu*c.f(lg/c.L)*coef*np.exp(-s*lg)
        for i in range(start, len(c.P)):
            q = c.P[i]; l2 = lg+log(q)
            if l2 >= c.lX: break
            rec(i+1, l2, coef*c.rho(log(q)/c.L), -mu)
    rec(0, 0.0, 1.0, 1); return complex(tot)
def F_free_primes(c, s, m, c0):  # sum over ordered prime m-tuples
    tot = 0j
    for tup in tuples(c.P, log, m, c.lX-c0):
        lg = c0+sum(log(q) for q in tup)
        tot += c.f(lg/c.L)*np.prod([c.rho(log(q)/c.L) for q in tup])*np.exp(-s*sum(log(q) for q in tup))
    return complex(tot)
def G_direct(c, s, r, c0):
    tot = 0j
    for tup in tuples(c.PP, lambda x: log(x[0]), r, c.lX-c0):
        lg = sum(log(x[0]) for x in tup)
        tot += c.f((c0+lg)/c.L)*np.prod([x[2] for x in tup])*np.exp(-s*lg)
    return complex(tot)
class ZeroG:
    def __init__(zg, c, s, Gam, du):
        t = s.imag; gam = ZEROS[np.abs(ZEROS-t) <= Gam]
        zg.nz = len(gam); zg.c = c; zg.du = du
        u = np.arange(0, c.b*c.L+du, du); ph = c.phi(u)
        Z = np.zeros(len(u), complex)
        for i in range(0, len(gam), 200):
            Z += np.exp(np.outer(u, (0.5+1j*gam[i:i+200]) - s)).sum(axis=1)
        g = -ph*Z; zg.pw = [np.array([1.0+0j])]
        for r in range(1, c.jmax+1): zg.pw.append(np.convolve(zg.pw[-1], g)*(du if r > 1 else 1.0))
        zg.pw[1] = g
    def G(zg, r, c0):
        c = zg.c
        if r == 0: return complex(c.f(c0/c.L))
        v = np.arange(len(zg.pw[r]))*zg.du
        return complex(np.sum(c.f((c0+v)/c.L)*zg.pw[r])*zg.du)
def expansion(c, s, Gfun, free_mode):
    tot = 0j
    for j in range(c.jmax+1):
        for lam in partitions(j):
            H = [h for h in lam if h >= 2]; m = lam.count(1); coef = (-1)**len(lam)/zlam(lam)
            for tup in tuples(c.P, log, len(H), c.lX):  # packet primes, ordered as in H
                c0 = sum(h*log(q) for h, q in zip(H, tup))
                if c0 + m*log(c.P[0]) >= c.lX: continue
                pk = np.prod([c.rho(log(q)/c.L)**h for h, q in zip(H, tup)])*np.exp(-s*c0)
                if free_mode == 'primes': Fm = F_free_primes(c, s, m, c0)
                else:
                    Fm = 0j
                    for jj in range(m+1):
                        for pt in tuples(c.PPh, lambda x: log(x[0]), jj, c.lX-c0):
                            lg = sum(log(x[0]) for x in pt)
                            Fm += math.comb(m, jj)*(-1)**jj*np.prod([x[2] for x in pt])*np.exp(-s*lg)*Gfun(m-jj, c0+lg)
                tot += coef*pk*Fm
    return complex(tot)
if __name__ == '__main__':
    cfgs = [(1e4, 0.56, 0.10, 0.66), (1.8e4, 0.60, 0.10, 0.70), (1.8e4, 0.52, 0.12, 0.60)]
    if len(sys.argv) > 1 and sys.argv[1] == 'tiny': cfgs = cfgs[:1]
    out = []
    for (T, th, eta, b) in cfgs:
        c = Case(T, th, eta, b)
        print(f"T={T:.0e} th={th} eta={eta} b={b}: |P|={len(c.P)} |PP|={len(c.PP)} pp-packets={[x[0] for x in c.PPh][:8]} jmax={c.jmax}", flush=True)
        for t in ([T*1.02] if len(sys.argv) > 1 else [T*1.02, T*1.37, 19650.3]):
            for sgn, y in [(1, 0.0), (1, 0.25), (-1, 0.25), (-1, 0.5)]:
                s = 0.5+sgn*y+1j*t
                Md = direct_M(c, s); A = expansion(c, s, None, 'primes')
                B = expansion(c, s, lambda r, c0: G_direct(c, s, r, c0), 'binom')
                res = dict(T=T, th=th, eta=eta, t=t, s=f"1/2{'+' if sgn>0 else '-'}{y}+it", M=abs(Md), errA=abs(A-Md), errB=abs(B-Md))
                for Gam in (60, 150):
                    zg = ZeroG(c, s, Gam, 0.004)
                    C = expansion(c, s, zg.G, 'binom'); res[f'errC_G{Gam}'] = abs(C-Md); res['nz'] = zg.nz
                out.append(res); print(json.dumps({k: (f"{v:.3e}" if isinstance(v, float) else v) for k, v in res.items()}), flush=True)
    json.dump(out, open('/Users/adriansager/playground/riemann/Riemann/Deep/UPAudit/logs/na_identity.json', 'w'), indent=1)
