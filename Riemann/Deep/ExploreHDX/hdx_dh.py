# hdx_dh.py (explore64): adversarial check. Davenport-Heilbronn's Weil-type form (no pole, atoms -b(n)/sqrt n, b from -L'/L) under the
# same localization: on Kronecker functions u(x) cos(Tx) with T ln p -> theta_p (all primes p < c, independent), the long-range atoms
# (n >= c', never inside a link) give  R(f_T)/|f_T|^2 -> -<u, P(theta) u>/|u|^2,  P(theta) = sum_{c'<=n<c} b(n) cos(k(n).theta)/sqrt n (T_n + T_n^*).
# Ceiling: Lam_long^DH = sup_theta lambda_max P(theta).  Random search over theta, then an arb Rayleigh quotient at rational theta, signed u.
import sys, math, numpy as np
import os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'ExploreFunctional', 'code'))
from xform import dh_b
from wmat import dh_coeffs
from flint import arb, fmpq, ctx
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh
ctx.prec = 128
rng = np.random.default_rng(64)

def fac(n):
    f = {}; p = 2
    while p*p <= n:
        while n % p == 0: f[p] = f.get(p, 0) + 1; n //= p
        p += 1
    if n > 1: f[n] = f.get(n, 0) + 1
    return f

def setup(c, cl, kind):
    b, _ = dh_b(int(c) + 2)
    if kind == 'zeta':
        from hdx_ceiling import ppow
        ns = [n for n in range(2, int(math.ceil(c))) if n < c and n >= cl and ppow(n)]
        wt = {n: math.log(ppow(n)[0])/math.sqrt(n) for n in ns}
    else:
        ns = [n for n in range(2, int(math.ceil(c))) if n < c and n >= cl and abs(b[n]) > 1e-14]
        wt = {n: b[n]/math.sqrt(n) for n in ns}
    ps = sorted(set(p for n in ns for p in fac(n)))
    return ns, wt, ps

def Pmat(c, ns, wt, th, ps, N):
    a = math.log(c)/2; h = 2*a/N; dg = {}
    for n in ns:
        beta = wt[n]*math.cos(sum(e*th[ps.index(p)] for p, e in fac(n).items())); s = math.log(n); k0 = int(math.floor(s/h))
        for k in (k0, k0 + 1):
            wgt = h - abs(k*h - s)
            if wgt > 0 and k < N: dg[k] = dg.get(k, np.zeros(N - k)) + beta*wgt/h
    if not dg: return None
    ks = sorted(dg); return diags([dg[k] for k in ks] + [dg[k] for k in ks], ks + [-k for k in ks], shape=(N, N), format='csr')

def top(A):
    w, V = eigsh(A, k=1, which='LA', tol=1e-9, maxiter=50000); return w[0], V[:, 0]

def certify(c, ns, th, ps, u, N, kind):
    a = arb(c).log()/2; h = 2*a/N; den = 10**6
    q = [fmpq(int(round(x/np.max(np.abs(u))*den)), den) for x in u]
    if kind == 'dh':
        A_, B_, _ = dh_coeffs(int(c) + 2); W = {n: B_[n]/arb(n).sqrt() for n in ns}
    else:
        from hdx_ceiling import ppow
        W = {n: arb(ppow(n)[0]).log()/arb(n).sqrt() for n in ns}
    thq = [fmpq(int(round(t*den)), den) for t in th]
    num = arb(0)
    for n in ns:
        beta = W[n]*sum((arb(e)*arb(thq[ps.index(p)]) for p, e in fac(n).items()), arb(0)).cos(); s = arb(n).log()
        k0 = int(math.floor(float(s.mid())/float(h.mid())))
        for k in (k0 - 1, k0, k0 + 1, k0 + 2):
            if k < 0 or k >= N: continue
            wgt = h - abs(arb(k)*h - s)
            if wgt.upper() <= 0: continue
            if wgt.lower() < 0: wgt = arb(0, wgt.upper())   # ball covering [0, upper] is safe either sign of the sum
            num += 2*beta*wgt*sum(arb(q[i])*arb(q[i + k]) for i in range(N - k))
    return num/(h*sum(arb(x)*arb(x) for x in q)), thq

if __name__ == '__main__':
    print('# kind, c, c_link, Lam_long = sup_theta lambda_max P(theta): float (N=3000), arb lower bound; zeta at theta = 0 for comparison')
    for c, cl in [(10, 10**0.5), (47, 47**0.5), (47, 31.98), (49, 47.0), (60, 47.0)]:
        for kind in ['zeta', 'dh']:
            ns, wt, ps = setup(c, cl, kind)
            if kind == 'zeta': cands = [np.zeros(len(ps))]
            else: cands = [rng.uniform(0, 2*math.pi, len(ps)) for _ in range(600)]
            best = (-1e9, None)
            for th in cands:
                A = Pmat(c, ns, wt, th, ps, 300)
                if A is None: continue
                v, _ = top(A)
                if v > best[0]: best = (v, th)
            th = best[1]; v, u = top(Pmat(c, ns, wt, th, ps, 3000))
            cv, thq = certify(c, ns, th, ps, u, 3000, kind)
            print(f'{kind:4s} {c:4d} {cl:7.3f}  {v:.4f} >= {float(cv.lower()):.4f}   atoms {ns[:12]}{"..." if len(ns) > 12 else ""}', flush=True)
