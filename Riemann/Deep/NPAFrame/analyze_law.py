#!/usr/bin/env python3
"""npaframe: structure of a saved finite-N law (cg_law.py / cg_sets.py npz): types (number of doubles), conditional form factors,
violations of every length-3/4 word by level (triples, diagonal quartics, off-diagonal (2,2), corner (3,1)), and the
between-type part of Cov(|p_a|^2, |p_b|^2). usage: analyze_law.py file.npz [file2.npz ...]"""
import sys, math, itertools
import numpy as np

def zl(lam):
    z = 1
    for j in set(lam):
        a = lam.count(j); z *= j**a * math.factorial(a)
    return z

def parts(n, mx=None):
    mx = n if mx is None else mx
    if n == 0: yield (); return
    for k in range(min(n, mx), 0, -1):
        for r in parts(n-k, k): yield (k,) + r

for fn in sys.argv[1:]:
    d = np.load(fn, allow_pickle=True)
    w = d['sw']; M = d['sm']; TH = d['sth']; N = M.shape[1]
    w = w / w.sum()
    k = np.arange(0, N+1)
    P = np.einsum('ci,cik->ck', M, np.exp(1j*TH[:, :, None]*k[None, None, :]))   # p_0..p_N per configuration
    nd = (M == 2).sum(1); nt = (M == 3).sum(1); simple = (M == 1).sum(1)/N
    print(f"== {fn}: N={N} support={len(w)} minSimple={float(w @ simple):.6f}")
    for t in sorted(set(zip(nd, nt))):
        sel = (nd == t[0]) & (nt == t[1]); ws = w[sel].sum()
        S = (w[sel, None]*np.abs(P[sel, 1:])**2).sum(0)/ws
        print(f"  type doubles={t[0]} triples={t[1]}: weight {ws:.4f}, configs {sel.sum()}, cond. form factor/k - 1 = " +
              " ".join(f"{x:+.2f}" for x in S/np.arange(1, N+1) - 1))
    rows = []
    for n in range(2, N+1):
        PP = list(parts(n)); worst = {}
        for i, lam in enumerate(PP):
            for mu in PP[i:]:
                L = len(lam) + len(mu)
                if L < 3 or L > 4: continue
                val = (w * np.real(np.prod(P[:, list(lam)], 1) * np.conj(np.prod(P[:, list(mu)], 1)))).sum()
                tgt = zl(lam) if lam == mu else 0.0
                dev = (val - tgt)/math.sqrt(zl(lam)*zl(mu))
                cls = 'tri' if L == 3 else ('diag' if lam == mu else ('corner' if 3 in (len(lam), len(mu)) else 'offd'))
                if abs(dev) > abs(worst.get(cls, (0,))[0]): worst[cls] = (dev, lam, mu)
        rows.append((n, worst))
    for n, wd in rows:
        print(f"  level {n:2d} (Sum|xi| = {2*n/N:.3f}): " + "  ".join(f"{c} {v[0]:+.3f} {v[1]}|{v[2]}" for c, v in sorted(wd.items())))
    # between-type covariance of |p_a|^2, |p_b|^2 (normalised by ab), a <= b, a + b <= N
    types = sorted(set(zip(nd, nt))); tot = []
    for a in range(1, N+1):
        for b in range(a, N+1-a):
            Sa = []; Sb = []; ws = []
            for t in types:
                sel = (nd == t[0]) & (nt == t[1]); wt = w[sel].sum(); ws.append(wt)
                Sa.append((w[sel]*np.abs(P[sel, a])**2).sum()/wt); Sb.append((w[sel]*np.abs(P[sel, b])**2).sum()/wt)
            ws = np.array(ws); Sa = np.array(Sa); Sb = np.array(Sb)
            between = float((ws*(Sa - a)*(Sb - b)).sum())/(a*b)
            total = float((w*np.abs(P[:, a])**2*np.abs(P[:, b])**2).sum() - a*b*(2 if a == b else 1))/(a*b)
            tot.append((a, b, between, total))
    print("  Cov(|p_a|^2,|p_b|^2)/(ab): (a,b) between-type | total excess over Gaussian")
    print("   " + "  ".join(f"({a},{b}) {bt:+.3f}|{tt:+.3f}" for a, b, bt, tt in tot))
