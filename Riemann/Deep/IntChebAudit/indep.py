#!/usr/bin/env python3
"""indep.py -- independent (non-arb, mpmath/numpy) recomputation for the integer-Chebyshev lower certificate.
Formula used (TeX thm:cert): A(x) = U_mu(x) - sum_j c_j log|Q_j(x)|, U_mu(x) = sum_k m_k psi_k(x),
psi_k(x) = [F(x-e_k) - F(x-e_{k+1})]/h_k, F(u) = u log|u| - u  (cell-by-cell form, not the beta-jump form of lowcert.c);
log|Q_j(x)| by direct high-precision Horner (not via roots) in the mpmath pass.
G_i (TeX form, normalised root measure) = L_i (M - B + c_i d_i) + mean_alpha U_mu(alpha) - sum_{j!=i} c_j mean_alpha log|Q_j(alpha)|
with roots from mpmath.polyroots and U_mu at complex alpha as Re sum_k m_k [F(alpha-e_k) - F(alpha-e_{k+1})]/h_k (principal log),
NOT via resultants (gcert.c uses resultants).
usage: indep.py CERT.json MODE [N|FILE]   MODE in exact, rand, scan, escan, G, leaves
"""
import json, sys, time, random, math
from fractions import Fraction as Fr
import mpmath as mp

path, mode = sys.argv[1], sys.argv[2]
N = int(sys.argv[3]) if len(sys.argv) > 3 else 0
C = json.load(open(path))
E = C['cell_edges']; Ms = C['masses']; lib = [[int(a) for a in q] for q in C['library']]; cs = C['multipliers']; V = C['V']
K = len(Ms)
act = [(q, c) for q, c in zip(lib, cs) if c > 0]
Mx = sum(Fr(x) for x in Ms); Bx = sum(Fr(c) * (len(q) - 1) for c, q in zip(cs, lib))
print('backend', mp.libmp.BACKEND, 'K', K, 'J', len(lib), 'active', len(act), 'V', V, flush=True)

_pre = {}
def _prep():
    d = mp.mp.dps
    if _pre.get('dps') != d:
        _pre['dps'] = d
        _pre['E'] = [mp.mpf(e) for e in E]
        _pre['rho'] = [mp.mpf(Ms[k]) / (mp.mpf(E[k + 1]) - mp.mpf(E[k])) for k in range(K)]
        with mp.workdps(d + 50):
            _pre['Q'] = [([mp.mpf(a) for a in q], mp.mpf(c)) for q, c in act]
def A_mp(x):
    # x: mpf. U_mu by cells (density rho_k = m_k/h_k); log|Q| by Horner at high precision.
    _prep(); Em = _pre['E']; rho = _pre['rho']; one = mp.mpf(1)
    Fv = []
    for e in Em:
        u = x - e
        Fv.append(mp.mpf(0) if u == 0 else u * (mp.log(abs(u)) - one))
    U = mp.fsum(rho[k] * (Fv[k] - Fv[k + 1]) for k in range(K))
    S = mp.mpf(0)
    with mp.workdps(mp.mp.dps + 50):
        for qm, c in _pre['Q']:
            v = mp.polyval(qm, x)
            if v == 0: return mp.inf
            S += c * mp.log(abs(v))
    return U - S

if mode == 'exact':
    # independent exact checks: support, positivity, M, B, symmetry (own reflection via Horner composition), irreducibility (sympy)
    import sympy as sp
    X = sp.Symbol('x')
    Ex = [Fr(e) for e in E]
    print('edges strictly increasing', all(Ex[i] < Ex[i + 1] for i in range(K)), 'e0', Ex[0], 'eK', Ex[-1])
    print('min mass', min(Fr(m) for m in Ms) >= 0, 'min c', min(Fr(c) for c in cs) >= 0)
    print('M-1 =', float(Mx - 1), ' B =', float(Bx), ' B-M =', float(Bx - Mx), ' B<M', Bx < Mx)
    sym_e = all(Ex[k] + Ex[K - k] == 1 for k in range(K + 1)); sym_m = all(Fr(Ms[k]) == Fr(Ms[K - 1 - k]) for k in range(K))
    def refl(q):
        P = sp.Poly(q, X); R = sp.Poly(P.as_expr().subs(X, 1 - X), X)
        co = [int(a) for a in R.all_coeffs()]
        return tuple(-a for a in co) if co[0] < 0 else tuple(co)
    idx = {tuple(q): i for i, q in enumerate(lib)}
    sym_q = all(refl(q) in idx and Fr(cs[idx[refl(q)]]) == Fr(cs[i]) for i, q in enumerate(lib))
    print('symmetry: edges', sym_e, 'masses', sym_m, 'library+multipliers', sym_q)
    t0 = time.time(); bad = []
    for q in lib:
        P = sp.Poly(q, X)
        cont, fl = sp.factor_list(P.as_expr())
        ok = (abs(cont) == 1 and len(fl) == 1 and fl[0][1] == 1 and q[0] >= 1)
        if not ok: bad.append(q)
    print('sympy irreducible+primitive+lc>=1: bad =', len(bad), '(%.0fs)' % (time.time() - t0))
    print('distinct', len({tuple(q) for q in lib}) == len(lib))
    mid = -sum(mp.mpf(c) * mp.log(q[0]) for q, c in zip(lib, cs))
    print('-sum c log a =', mp.nstr(mid, 12), ' /M =', mp.nstr(mid / mp.mpf(Mx.numerator) * Mx.denominator, 12))
    mp.mp.dps = 40
    Mm = mp.mpf(Mx.numerator) / Mx.denominator
    print('V as exact double =', mp.nstr(mp.mpf(Fr(V).numerator) / Fr(V).denominator, 25))
    print('exp(V/M) =', mp.nstr(mp.exp(mp.mpf(Fr(V).numerator) / Fr(V).denominator / Mm), 25), ' exp(V) =', mp.nstr(mp.exp(mp.mpf(V)), 25))
    print('largest active root count (sum deg active) =', sum(len(q) - 1 for q, c in act))

if mode == 'rand':
    mp.mp.dps = 30
    random.seed(20261119)
    lo = mp.inf; arg = None; t0 = time.time(); below = 0
    for n in range(N):
        x = mp.mpf(random.random())
        a = A_mp(x)
        if a < V: below += 1
        if a < lo: lo, arg = a, x
        if n == 99: print('100 pts in %.1fs' % (time.time() - t0), flush=True)
    print('rand N=%d: min A = %s at x=%s ; points with A < V: %d ; (%.0fs)' % (N, mp.nstr(lo, 15), mp.nstr(arg, 12), below, time.time() - t0))
    # symmetry spot check A(x) = A(1-x)
    mx = 0
    for n in range(50):
        x = mp.mpf(random.random()) / 2
        mx = max(mx, abs(A_mp(x) - A_mp(1 - x)))
    print('max |A(x)-A(1-x)| over 50 pts =', mp.nstr(mx, 5))

def roots_mp(q, dps=40):
    with mp.workdps(dps):
        rts = mp.polyroots([mp.mpf(a) for a in q], maxsteps=800, extraprec=400)
    out = []
    for r in rts:
        r = mp.mpc(r)
        out.append(mp.mpf(r.real) if abs(r.imag) < mp.mpf(10) ** (-dps + 5) else r)
    return out

if mode == 'G':
    mp.mp.dps = 40
    Em = [mp.mpf(e) for e in E]; rho = [mp.mpf(Ms[k]) / (Em[k + 1] - Em[k]) for k in range(K)]
    Mm = mp.mpf(Mx.numerator) / Mx.denominator; Bm = mp.mpf(Bx.numerator) / Bx.denominator
    def Fc(u):
        if u == 0: return mp.mpf(0)
        if isinstance(u, mp.mpf): return u * (mp.log(abs(u)) - 1)      # real: u log|u| - u
        return mp.re(u * mp.log(u) - u)                                  # principal branch
    def U(z):
        Fv = [Fc(z - e) for e in Em]
        return mp.fsum(rho[k] * (Fv[k] - Fv[k + 1]) for k in range(K))
    actidx = [j for j, c in enumerate(cs) if c > 0]
    t0 = time.time(); res = []
    for i, q in enumerate(lib):
        d = len(q) - 1; rts = roots_mp(q)
        # residual check of the roots
        with mp.workdps(80):
            resid = max(abs(mp.polyval([mp.mpf(a) for a in q], r)) for r in rts)
        meanU = mp.fsum(U(r) for r in rts) / d
        cross = mp.mpf(0)
        with mp.workdps(90):
            for j in actidx:
                if j == i: continue
                qj = [mp.mpf(a) for a in lib[j]]
                cross += mp.mpf(cs[j]) * mp.fsum(mp.log(abs(mp.polyval(qj, r))) for r in rts) / d
        Li = mp.log(q[0]) / d
        G = Li * (Mm - Bm + mp.mpf(cs[i]) * d) + meanU - cross
        res.append((G, i)); print('G %d %s resid %s' % (i, mp.nstr(G, 17), mp.nstr(resid, 3)), flush=True)
    G, i = min(res)
    print('MIN G = %s at i=%d (%.0fs)' % (mp.nstr(G, 15), i, time.time() - t0))

if mode == 'scan':
    import numpy as np
    mp.mp.dps = 30
    rts = []
    for q, c in act:
        rr = roots_mp(q)
        rts.append((c, math.log(q[0]), np.array([complex(r) for r in rr])))
    Ea = np.array(E, dtype=float); ha = np.diff(Ea); rho = np.array(Ms) / ha
    X1 = 1.0
    xs = np.linspace(0.0, X1, N + 1)
    Avals = np.empty_like(xs)
    for s in range(0, len(xs), 500):
        x = xs[s:s + 500]
        D = x[:, None] - Ea[None, :]
        with np.errstate(divide='ignore', invalid='ignore'):
            Fv = np.where(D == 0, 0.0, D * (np.log(np.abs(D)) - 1))
        U = ((Fv[:, :-1] - Fv[:, 1:]) * rho[None, :]).sum(axis=1)
        S = np.zeros_like(x)
        for c, la, rr in rts:
            with np.errstate(divide='ignore'):
                S += c * (la + np.log(np.abs(x[:, None] - rr[None, :])).sum(axis=1))
        Avals[s:s + 500] = U - S
    i0 = int(np.argmin(Avals))
    print('grid N=%d on [0,%g]: min A(float64) = %.12f at x=%.9f' % (N, X1, Avals[i0], xs[i0]), flush=True)
    loc = [i for i in range(1, len(xs) - 1) if Avals[i] <= Avals[i - 1] and Avals[i] <= Avals[i + 1]]
    loc.sort(key=lambda i: Avals[i])
    gr = (mp.sqrt(5) - 1) / 2; best = []
    for i in loc[:12]:
        a, b = mp.mpf(xs[i - 1]), mp.mpf(xs[i + 1])
        c1 = b - gr * (b - a); c2 = a + gr * (b - a); f1 = A_mp(c1); f2 = A_mp(c2)
        for it in range(45):
            if f1 < f2: b, c2, f2 = c2, c1, f1; c1 = b - gr * (b - a); f1 = A_mp(c1)
            else: a, c1, f1 = c1, c2, f2; c2 = a + gr * (b - a); f2 = A_mp(c2)
        xm = (a + b) / 2; fm = A_mp(xm)
        best.append((fm, xm)); print('local min x=%s A=%s (grid %.12f)' % (mp.nstr(xm, 15), mp.nstr(fm, 15), Avals[i]), flush=True)
    fm, xm = min(best)
    print('REFINED MIN A = %s at x=%s ; V = %r ; margin = %s' % (mp.nstr(fm, 15), mp.nstr(xm, 15), V, mp.nstr(fm - V, 5)))

if mode == 'escan':
    # edge-adapted scan: uniform grid + log-spaced offsets on both sides of every edge + 16 Chebyshev points per cell
    import numpy as np
    mp.mp.dps = 30
    rts = []
    for q, c in act:
        rr = roots_mp(q)
        rts.append((c, math.log(q[0]), np.array([complex(r) for r in rr])))
    Ea = np.array(E, dtype=float); ha = np.diff(Ea); rho = np.array(Ms) / ha
    offs = 10.0 ** np.linspace(-13, -3.5, 60)
    pts = [np.linspace(0.0, 1.0, N + 1)]
    for e in Ea:
        pts.append(e + offs); pts.append(e - offs)
    th = np.cos(np.pi * (np.arange(16) + 0.5) / 16)
    for k in range(K):
        pts.append(0.5 * (Ea[k] + Ea[k + 1]) + 0.5 * ha[k] * th)
    xs = np.unique(np.concatenate(pts)); xs = xs[(xs >= 0) & (xs <= 1)]
    Avals = np.empty_like(xs)
    def Anp(x):
        D = x[:, None] - Ea[None, :]
        with np.errstate(divide='ignore', invalid='ignore'):
            Fv = np.where(D == 0, 0.0, D * (np.log(np.abs(D)) - 1))
        U = ((Fv[:, :-1] - Fv[:, 1:]) * rho[None, :]).sum(axis=1)
        S = np.zeros_like(x)
        for c, la, rr in rts:
            with np.errstate(divide='ignore'):
                S += c * (la + np.log(np.abs(x[:, None] - rr[None, :])).sum(axis=1))
        return U - S
    for s in range(0, len(xs), 500):
        Avals[s:s + 500] = Anp(xs[s:s + 500])
    i0 = int(np.nanargmin(Avals))
    print('escan: %d points; min A(float64) = %.12f at x=%.15f' % (len(xs), Avals[i0], xs[i0]), flush=True)
    loc = [i for i in range(1, len(xs) - 1) if Avals[i] <= Avals[i - 1] and Avals[i] <= Avals[i + 1]]
    loc.sort(key=lambda i: Avals[i])
    gr = (mp.sqrt(5) - 1) / 2; best = []
    for i in loc[:16]:
        a, b = mp.mpf(xs[i - 1]), mp.mpf(xs[i + 1])
        c1 = b - gr * (b - a); c2 = a + gr * (b - a); f1 = A_mp(c1); f2 = A_mp(c2)
        for it in range(50):
            if f1 < f2: b, c2, f2 = c2, c1, f1; c1 = b - gr * (b - a); f1 = A_mp(c1)
            else: a, c1, f1 = c1, c2, f2; c2 = a + gr * (b - a); f2 = A_mp(c2)
        xm = (a + b) / 2; fm = A_mp(xm)
        kk = int(np.argmin(np.abs(Ea - float(xm))))
        best.append((fm, xm)); print('local min x=%s A=%s (float64 %.12f) nearest edge %d at dist %.3e' % (mp.nstr(xm, 17), mp.nstr(fm, 15), Avals[i], kk, abs(Ea[kk] - float(xm))), flush=True)
    fm, xm = min(best)
    print('REFINED MIN A = %s at x=%s ; V = %r ; margin = %s' % (mp.nstr(fm, 15), mp.nstr(xm, 17), V, mp.nstr(fm - V, 5)))

if mode == 'leaves':
    # evaluate A (mpmath, dps 30) on 41 points of each leaf listed as 'TIGHT lo hi lb depth' in file sys.argv[3]
    mp.mp.dps = 30
    for line in open(sys.argv[3]):
        if not line.startswith('TIGHT'): continue
        _, lo, hi, lb, dep = line.split()
        lo, hi = mp.mpf(lo), mp.mpf(hi)
        vals = [(A_mp(lo + (hi - lo) * t / 40), t) for t in range(41)]
        vmin, tmin = min(vals)
        print('leaf [%s, %s] lb=%s : min A on 41 pts = %s (t=%d/40), A(lo)=%s A(hi)=%s' % (mp.nstr(lo, 12), mp.nstr(hi, 12), lb, mp.nstr(vmin, 13), tmin, mp.nstr(vals[0][0], 13), mp.nstr(vals[-1][0], 13)), flush=True)

if mode == 'gscan':
    # all local minima of A on a cell-adapted grid, every one refined by vectorised float64 golden section, lowest 8 by mpmath
    import numpy as np
    mp.mp.dps = 30
    rts = []
    for q, c in act:
        rts.append((c, math.log(q[0]), np.array([complex(r) for r in roots_mp(q)])))
    Ea = np.array(E, dtype=float); ha = np.diff(Ea); rho = np.array(Ms) / ha
    def Anp(x):
        out = np.empty_like(x)
        for s in range(0, len(x), 400):
            xx = x[s:s + 400]
            D = xx[:, None] - Ea[None, :]
            with np.errstate(divide='ignore', invalid='ignore'):
                Fv = np.where(D == 0, 0.0, D * (np.log(np.abs(D)) - 1))
            U = ((Fv[:, :-1] - Fv[:, 1:]) * rho[None, :]).sum(axis=1)
            S = np.zeros_like(xx)
            for c, la, rr in rts:
                with np.errstate(divide='ignore'):
                    S += c * (la + np.log(np.abs(xx[:, None] - rr[None, :])).sum(axis=1))
            out[s:s + 400] = U - S
        return out
    offs = 10.0 ** np.linspace(-13, -3.5, 60)
    pts = [np.linspace(0.0, 1.0, 100001)]
    for e in Ea: pts.append(e + offs); pts.append(e - offs)
    for k in range(K): pts.append(Ea[k] + ha[k] * (np.arange(1, N) / N))
    xs = np.unique(np.concatenate(pts)); xs = xs[(xs >= 0) & (xs <= 1)]
    t0 = time.time(); Av = Anp(xs)
    loc = np.where((Av[1:-1] <= Av[:-2]) & (Av[1:-1] <= Av[2:]))[0] + 1
    print('gscan: %d points (%d per cell), %d local minima, sampled min %.12f (%.0fs)' % (len(xs), N, len(loc), Av.min(), time.time() - t0), flush=True)
    a = xs[loc - 1].copy(); b = xs[loc + 1].copy(); g = (math.sqrt(5) - 1) / 2
    c1 = b - g * (b - a); c2 = a + g * (b - a); f1 = Anp(c1); f2 = Anp(c2)
    for it in range(55):
        L = f1 < f2
        b = np.where(L, c2, b); a = np.where(L, a, c1)
        n1 = np.where(L, b - g * (b - a), c2); n2 = np.where(L, c1, a + g * (b - a))
        fn = Anp(np.where(L, n1, n2))
        # L: new c2 = old c1 (f2 = old f1), new c1 = n1 (f1 = fn); ~L: new c1 = old c2 (f1 = old f2), new c2 = n2 (f2 = fn)
        f1n = np.where(L, fn, f2); f2n = np.where(L, f1, fn)
        f1, f2, c1, c2 = f1n, f2n, n1, n2
    xm = 0.5 * (a + b); fm = Anp(xm)
    order = np.argsort(fm)
    print('float64-refined: global min %.12f at x=%.15f ; minima below V+1e-7: %d ; below V: %d (%.0fs)' % (fm[order[0]], xm[order[0]], int((fm < V + 1e-7).sum()), int((fm < V).sum()), time.time() - t0), flush=True)
    gr = (mp.sqrt(5) - 1) / 2; best = []
    for i in order[:8]:
        aa, bb = mp.mpf(xm[i]) - mp.mpf('2e-7'), mp.mpf(xm[i]) + mp.mpf('2e-7')
        d1 = bb - gr * (bb - aa); d2 = aa + gr * (bb - aa); g1 = A_mp(d1); g2 = A_mp(d2)
        for it in range(45):
            if g1 < g2: bb, d2, g2 = d2, d1, g1; d1 = bb - gr * (bb - aa); g1 = A_mp(d1)
            else: aa, d1, g1 = d1, d2, g2; d2 = aa + gr * (bb - aa); g2 = A_mp(d2)
        x_ = (aa + bb) / 2; f_ = A_mp(x_); best.append((f_, x_))
        print('mp-refined x=%s A=%s (float64 %.12f)' % (mp.nstr(x_, 15), mp.nstr(f_, 15), fm[i]), flush=True)
    f_, x_ = min(best)
    print('GLOBAL MIN (independent) A = %s at x=%s ; A-V = %s' % (mp.nstr(f_, 15), mp.nstr(x_, 15), mp.nstr(f_ - V, 6)))
