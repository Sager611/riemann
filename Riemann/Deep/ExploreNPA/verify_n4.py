#!/usr/bin/env python3
"""explore52 — rigorous checks at N = 4 on the lattice Z_M, for ALL off-line depths.
(a) primal (upper bound): the CG optimum for a set (e.g. P3, off world, depth grid) re-solved on its support at 50 digits (mpmath):
    every Haar row and the normalisation hold to < 1e-40 and all weights are >= 0  =>  min E[#simple on-line]/N <= value.
(b) dual (lower bound, every law, every depth s in [0, oo)): y = LP duals; g(C) = ns/N - y_norm - sum_r y_r f_r(C) must be >= -eps for
    every pure on-line configuration (float scan, mpmath recheck of the 50 smallest slacks), every configuration
    'one pair + two on-line points' (h(u) = u^8 g, u = e^{-s} in [0, 1]: arb bisection) and every 'two pairs' configuration
    (h(u1, u2) = u1^8 u2^8 g on [0, 1]^2: 2-D arb bisection).  Certified bound = dual objective - max eps.
usage: verify_n4.py M depths SET_primal SET_dual
"""
import sys, itertools, math, time
import numpy as np, mpmath as mp
from flint import arb
import npa_cg as cg

mp.mp.dps = 50
N = 4

def run(M, depths, sets, world):
    spec = cg.rowspec(N, sets); b = np.array([t for (_, _, t, _) in spec])
    cg.RAYS = cg.rays(N, spec) if ('-rays' in sys.argv and world == 'off') else None
    g = cg.blocks(N, M, depths, 'on'); P0, ns0, lab0 = next(g)
    A = cg.feats(P0, spec); c = ns0 / N; labs = [lab0(i) for i in range(len(ns0))]
    while True:
        val, x, y, res = cg.solve(A, b, c)
        best = []
        for P, ns, lab in cg.blocks(N, M, depths, world):
            F = cg.feats(P, spec)
            cs = np.maximum(1.0, np.abs(np.vstack([F, np.ones((1, F.shape[1]))])).max(axis=0))
            rc = (ns / N - y[:-1] @ F - y[-1]) / cs
            idx = np.where(rc < -1e-10)[0]
            best += [(rc[i], F[:, i], ns[i] / N, lab(i)) for i in idx[np.argsort(rc[idx])[:200]]]
        if not best: break
        best.sort(key=lambda z: z[0]); best = best[:400]
        A = np.hstack([A, np.array([z[1] for z in best]).T]); c = np.concatenate([c, [z[2] for z in best]]); labs += [z[3] for z in best]
    return val, x, y, labs, spec, b

def mp_p(lab, M):
    on, prs = lab
    p = [mp.mpc(0)] * (N + 1)
    for k in range(N + 1):
        z = mp.mpc(0)
        for xx in on: z += mp.expjpi(mp.mpf(2 * k * xx) / M)
        for (xx, s) in prs: z += 2 * mp.cosh(k * mp.mpf(s)) * mp.expjpi(mp.mpf(2 * k * xx) / M)
        p[k] = z
    return p

def mp_feats(lab, M, spec):
    p = mp_p(lab, M); out = []
    for (lam, mu, t, s) in spec:
        a = mp.mpc(1); bb = mp.mpc(1)
        for k in lam: a *= p[k]
        for k in mu: bb *= p[k]
        out.append(mp.re(a * mp.conj(bb)) / mp.sqrt(mp.mpf(s * s)))
    return out

def simple_on(lab):
    on = lab[0]; return sum(1 for v in set(on) if on.count(v) == 1)

def primal_check(M, depths, sets):
    val, x, y, labs, spec, b = run(M, depths, sets, 'off')
    supp = [i for i in range(len(x)) if x[i] > 0]
    B = mp.matrix(len(spec) + 1, len(supp))
    for j, i in enumerate(supp):
        f = mp_feats(labs[i], M, spec)
        for r in range(len(spec)): B[r, j] = f[r]
        B[len(spec), j] = 1
    rhs = mp.matrix([mp.mpf(t) for t in b] + [1])
    w = mp.lu_solve(B, rhs) if len(supp) == len(spec) + 1 else mp.qr_solve(B, rhs)[0]
    res = mp.norm(B * w - rhs, mp.inf)
    obj = sum(w[j] * simple_on(labs[i]) for j, i in enumerate(supp)) / N
    print(f"[primal {sets}] M={M} float={val:.6f} support={len(supp)} mp-value={mp.nstr(obj, 12)} resid={mp.nstr(res, 3)} "
          f"min w={mp.nstr(min(w), 3)} off-line pair cols={sum(1 for i in supp if labs[i][1])}")
    for j, i in enumerate(supp):
        if labs[i][1]: print(f"     w={mp.nstr(w[j], 4)} {labs[i]}")
    return obj

def laurent_mul(a, b):  # dict power -> mpc
    out = {}
    for i, u in a.items():
        for j, v in b.items(): out[i + j] = out.get(i + j, 0) + u * v
    return out

def dual_check(M, depths, sets):
    val, x, y, labs, spec, b = run(M, depths, sets, 'off')
    yr = [mp.mpf(float(v)) for v in y[:-1]]; y0 = mp.mpf(float(y[-1]))
    dual = sum(yr[r] * mp.mpf(b[r]) for r in range(len(spec))) + y0
    t0 = time.time(); worst = mp.mpf(0)
    # pure on-line (float scan + mp recheck)
    from offline_lp import online_configs
    X = online_configs(N, M); slack = []
    for P, ns, lab in cg.blocks(N, M, depths, 'on'):
        F = cg.feats(P, spec); gv = ns / N - y[:-1] @ F - y[-1]
        for i in np.argsort(gv)[:50]:
            f = mp_feats(lab(i), M, spec); gg = mp.mpf(simple_on(lab(i))) / N - y0 - sum(yr[r] * f[r] for r in range(len(spec)))
            worst = min(worst, gg)
    print(f"[dual {sets}] float dual objective={val:.6f} mp={mp.nstr(dual, 10)}; on-line min slack={mp.nstr(worst, 3)} t={time.time()-t0:.0f}s", flush=True)
    # one pair at 0 (depth s) + on-line points x1 <= x2, s in [0, oo):  t = e^s, u = 1/t
    w1 = 0.0; nbox = 0; lead = mp.mpf(10 ** 9)
    for x1 in range(M):
        for x2 in range(x1, M):
            ck = [mp.expjpi(mp.mpf(2 * k * x1) / M) + mp.expjpi(mp.mpf(2 * k * x2) / M) for k in range(N + 1)]
            pk = [{k: mp.mpc(1), -k: mp.mpc(1)} if k else {0: mp.mpc(2)} for k in range(N + 1)]
            for k in range(N + 1): pk[k] = dict(pk[k]); pk[k][0] = pk[k].get(0, 0) + ck[k]
            pc = [{e: mp.conj(v) for e, v in d.items()} for d in pk]
            gl = {0: mp.mpf(2 if x1 != x2 else 0) / N - y0}
            for r, (lam, mu, t, s) in enumerate(spec):
                a = {0: mp.mpc(1)}
                for k in lam: a = laurent_mul(a, pk[k])
                for k in mu: a = laurent_mul(a, pc[k])
                for e, v in a.items(): gl[e] = gl.get(e, 0) - yr[r] * mp.re(v) / s
            coef = [gl.get(8 - j, mp.mpf(0)) for j in range(17)]   # h(u) = sum coef[j] u^j
            lead = min(lead, coef[0]) if 'lead' in dir() else coef[0]
            # uniform in g: on a box [a, b] of u = e^{-s}, g = h/u^8 >= h/a^8; accept when lo(h) >= -TOL * a^8 (s in [0, SMAX])
            m, nb = bisect1([arb(mp.nstr(c, 40)) for c in coef], math.exp(-SMAX), 1.0, 0); nbox += nb
            w1 = min(w1, m)
    print(f"  one pair + 2 on-line, depths s in [0, {SMAX}]: certified g >= {w1:.3e} (tol {TOL}; {nbox} boxes); "
          f"min leading coefficient (s -> oo, g ~ coef * e^(8s)) = {mp.nstr(lead, 6)} t={time.time()-t0:.0f}s", flush=True)
    # two pairs at angles 0 and x, depths (s1, s2): float scan on [0, 15]^2 (301 x 301 per angle; not interval-certified)
    w2 = 1.0; arg = None
    sg = np.linspace(0, 15, 301); S1, S2 = np.meshgrid(sg, sg, indexing='ij'); kk = np.arange(N + 1)
    yv = np.array([float(v) for v in yr]); y0f = float(y0)
    for xx in range(M // 2 + 1):
        ph = np.exp(2j * np.pi * kk * xx / M)
        Pk = 2 * np.cosh(kk[None, None, :] * S1[:, :, None]) + 2 * np.cosh(kk[None, None, :] * S2[:, :, None]) * ph[None, None, :]
        gv = -y0f * np.ones(S1.shape)
        for r, (lam, mu, t, s) in enumerate(spec):
            a = np.prod([Pk[:, :, k] for k in lam], axis=0) * np.conj(np.prod([Pk[:, :, k] for k in mu], axis=0))
            gv -= yv[r] * a.real / s
        # relative slack: g / (1 + |largest term|) guards against float cancellation at large depth
        scale = 1 + np.abs(Pk[:, :, N]) ** 2
        i = np.unravel_index(np.argmin(gv / scale), gv.shape)
        if gv[i] / scale[i] < w2: w2 = gv[i] / scale[i]; arg = (xx, sg[i[0]], sg[i[1]], gv[i])
    print(f"  two pairs, depths in [0,15]^2 (float scan): min g/(1+|p_N|^2) = {w2:.3e} at angle,s1,s2,g = {arg} t={time.time()-t0:.0f}s", flush=True)
    w2 = min(0.0, arg[3])
    cert = dual + min(worst, mp.mpf(w1), 0)
    print(f"BOUND (lattice Z_{M}; on-line + one-pair parts arb-certified for every depth, two-pair part float-scanned): "
          f"min E[#simple on-line]/N >= {mp.nstr(cert, 10)} (two-pair float min g = {w2:.2e})", flush=True)
    return dual, worst, w1, y, spec, M

def mul2(a, b):
    out = {}
    for (i1, i2), u in a.items():
        for (j1, j2), v in b.items():
            kk = (i1 + j1, i2 + j2); out[kk] = out.get(kk, 0) + u * v
    return out

def eval2(Ca, U1, U2):
    v = arb(0)
    for row in reversed(Ca):
        w = arb(0)
        for c in reversed(row): w = w * U2 + c
        v = v * U1 + w
    return v

def bisect2(Ca, a1, b1, a2, b2, depth):
    v = eval2(Ca, arb((a1 + b1) / 2, (b1 - a1) / 2), arb((a2 + b2) / 2, (b2 - a2) / 2))
    lo = float(v.mid()) - float(v.rad())
    if lo >= 0 or depth >= 18: return min(lo, 0.0), 1
    m1, m2 = (a1 + b1) / 2, (a2 + b2) / 2; best = 0.0; nb = 1
    for (p, q, r_, s_) in [(a1, m1, a2, m2), (m1, b1, a2, m2), (a1, m1, m2, b2), (m1, b1, m2, b2)]:
        m, n = bisect2(Ca, p, q, r_, s_, depth + 1); best = min(best, m); nb += n
    return best, nb

def horner(cs, U):
    v = arb(0)
    for c in reversed(cs): v = v * U + c
    return v

TOL = 1e-9; SMAX = 12.0
def bisect1(cs, a, b, depth):
    """returns (certified lower bound of g = h/u^8 on [a, b], boxes); unresolved boxes report their g-lower bound."""
    v = horner(cs, arb((a + b) / 2, (b - a) / 2)); lo = float(v.mid()) - float(v.rad())
    glo = lo / a ** 8 if lo < 0 else lo / b ** 8
    if glo >= -TOL or depth >= 40: return min(glo, 0.0), 1
    m1, n1 = bisect1(cs, a, (a + b) / 2, depth + 1); m2, n2 = bisect1(cs, (a + b) / 2, b, depth + 1)
    return min(m1, m2), n1 + n2 + 1

if __name__ == '__main__':
    M = int(sys.argv[1]); depths = [float(t) for t in sys.argv[2].split(',')]
    if sys.argv[3] != '-': primal_check(M, depths, sys.argv[3])
    if len(sys.argv) > 4: dual_check(M, depths, sys.argv[4])
