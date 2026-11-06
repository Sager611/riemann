#!/usr/bin/env python3
# check.py -- independent floating-point cross-check of Ingham4's cert.py bounds (NOT the certificate; numpy, log-space trapezoid rules).
# Adapted from Riemann/Deep/Ingham3/check.py: recomputes Z = 2[x^{s0-1} S + log x int_{s0}^{B} x^{s-1} F(s) ds] directly, with
#   F(s) = N(s,T0) G(delta T0) + int_{T0}^oo G(delta t) dN(s,t), N = KLN with its exact (log kt)^{2s}(log t)^{5-4s} t^{8(1-s)/3} shape
#   (best admissible row at each s), or min(KLN, CS) with CS's exact B1 t^{3(1-s)/(2-s)} (log t)^{(7-5s)/(2-s)} + B2 log^2 t + B3 log t shape,
#   S = 14 log 14 G(14 delta) + int_14^oo G(delta t)(1 + log t) dt, G = the kernel envelope (trapezoid 1_[0,a]*1_[0,1-a], a = 1/10, pipeline A1;
#   PIPE=A0: the triangle a = 1/2 and KLN at H0 = 3.0610046e10); plus PP and 2x^{-3} as in cert.py; and a random test of the kernel bound.
# usage: uv run --no-project --with python-flint --with mpmath --with sympy --with numpy python3 Riemann/Deep/Ingham4/check.py   (reads logs/cert_result.json)
import sys, os, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import cert
from flint import arb
PIPE = os.environ.get('PIPE', 'A1')
cert.main_setup(('trap',) if PIPE == 'A1' else ('tri',))
S0 = cert.build(False, pipe=PIPE); S1 = cert.build(True, pipe=PIPE)
A_K = float(cert.A_TRAP) if PIPE == 'A1' else 0.5      # kernel: 1_[0,a] * 1_[0,1-a] (a = 1/2: the triangle)
B = 0.8733819666132076; LT0 = math.log(3e12); S1c = 0.2173
ROWS = [(float(P.c0.mid()), float(P.c1.mid()), float(P.smin.mid()), float(P.k.mid()), float(P.q0.mid()), float(P.q1.mid())) for P in S0.P]
CSR = [(float(r.s1.mid()), float(r.s2.mid()), float(r.B1.mid()), float(arb(v[3]).mid()), float(arb(v[4]).mid())) for r, v in zip([cert.CSRow(v) for v in cert.CS_ROWS], cert.CS_ROWS)]
def best_row(s):
    return min((r[0] + r[1]*s, r) for r in ROWS if s > r[2])[1]
def cs_row(s):
    return max((r for r in CSR if r[0] <= s <= r[1]), key=lambda r: r[2])
def Senv(y):
    with np.errstate(divide='ignore', invalid='ignore'):
        sinc = np.where(y > 0, np.sin(y)/np.where(y > 0, y, 1), 1.0)
    return np.where(y <= math.pi, np.maximum(sinc, S1c), np.minimum(S1c, 1/np.maximum(y, 1e-300)))
def logG(lam, dmax, a=None):
    a = A_K if a is None else a; b = 1 - a
    il = 1/np.maximum(lam, 1e-300)
    return np.log(Senv(a*lam/2) + 2*dmax*np.minimum(a, il)) + np.log(Senv(b*lam/2) + 2*dmax*np.minimum(b, il))
def kernel_random_test(a, ntest=20000, seed=1):
    """|Phi~(z)|/int Phi <= G(Im z) for Re z in [0, dmax]; Phi = 1_[0,a] * 1_[0,1-a], Phi~(z) = ((e^{az}-1)/z)((e^{(1-a)z}-1)/z)."""
    rng = np.random.default_rng(seed); worst = 0.0; b = 1 - a
    for dmax in (1e-9, 0.025, 0.15):
        lam = np.exp(rng.uniform(math.log(1e-3), math.log(1e5), ntest)); re = rng.uniform(0, dmax, ntest)
        z = re + 1j*lam
        ratio = np.abs((np.exp(a*z) - 1)/z*(np.exp(b*z) - 1)/z)/(a*b)
        worst = max(worst, float(np.max(ratio/np.exp(logG(lam, dmax, a)))))
    return worst
def lse(a): m = np.max(a); return m + math.log(np.sum(np.exp(a - m)))
def trap_log(v, lf):
    w = np.diff(v); m = np.maximum(lf[:-1], lf[1:])
    return lse(np.log(w) + m + np.log(0.5*(np.exp(lf[:-1] - m) + np.exp(lf[1:] - m))))
def grid(ldelta):
    Vc = -ldelta
    v = np.unique(np.concatenate([np.linspace(LT0, LT0 + 2, 400), np.linspace(LT0, max(LT0 + 1, Vc + 3), 3000),
                                  np.linspace(max(LT0, Vc - 2), max(LT0, Vc) + 8, 4000), np.linspace(LT0, max(LT0, Vc) + 150, 3000)]))
    return v[v >= LT0]
def logF_kln(s, ldelta, dmax):
    c0, c1, smin, k, q0, q1 = best_row(s)
    lA = c0 + c1*s; Bc = q0 - q1*s; a = 8*(1 - s)/3
    lt = grid(ldelta); lkt = math.log(k) + lt
    lmain = lA + a*lt + 2*s*np.log(lkt) + (5 - 4*s)*np.log(lt) + np.log(a + 2*s/lkt + (5 - 4*s)/lt)
    ldN = np.logaddexp(lmain, math.log(Bc) + math.log(2) + np.log(lt))
    lI = trap_log(lt, logG(np.exp(lt + ldelta), dmax) + ldN)
    lN0 = math.log(math.exp(lA + a*LT0 + 2*s*math.log(math.log(k) + LT0) + (5 - 4*s)*math.log(LT0)) + Bc*LT0**2)
    return float(np.logaddexp(lN0 + float(logG(np.array([math.exp(LT0 + ldelta)]), dmax)[0]), lI))
def logF_cs(s, ldelta, dmax):
    s1, s2, B1, B2, B3 = cs_row(s)
    a = 3*(1 - s)/(2 - s); e = (7 - 5*s)/(2 - s)
    lt = grid(ldelta)
    lmain = math.log(B1) + a*lt + e*np.log(lt) + np.log(a + e/lt)
    ldN = np.logaddexp(lmain, np.log(2*B2*lt + B3))
    lI = trap_log(lt, logG(np.exp(lt + ldelta), dmax) + ldN)
    lN0 = math.log(B1*math.exp(a*LT0)*LT0**e + B2*LT0**2 + B3*LT0)
    return float(np.logaddexp(lN0 + float(logG(np.array([math.exp(LT0 + ldelta)]), dmax)[0]), lI))
def logS(ldelta, dmax):
    Vc = -ldelta
    v = np.unique(np.concatenate([np.linspace(math.log(14), Vc + 10, 6000), np.linspace(Vc - 3, Vc + 160, 4000)]))
    v = v[v >= math.log(14)]
    lf = logG(np.exp(v + ldelta), dmax) + np.log(1 + v) + v
    return float(np.logaddexp(math.log(14*math.log(14)) + float(logG(np.array([14*math.exp(ldelta)]), dmax)[0]), trap_log(v, lf)))
def R_direct(delta, L, use_cs):
    dmax = [float(K.dmax.mid()) for K in S0.kernels[2] if delta < float(K.dmax.mid())][0]
    ldelta = math.log(delta); s0 = 0.515
    ss = np.linspace(s0, B, 500)
    lf = []
    for s in ss:
        f = logF_kln(s, ldelta, dmax)
        if use_cs: f = min(f, logF_cs(s, ldelta, dmax))
        lf.append(math.log(L) + (s - 1)*L + f)
    Z = 2*(math.exp((s0 - 1)*L + logS(ldelta, dmax)) + math.exp(trap_log(ss, np.array(lf))))
    J = math.floor((L + delta)/math.log(2))
    PP = 2*(L + delta)*(math.exp(delta/2 - L/2)/4 + J*math.exp(delta - 2*L/3)/9 + math.log(J)*math.exp(-L)/delta)
    return Z + PP + 2*math.exp(-3*L)
def Nfun(s, lt, use_cs):
    c0, c1, smin, kk, q0, q1 = best_row(s)
    v = math.exp(c0 + c1*s + 8*(1-s)/3*lt + 2*s*math.log(math.log(kk) + lt) + (5-4*s)*math.log(lt)) + (q0 - q1*s)*lt**2
    if use_cs:
        s1, s2, B1, B2, B3 = cs_row(s)
        v = min(v, B1*math.exp(3*(1-s)/(2-s)*lt)*lt**((7-5*s)/(2-s)) + B2*lt**2 + B3*lt)
    return v
def band_breakdown(delta, L, use_cs, dmax=1e-9):
    """contributions to Z (float) by sigma band and height band: 2 log x int x^{s-1} int N(s,t) d(-G(delta t)) ds (right-end Stieltjes sums)"""
    lines = ['  L=%.2f delta*T0=%.3f (%s)' % (L, delta*3e12, 'KLN+CS' if use_cs else 'KLN')]
    for (a, b) in [(0.515, 0.80), (0.80, 0.83), (0.83, 0.86), (0.86, 0.87), (0.87, B)]:
        row = []
        for (t1, t2) in [(3e12, 1e13), (1e13, 1e15), (1e15, 1e40)]:
            vals = []; ss = np.linspace(a, b, 100)
            for s in ss:
                v = np.linspace(math.log(t1), math.log(t2), 2400); G = np.exp(logG(np.exp(v)*delta, dmax))
                Nv = np.interp(v, v[::40], np.array([Nfun(s, x, use_cs) for x in v[::40]]))
                vals.append(L*math.exp((s - 1)*L)*np.sum(Nv[1:]*(-np.diff(G))))
            row.append(2*np.trapezoid(vals, ss))
        lines.append('    sigma [%.3f,%.4f]: T in [3e12,1e13] %.3e  [1e13,1e15] %.3e  [1e15,oo) %.3e' % (a, b, *row))
    return lines
CASES = None
if __name__ == '__main__':
    out = []
    for a in sorted({0.5, A_K}):
        w = kernel_random_test(a)
        out.append('kernel a=%.3f: max over 60000 random z of (|Phi~(z)|/int Phi)/G(Im z) = %.12f (must be <= 1)' % (a, w)); print(out[-1], flush=True)
        assert w <= 1
    import json
    res = json.load(open(os.path.join(HERE, 'logs', 'cert_result.json')))
    km = 'main' if PIPE == 'A1' else 'A0'          # the certified theta of this pipeline
    thm, thc = float(res[km]['theta']), float(res['cs']['theta'])
    Lm, Lc, L6 = res[km]['worstL'], res['cs']['worstL'], res['cs_k6_alln']['worstL']
    for (th, L, cs) in ([(thm, Lm, False), (thc, Lc, True)] if PIPE == 'A1' else [(thm, Lm, False)]):
        out.append('band breakdown theta=%s:' % th); out += band_breakdown(math.log1p(math.exp((th - 1)*L)), L, cs); print(chr(10).join(out[-7:]), flush=True)
    cases = [('theta', thm, False, L) for L in (42, 100, 150, 185, Lm, 200, 230, 280)]
    if PIPE == 'A1':                               # the conditional runs of the certificate use pipeline A1
        cases += [('theta', thc, True, L) for L in (100, Lc, 200, 250)] + [('k', 6, True, L) for L in (40, 120, L6, 250)] + \
                 [('k', k, True, res['cs_k%d' % k]['L0']) for k in (5, 4, 3)]
    for (kind, par, cs, L) in cases:
        S = S1 if cs else S0
        if kind == 'theta':
            dfun = cert.delta_theta(arb(int(round(par*1000)))/1000); delta = math.log1p(math.exp((par - 1)*L))
        else:
            dfun = cert.delta_pow(par); delta = par*math.log1p(math.exp(-L/par))
        Rc = float(cert.cell_bound(S, dfun, L, L + 1e-6, 2).mid())
        Rd = R_direct(delta, L, cs)
        line = '%s=%s %s L=%.2f  cert R=%.5g   direct R=%.5g   ratio cert/direct=%.4f' % (kind, par, 'KLN+CS' if cs else 'KLN', L, Rc, Rd, Rc/Rd)
        print(line, flush=True); out.append(line)
    open(os.path.join(HERE, 'logs', 'check.log' if PIPE == 'A1' else 'check_A0.log'), 'w').write('\n'.join(out) + '\n')
