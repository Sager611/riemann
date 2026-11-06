#!/usr/bin/env python3
# check.py -- independent floating-point cross-check of cert.py's bound (NOT the certificate; numpy, log-space, trapezoid rules).
# Recomputes Z = 2[x^{s0-1} S + log x int_{s0}^{B} x^{s-1} F(s) ds] by direct quadrature, where
#   F(s) = KLN(s,T0) G(delta T0) + int_{T0}^oo G(delta t) dKLN(s,t)   (KLN with its exact (log kt)^{2s}(log t)^{5-4s} t^{8(1-s)/3} shape,
#   best admissible row at each s), S = 14 log 14 G(14 delta) + int_14^oo G(delta t)(1 + log t) dt, G = box-spline envelope (m = 2).
# usage: uv run --no-project --with python-flint --with mpmath --with sympy --with numpy python3 Riemann/Deep/Ingham3/check.py
import sys, os, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.argv.append('quick')
import cert
from flint import arb
S = cert.main()
B = 0.8733819666132076; LT0 = math.log(3e12); S1 = 0.2173
ROWS = [(float(P.c0.mid()), float(P.c1.mid()), float(P.smin.mid()), float(P.k.mid()), float(P.q0.mid()), float(P.q1.mid())) for P in S.P]
def best_row(s):
    c = [(c0 + c1*s, r) for r in ROWS for (c0, c1) in [(r[0], r[1])] if s > r[2]]
    return min(c)[1]
YC = 2.7328      # sin(y)/y = S1 (main-lobe end of the envelope); refined below
for _ in range(60):
    YC = YC - (math.sin(YC)/YC - S1)/((YC*math.cos(YC) - math.sin(YC))/YC**2)
def logG(lam, dmax):
    y = lam/4.0
    with np.errstate(divide='ignore', invalid='ignore'):
        sinc = np.where(y > 0, np.sin(y)/np.where(y > 0, y, 1), 1.0)
    env = np.where(y <= math.pi, np.maximum(sinc, S1), np.minimum(S1, 1/np.maximum(y, 1e-300)))
    corr = 2*dmax*np.minimum(0.5, 1/np.maximum(lam, 1e-300))
    return 2*np.log(env + corr)
def lse(a): m = np.max(a); return m + math.log(np.sum(np.exp(a - m)))
def trap_log(v, lf):   # log of int exp(lf) dv, trapezoid on grid v
    w = np.diff(v); m = np.maximum(lf[:-1], lf[1:])
    return lse(np.log(w) + m + np.log(0.5*(np.exp(lf[:-1] - m) + np.exp(lf[1:] - m))))
def logF(s, ldelta, dmax):
    c0, c1, smin, k, q0, q1 = best_row(s)
    lA = c0 + c1*s; Bc = q0 - q1*s; a = 8*(1 - s)/3
    Vc = -ldelta
    v = np.unique(np.concatenate([np.linspace(LT0, LT0 + 2, 400), np.linspace(LT0, max(LT0 + 1, Vc + 3), 3000),
                                  np.linspace(max(LT0, Vc - 2), max(LT0, Vc) + 8, 4000), np.linspace(LT0, max(LT0, Vc) + 150, 3000)]))
    v = v[v >= LT0]
    lt = v; lkt = np.log(np.exp(lt*0) * k) + lt   # log(kt) = log k + log t
    # d/dv [A (log kt)^{2s} (log t)^{5-4s} t^a] = A t^a (log kt)^{2s} (log t)^{5-4s} [a + 2s/log kt + (5-4s)/log t]
    lmain = lA + a*lt + 2*s*np.log(lkt) + (5 - 4*s)*np.log(lt) + np.log(a + 2*s/lkt + (5 - 4*s)/lt)
    lB = math.log(Bc) + math.log(2) + np.log(lt)
    ldN = np.logaddexp(lmain, lB)
    lg = logG(np.exp(lt + ldelta), dmax)
    lI = trap_log(v, lg + ldN)
    lN0 = math.log(math.exp(lA + a*LT0 + 2*s*math.log(math.log(k) + LT0) + (5 - 4*s)*math.log(LT0)) + Bc*LT0**2)
    lb = float(logG(np.array([math.exp(LT0 + ldelta)]), dmax)[0])
    return float(np.logaddexp(lN0 + lb, lI))
def logS(ldelta, dmax):
    Vc = -ldelta
    v = np.unique(np.concatenate([np.linspace(math.log(14), Vc + 10, 6000), np.linspace(Vc - 3, Vc + 160, 4000)]))
    v = v[v >= math.log(14)]
    lf = logG(np.exp(v + ldelta), dmax) + np.log(1 + v) + v
    return float(np.logaddexp(math.log(14*math.log(14)) + float(logG(np.array([14*math.exp(ldelta)]), dmax)[0]), trap_log(v, lf)))
def Z_direct(k, L, dmax):
    n = math.exp(L/k); ldelta = math.log(k*math.log1p(1/n)); s0 = 0.515
    ss = np.linspace(s0, B, 600)
    lf = np.array([math.log(L) + (s - 1)*L + logF(s, ldelta, dmax) for s in ss])
    lI = trap_log(ss, lf)
    return 2*(math.exp((s0 - 1)*L + logS(ldelta, dmax)) + math.exp(lI))
def kernel_random_test(ntest=20000, seed=1):
    """|Phi~(z)|/int Phi <= G(Im z) for Re z in [0, dmax] (m = 2, Phi = 1_[0,1/2]*1_[0,1/2]); exact Phi~(z) = ((e^{z/2}-1)/z)^2."""
    rng = np.random.default_rng(seed); worst = 0.0
    for dmax in (1e-9, 0.025, 0.15):
        lam = np.exp(rng.uniform(math.log(1e-3), math.log(1e5), ntest)); re = rng.uniform(0, dmax, ntest)
        z = re + 1j*lam
        ratio = np.abs(((np.exp(z/2) - 1)/z)**2)/0.25
        bound = np.exp(logG(lam, dmax))
        worst = max(worst, float(np.max(ratio/bound)))
    return worst
def ramare_b1_b2(N=10**7):
    """numerical sanity check (not a proof) of KLN's b1, b2 (Ramare Lemmas 3.1, 3.4) for X <= 1e7"""
    sq = np.ones(N + 1, dtype=bool); sq[0] = False
    for p in range(2, int(N**0.5) + 1): sq[p*p::p*p] = False
    X = np.arange(N + 1); cs = np.cumsum(sq)
    inv = np.where(sq, 1.0/np.maximum(X, 1), 0.0); inv[0] = 0; cinv = np.cumsum(inv)
    return float((cs[1700:]/X[1700:]).max()), float((cinv[1002:] - 6/math.pi**2*np.log(X[1002:])).max())
def kln_vs_trivial():
    lines = []
    def kln(s, lt):
        c0, c1, smin, k, q0, q1 = best_row(s)
        return math.exp(c0 + c1*s + 8*(1-s)/3*lt + 2*s*math.log(math.log(k) + lt) + (5-4*s)*math.log(lt)) + (q0 - q1*s)*lt**2
    def triv(lt):
        def Nm(l, sgn):
            t = math.exp(l); return t/(2*math.pi)*math.log(t/(2*math.pi*math.e)) + sgn*(0.1038*l + 0.2573*math.log(l) + 9.3675)
        return (Nm(lt, 1) - Nm(LT0, -1))/2
    for s in (0.75, 0.78, 0.80, 0.83, 0.85, 0.86, 0.87, 0.8734):
        lines.append('  sigma %.4f: ' % s + '  '.join('T=%.0e KLN %.3e trivial %.3e' % (math.exp(lt), kln(s, lt), triv(lt)) for lt in (math.log(1e13), math.log(1e14))))
    return lines
def band_breakdown(k, L, dmax=1e-9):
    """contributions to Z (float) by sigma band and by height band, int KLN d(-G) with right-end Stieltjes sums"""
    n = math.exp(L/k); delta = k*math.log1p(1/n); lines = ['k=%d L=%.2f n=%.4e lambda0=delta*T0=%.3f' % (k, L, n, delta*3e12)]
    def kln(s, lt):
        c0, c1, smin, kk, q0, q1 = best_row(s)
        return math.exp(c0 + c1*s + 8*(1-s)/3*lt + 2*s*math.log(math.log(kk) + lt) + (5-4*s)*math.log(lt)) + (q0 - q1*s)*lt**2
    for (a, b) in [(0.515, 0.80), (0.80, 0.83), (0.83, 0.86), (0.86, 0.87), (0.87, B)]:
        row = []
        for (t1, t2) in [(3e12, 1e13), (1e13, 1e15), (1e15, 1e40)]:
            vals = []; ss = np.linspace(a, b, 120)
            for s in ss:
                v = np.linspace(math.log(t1), math.log(t2), 2400); G = np.exp(logG(np.exp(v)*delta, dmax))
                Nv = np.interp(v, v[::40], np.array([kln(s, x) for x in v[::40]]))
                vals.append(L*math.exp((s - 1)*L)*np.sum(Nv[1:]*(-np.diff(G))))
            row.append(2*np.trapezoid(vals, ss))
        lines.append('  sigma [%.3f,%.4f]: T in [3e12,1e13] %.3e  [1e13,1e15] %.3e  [1e15,oo) %.3e' % (a, b, *row))
    return lines
if __name__ == '__main__':
    out = []
    w = kernel_random_test()
    out.append('kernel: max over 60000 random z of (|Phi~(z)|/int Phi)/G(Im z) = %.12f (must be <= 1)' % w); print(out[-1], flush=True)
    r1, r2 = ramare_b1_b2(); out.append('Ramare b1: max_{1700<=X<=1e7} sum mu^2(n)/X = %.6f <= 0.62;  b2: max_{1002<=X<=1e7} sum mu^2(n)/n - (6/pi^2)log X = %.6f <= 1.048' % (r1, r2)); print(out[-1], flush=True)
    out.append('KLN (our evaluation, best admissible row) vs trivial count (N(T)-N(T0))/2:'); out += kln_vs_trivial(); print(chr(10).join(out[-9:]), flush=True)
    for (k, L) in [(7, 201.87), (7, 205.5)]:
        out += band_breakdown(k, L); print(chr(10).join(out[-6:]), flush=True)
    for (k, L) in [(7, 40), (7, 100), (7, 180), (7, 201.8), (7, 205.5), (7, 230), (6, 217.45), (5, 288.73), (4, 468.84), (3, 1881.03)]:
        Rc = cert.cell_bound(S, k, L, L + 1e-6, 2)
        n = math.exp(L/k); delta = k*math.log1p(1/n)
        dmax = [float(K.dmax.mid()) for K in S.kernels[2] if delta < float(K.dmax.mid())][0]
        Zd = Z_direct(k, L, dmax)
        line = 'k=%d L=%.2f  cert R=%.5g   direct Z=%.5g   ratio cert/direct=%.4f' % (k, L, float(Rc.mid()), Zd, float(Rc.mid())/Zd)
        print(line, flush=True); out.append(line)
    open(os.path.join(HERE, 'logs', 'check.log'), 'w').write('\n'.join(out) + '\n')
