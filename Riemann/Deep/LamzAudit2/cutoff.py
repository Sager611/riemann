# cutoff.py -- lamz2audit: concrete smooth admissible psi_d for which C(f_d) - C_MT > d (quote of Lemma 3.2 in 2610.26841),
# and the gentle cutoffs of 2609.02882's proof converging to C_MT from above. Float (numpy Gauss-Legendre panels).
# Run: uv run --no-project --with numpy python3 cutoff.py > logs/cutoff.log
import numpy as np
CMT = 1.32749929632058835
s2 = np.sqrt(2.0); A0 = 1/(s2*np.sin(1/s2))
def S(t):
    t = np.clip(t, 0, 1); out = np.zeros_like(t)
    m = (t > 0) & (t < 1); tm = t[m]
    a = np.exp(-1/tm); b = np.exp(-1/(1-tm)); out[m] = a/(a+b); out[t >= 1] = 1.0
    return out
def make(d, kap):
    b = 0.5 - d + kap*d          # psi = 0 for |x| >= b < 1/2: compact support in (-1/2,1/2)
    def psi(x): return S((b - np.abs(x))/(kap*d))
    return psi, b
xg, wg = np.polynomial.legendre.leggauss(30)
def panels(bps, n_each):
    edges = []
    for i in range(len(bps)-1):
        edges += list(np.linspace(bps[i], bps[i+1], n_each+1)[:-1])
    edges.append(bps[-1]); return np.array(edges)
def C_of(fun, edges):
    # nodes, weights
    X = []; W = []
    for i in range(len(edges)-1):
        lo, hi = edges[i], edges[i+1]; X.append((hi-lo)/2*xg + (hi+lo)/2); W.append((hi-lo)/2*wg)
    X = np.concatenate(X); W = np.concatenate(W); fx = fun(X)
    M = np.sum(W*fx); fx = fx/M; scale = 1/M
    # cumulative F(u), G(u) at nodes: previous panels + partial panel integral
    F = np.zeros_like(X); G = np.zeros_like(X); k = 0; accF = 0.0; accG = 0.0
    for i in range(len(edges)-1):
        lo, hi = edges[i], edges[i+1]
        for j in range(len(xg)):
            u = X[k]; xs = (u-lo)/2*xg + (u+lo)/2; ws = (u-lo)/2*wg; fv = fun(xs)*scale
            F[k] = accF + np.sum(ws*fv); G[k] = accG + np.sum(ws*xs*fv); k += 1
        xs = (hi-lo)/2*xg + (hi+lo)/2; ws = (hi-lo)/2*wg; fv = fun(xs)*scale
        accF += np.sum(ws*fv); accG += np.sum(ws*xs*fv)
    Gtot = accG
    H = 2*X*F - X*1.0 - 2*G + Gtot          # int |u-v| f(v) dv for total mass 1
    return np.sum(W*fx*fx) + np.sum(W*fx*H), M
f0 = lambda x: np.where(np.abs(x) <= 0.5, A0*np.cos(s2*x), 0.0)
# sanity: f0 itself
C, M = C_of(f0, panels([-0.5, 0.5], 40)); print("sanity C(f0) =", "%.15f" % C, " C_MT =", "%.15f" % CMT, " diff %.2e" % (C-CMT))
# indicator cutoff a=1/2-d via quadrature vs closed form a + cot(sqrt2 a)/sqrt2
for d in [0.1, 0.01]:
    a = 0.5-d; Cind, _ = C_of(lambda x: np.where(np.abs(x) <= a, A0*np.cos(s2*x), 0.0), panels([-a, a], 40))
    print("indicator d=%g: quad C-C_MT = %.12f closed %.12f" % (d, Cind-CMT, a + 1/(np.tan(s2*a)*s2) - CMT))
print("smooth admissible psi_d (=1 on |x|<=1/2-d, 0 for |x|>=1/2-d+kap*d, C^inf step S):")
for d in [0.09, 0.05, 0.01, 0.003]:
    for kap in [0.9, 0.5, 0.1, 0.02]:
        psi, b = make(d, kap)
        fun = lambda x, psi=psi: f0(x)*psi(x)**2
        bps = [-b, -(0.5-d), 0.5-d, b]
        C, M = C_of(fun, panels(bps, 30))
        D = C - CMT
        print("  d=%-6g kap=%-5g A_d=%.6f  C-C_MT=%.8f  (C-C_MT)/d=%.5f  %s" % (d, kap, M, D, D/d, "VIOLATES <d" if D > d else "ok <d"))
