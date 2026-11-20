# nb2fals f3: robustness of the GP fit (f2) for C2 + C3/L and Q on N in [100,4000]:
#  (a) fitted shape s_ab = (gamma_1 gamma_2/(gamma_a gamma_b))^alpha (g_12/g_ab)^beta (alpha, beta by REML), K = 50, 100;
#  (b) extra 1/u-modulated copies of the pair lines (next-order AP term C3(N)/L), scale tau_m by REML;
#  (c) trend vs no trend by ML (likelihood ratio), Q and C2.
import sys, math, numpy as np
from scipy.linalg import cho_factor, cho_solve, toeplitz
from scipy.optimize import minimize
O = '/Users/adriansager/playground/riemann/.scratch/nb2fals/'
import os; NMX = int(os.environ.get('NB2F_NMAX', '4000')); P = np.load(O + f'perN{NMX}.npz')['perN']; Z = np.load(O + 'zeros100.npz'); gam = Z['gam']; zp = np.abs(Z['zp'])
Nall, uall = P[:, 0], P[:, 1]; Y = {'C2': P[:, 2], 'Q': P[:, 2] - P[:, 4]}
d = 0.005
def grid(y, lo, hi):
    m = (Nall >= lo) & (Nall <= hi); ug = np.arange(uall[m][0], uall[m][-1] + 1e-12, d); return ug, np.interp(ug, uall[m], y[m])
CACHE = {}
def basis(ug, K):
    key = (len(ug), ug[0], K)
    if key not in CACHE:
        n = len(ug); lag = np.arange(n) * d; i, j = np.triu_indices(K, 1); om = gam[j] - gam[i]
        CACHE[key] = (np.cos(np.outer(om, lag)), np.log(gam[0] * gam[1] / (gam[i] * gam[j])), np.log((gam[1] - gam[0]) / om),
                      toeplitz(((1 / (gam[:K] ** 2 * zp[:K])) ** 2) @ np.cos(np.outer(gam[:K], lag))) * np.outer(np.exp(-ug / 2), np.exp(-ug / 2)))
    return CACHE[key]
def covm(ug, K, lt_p, lt_s, ls0, al, be, lt_m=None):
    Cp, la, lb, Ks = basis(ug, K); s2 = np.exp(2 * (lt_p + al * la + be * lb)) / (gam[0] * gam[1]) ** 2
    Kp = toeplitz(s2 @ Cp); Km = Kp + math.exp(2 * lt_s) * Ks
    if lt_m is not None: Km = Km + math.exp(2 * (lt_m - lt_p)) * Kp / np.outer(ug, ug) * ug.mean() ** 2
    Km[np.diag_indices(len(ug))] += math.exp(2 * ls0)
    return Km
def gls(ug, yg, Km, trend, ml=False):
    X = np.column_stack([np.ones_like(ug)] + [ug ** (-p) for p in range(1, int(trend) + 1)])
    cf = cho_factor(Km, lower=True); KX = cho_solve(cf, X); A = X.T @ KX; Ai = np.linalg.inv(A); beta = Ai @ (KX.T @ yg)
    r = yg - X @ beta; q = r @ cho_solve(cf, r); ld = 2 * np.sum(np.log(np.diag(cf[0])))
    return beta, Ai, -0.5 * (ld + q + (0 if ml else np.linalg.slogdet(A)[1]))
def fit(tgt, K, trend, free_shape, modul=False, ml=False, lo=100, hi=None, p0=None, beta_fix=None):
    hi = hi or NMX
    ug, yg = grid(Y[tgt], lo, hi)
    def unpack(p):
        lt_p, lt_s, ls0 = p[0], p[1], p[2]; al, be = (p[3], p[4]) if free_shape else (1.0, 1.0 if beta_fix is None else beta_fix); lt_m = p[-1] if modul else None
        return lt_p, lt_s, ls0, al, be, lt_m
    def f(p):
        try: return -gls(ug, yg, covm(ug, K, *unpack(p)), trend, ml)[2]
        except np.linalg.LinAlgError: return 1e30
    p = list(p0) if p0 is not None else [math.log(0.7), math.log(8.0), math.log(1.6e-4)] + ([1.0, 1.0] if free_shape else []) + ([math.log(0.3)] if modul else [])
    sim = np.array(p) + np.vstack([np.zeros(len(p)), np.eye(len(p)) * 0.5])
    r = minimize(f, p, method='Nelder-Mead', options={'initial_simplex': sim, 'xatol': 1e-2, 'fatol': 1e-2, 'maxiter': 150 * len(p)})
    beta, Ai, ll = gls(ug, yg, covm(ug, K, *unpack(r.x)), trend, ml)
    return beta, np.sqrt(np.diag(Ai)), ll, r.x
mode = sys.argv[1] if len(sys.argv) > 1 else 'a'
if mode == 'a':
    for K in (50, 100):
        for tgt in ('C2',):
            b, s, ll, x = fit(tgt, K, True, True)
            print(f'(a) {tgt}+C3/L K={K} fitted shape: c = {b[0]:+.5f} +- {s[0]:.5f}, C3 = {b[1]:+.4f} +- {s[1]:.4f}; REML {ll:.1f}; tau_p {math.exp(x[0]):.3g} tau_s {math.exp(x[1]):.3g} s0 {math.exp(x[2]):.2e} alpha {x[3]:.2f} beta {x[4]:.2f}', flush=True)
            b, s, ll, x = fit(tgt, K, True, False)
            print(f'(a) {tgt}+C3/L K={K} shape B:      c = {b[0]:+.5f} +- {s[0]:.5f}, C3 = {b[1]:+.4f} +- {s[1]:.4f}; REML {ll:.1f}', flush=True)
if mode == 'b':
    for tgt, trend in (('C2', True), ('Q', False), ('Q', True)):
        b, s, ll, x = fit(tgt, 50, trend, False, modul=True)
        print(f'(b) {tgt}{"+C3/L" if trend else ""} K=50 shape B + 1/u-modulated pair lines: c = {b[0]:+.5f} +- {s[0]:.5f}' + (f', C3 = {b[1]:+.4f} +- {s[1]:.4f}' if trend else '') + f'; REML {ll:.1f}; tau_p {math.exp(x[0]):.3g} tau_m {math.exp(x[-1]):.3g} s0 {math.exp(x[2]):.2e}', flush=True)
if mode == 'c':
    for tgt in ('Q', 'C2'):
        out = []
        for trend in (False, True):
            b, s, ll, x = fit(tgt, 50, trend, False, ml=True); out.append(ll)
            print(f'(c) {tgt} K=50 shape B, ML fit, trend={trend}: c = {b[0]:+.5f} +- {s[0]:.5f}' + (f', C3 = {b[1]:+.4f} +- {s[1]:.4f}' if trend else '') + f'; logL_ML {ll:.2f}', flush=True)
        print(f'(c) {tgt}: 2 Delta logL (trend vs none) = {2*(out[1]-out[0]):.2f} (chi2_1: 3.84 at 5%, 9.0 at 3sigma)', flush=True)

if mode == 'd':   # next order of the 1/L expansion: C2 + C3/L + C4/L^2 (trend=2); also Q
    for tgt in ('C2', 'Q'):
        for fs in (False, True):
            b, s, ll, x = fit(tgt, 50, 2, fs)
            print(f'(d) {tgt}+C3/L+C4/L^2 K=50 shape {"fitted" if fs else "B"}: c = {b[0]:+.5f} +- {s[0]:.5f}, C3 = {b[1]:+.4f} +- {s[1]:.4f}, C4 = {b[2]:+.4f} +- {s[2]:.4f}; REML {ll:.1f}' + (f'; alpha {x[3]:.2f} beta {x[4]:.2f}' if fs else ''), flush=True)
        b, s, ll, x = fit(tgt, 50, 2, False, ml=True); b1, s1, ll1, x1 = fit(tgt, 50, 1, False, ml=True)
        print(f'(d) {tgt}: ML 2 Delta logL (1/L + 1/L^2 vs 1/L) = {2*(ll-ll1):.2f}', flush=True)

if mode == 'e':   # profile over the gap exponent beta (alpha = 1): REML-weighted mixture of c
    rows = []
    for be in (0.5, 0.65, 0.8, 0.9, 1.0, 1.2, 1.5):
        b, s, ll, x = fit('C2', 50, 1, False, beta_fix=be); rows.append((be, b[0], s[0], ll))
        print(f'(e) C2+C3/L K=50 alpha=1 beta={be}: c = {b[0]:+.5f} +- {s[0]:.5f}; REML {ll:.1f}', flush=True)
    R = np.array(rows); w = np.exp(R[:, 3] - R[:, 3].max()); w /= w.sum(); m = w @ R[:, 1]; v = w @ (R[:, 2] ** 2 + (R[:, 1] - m) ** 2)
    print(f'(e) REML-weighted mixture over beta: c = {m:+.5f} +- {math.sqrt(v):.5f}; weights {np.round(w, 3).tolist()}', flush=True)
