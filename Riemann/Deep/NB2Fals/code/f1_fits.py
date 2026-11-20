# nb2fals f1: regression of per-N C2(N) (and Q = C2 - (J - Sigma_N), kappa) on const [+ C3/L] + almost periodic zero terms:
# single-zero N^{-1/2} cos/sin(gamma_k log N), zero-pair cos/sin((gamma_j - gamma_i) log N), first K zeros (NB2nd derivation.md sec. 6, NBRate D(N)).
# Weights 1/N (uniform in u = log N). c_hat = sum_N a_N y_N (linear functional) -> errors:
#  (i) moving-block bootstrap of the residual function in u (block length b), (ii) out-of-sample, (iii) window stability,
#  (iv) leakage of unmodelled pair/single terms of zeros <= 100 through the transfer function T(w) = sum_N a_N e^{i w u_N},
#       amplitudes from the theory shape s_ab = 1/(gamma_a gamma_b) (pairs), s_k = 1/(gamma_k^2 |zeta'(rho_k)|) (singles), scale tau fitted at K = 5.
# usage: f1_fits.py NMAX
import sys, math, numpy as np
O = '/Users/adriansager/playground/riemann/.scratch/nb2fals/'
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 4000
P = np.load(O + f'perN{NMAX}.npz')['perN']
Z = np.load(O + 'zeros100.npz'); gam = Z['gam']; zp = np.abs(Z['zp'])
Yall = {'C2': P[:, 2], 'Q': P[:, 2] - P[:, 4], 'kappa': P[:, 3]}
Nall, uall = P[:, 0], P[:, 1]
rng = np.random.default_rng(12345)
def lines(K, Kz=None):
    sg = [(gam[k], 1 / (gam[k] ** 2 * zp[k])) for k in range(K)]
    pr = [(gam[j] - gam[i], 1 / (gam[i] * gam[j])) for i in range(K) for j in range(i + 1, K)]
    return sg, pr
def design(uu, K, trend):
    cols = [np.ones_like(uu)] + ([1 / uu] if trend else [])
    sg, pr = lines(K); e = np.exp(-uu / 2)
    for g, _ in sg: cols += [e * np.cos(g * uu), e * np.sin(g * uu)]
    for w, _ in pr: cols += [np.cos(w * uu), np.sin(w * uu)]
    return np.column_stack(cols)
def ols(uu, NN, y, K, trend):
    w = (1 / NN) / np.sum(1 / NN); sw = np.sqrt(w); X = design(uu, K, trend); Xw = X * sw[:, None]
    Qm, Rm = np.linalg.qr(Xw); e0 = np.zeros(X.shape[1]); e0[0] = 1
    beta = np.linalg.solve(Rm, Qm.T @ (sw * y)); a = sw * (Qm @ np.linalg.solve(Rm.T, e0))
    sv = np.linalg.svd(Rm, compute_uv=False)
    return beta, a, y - X @ beta, sv[0] / sv[-1], X
def boot(uu, res, a, b, B=400):
    lo, hi = uu[0], uu[-1]; nb = int(math.ceil((hi - lo) / b)); ib = np.minimum(((uu - lo) / b).astype(int), nb - 1)
    st = rng.uniform(lo, hi - b, size=(B, nb)); src = st[:, ib] + (uu - lo - ib * b)[None, :]
    return float(np.std(np.interp(src, uu, res) @ a))
def leak(a, uu, K, tau_s, tau_p, Kmax=100):
    # rms leakage into c_hat of the single/pair lines of zeros <= Kmax not in the model (random phases)
    tot = 0.0
    for k in range(K, Kmax):
        T = abs(np.sum(a * np.exp(-uu / 2) * np.exp(1j * gam[k] * uu))); tot += (tau_s / (gam[k] ** 2 * zp[k])) ** 2 * T ** 2
    for j in range(1, Kmax):
        for i in range(j):
            if j < K: continue
            T = abs(np.sum(a * np.exp(1j * (gam[j] - gam[i]) * uu))); tot += (tau_p / (gam[i] * gam[j])) ** 2 * T ** 2
    return math.sqrt(tot)
def sel(lo, hi): return (Nall >= lo) & (Nall <= hi)
def tau_from(beta, K, trend):
    o = 2 if trend else 1; sg, pr = lines(K); cs = beta[o:o + 2 * K].reshape(-1, 2); cp = beta[o + 2 * K:].reshape(-1, 2)
    ts = math.sqrt(np.mean([(c[0] ** 2 + c[1] ** 2) / (2 * s ** 2) for c, (g, s) in zip(cs, sg)])) if K else 0
    tp = math.sqrt(np.mean([(c[0] ** 2 + c[1] ** 2) / (2 * s ** 2) for c, (w, s) in zip(cp, pr)])) if K > 1 else 0
    return ts, tp
b20 = lambda m: 20 / np.median(Nall[m])   # "20-N blocks" in u at the window's median N
print(f'NMAX = {NMAX}; weights 1/N; bootstrap B = 400; SE_b = moving-block bootstrap s.e. with block length b in u')
for tgt, trend in (('C2', True), ('Q', False), ('Q', True), ('kappa', False), ('kappa', True)):
    y_all = Yall[tgt]
    for (lo, hi) in ((100, 640), (100, NMAX)):
        m = sel(lo, hi); uu, NN, y = uall[m], Nall[m], y_all[m]
        _, _, _, _, _ = (None,) * 5
        b5, _, _, _, _ = ols(uu, NN, y, 5, trend); ts, tp = tau_from(b5, 5, trend)
        print(f'\n== {tgt} {"+ C3/L " if trend else ""}on N in [{lo},{hi}] ({m.sum()} N, u-length {uu[-1]-uu[0]:.2f}); theory-shape scales from K=5: tau_single {ts:.3f}, tau_pair {tp:.3f}')
        print(' K    p   const     C3/L   | SE_b(20N)  SE_b(0.1) SE_b(0.3) SE_b(0.6) | leak(z<=100) | cond(Xw)  | in-rms   | const/SE(0.6)')
        for K in (0, 5, 10, 20, 50):
            p = 1 + trend + 2 * K + K * (K - 1)
            if p > 0.8 * m.sum(): print(f'{K:2d} {p:5d}  not identifiable (p > 0.8 n)'); continue
            if K == 50: print(f'{K:2d} {p:5d}  OLS skipped (1225 pair lines, spacing ~0.1 << resolution {2*math.pi/(uu[-1]-uu[0]):.2f}); see f2 ridge'); continue
            beta, a, res, cond, X = ols(uu, NN, y, K, trend)
            ses = [boot(uu, res, a, b) for b in (b20(m), 0.1, 0.3, 0.6)]
            lk = leak(a, uu, K, ts, tp) if K < 50 else float('nan')
            win = np.sqrt(np.sum(res ** 2 / NN) / np.sum(1 / NN))
            c3 = f'{beta[1]:+.4f}' if trend else '   -   '
            print(f'{K:2d} {p:5d} {beta[0]:+.5f} {c3} | {ses[0]:.5f}   {ses[1]:.5f}   {ses[2]:.5f}   {ses[3]:.5f}  | {lk:.5f}     | {cond:.1e} | {win:.2e} | {beta[0]/math.hypot(ses[3], lk):+.1f}')
