# nb2fals f2: Gaussian-process (random-amplitude) version of the joint fit: y(u) = c [+ C3/u] + sum of zero lines with random amplitudes
#   pairs (a<b<=K): cos/sin((gamma_b - gamma_a) u), coefficient variance (tau_p s_ab)^2; singles: e^{-u/2} cos/sin(gamma_k u), (tau_s s_k)^2;
#   + white nugget s0^2.  c, C3 flat (GLS). Hyperparameters by REML. Data linearly interpolated onto a uniform u grid (step d).
# Shapes: A: s_ab = 1/(gamma_a gamma_b);  B: s_ab = (gamma_2-gamma_1)/(gamma_a gamma_b |gamma_b-gamma_a|) (close pairs enhanced). s_k = 1/(gamma_k^2|zeta'(rho_k)|).
# Outputs: c_hat +- posterior sd; REML; out-of-sample z-scores; posterior sd of c for longer (hypothetical) windows (data-free).
import sys, math, numpy as np
from scipy.linalg import cho_factor, cho_solve, toeplitz
from scipy.optimize import minimize
O = '/Users/adriansager/playground/riemann/.scratch/nb2fals/'
import os; NMX = int(os.environ.get('NB2F_NMAX', '4000')); P = np.load(O + f'perN{NMX}.npz')['perN']; Z = np.load(O + 'zeros100.npz'); gam = Z['gam']; zp = np.abs(Z['zp'])
Nall, uall = P[:, 0], P[:, 1]; Y = {'C2': P[:, 2], 'Q': P[:, 2] - P[:, 4]}
def lines(K, shape):
    i, j = np.triu_indices(K, 1); om = gam[j] - gam[i]; s = 1 / (gam[i] * gam[j])
    if shape == 'B': s = s * (gam[1] - gam[0]) / om
    return om, s, gam[:K], 1 / (gam[:K] ** 2 * zp[:K])
def covm(ug, d, K, shape, tp, ts, s0):
    n = len(ug); lag = np.arange(n) * d; om, s, gs, ss = lines(K, shape)
    kp = ((tp * s) ** 2) @ np.cos(np.outer(om, lag)) if K > 1 else np.zeros(n)
    Km = toeplitz(kp)
    if K > 0:
        e = np.exp(-ug / 2); Km += toeplitz(((ts * ss) ** 2) @ np.cos(np.outer(gs, lag))) * np.outer(e, e)
    Km[np.diag_indices(n)] += s0 ** 2
    return Km
def gls(ug, yg, Km, trend):
    X = np.column_stack([np.ones_like(ug)] + ([1 / ug] if trend else []))
    cf = cho_factor(Km, lower=True); Ki_X = cho_solve(cf, X); A = X.T @ Ki_X; Ai = np.linalg.inv(A)
    beta = Ai @ (Ki_X.T @ yg); r = yg - X @ beta; q = r @ cho_solve(cf, r)
    ld = 2 * np.sum(np.log(np.diag(cf[0])))
    reml = -0.5 * (ld + np.linalg.slogdet(A)[1] + q)
    return beta, Ai, reml, cf, X
def grid(y, lo, hi, d):
    m = (Nall >= lo) & (Nall <= hi); ug = np.arange(uall[m][0], uall[m][-1] + 1e-12, d); return ug, np.interp(ug, uall[m], y[m])
def fit(y, lo, hi, K, shape, trend, d=0.005, x0=(0.0, 0.0, math.log(3e-4))):
    ug, yg = grid(y, lo, hi, d)
    if K == 0:
        f = lambda p: -gls(ug, yg, covm(ug, d, 0, shape, 0, 0, math.exp(p[0])), trend)[2]
        r = minimize(f, [x0[2]], method='Nelder-Mead', options={'initial_simplex': [[x0[2]], [x0[2] + 1.0]], 'xatol': 1e-3, 'fatol': 1e-3}); hp = (0, 0, math.exp(r.x[0]))
    else:
        def f(p):
            try: return -gls(ug, yg, covm(ug, d, K, shape, *np.exp(p)), trend)[2]
            except np.linalg.LinAlgError: return 1e30
        best = min(((f([a, b, c]), (a, b, c)) for a in np.arange(-3, 5.1, 1.0) for b in (-2.0, 1.0, 4.0) for c in np.log([1e-5, 1e-4, 1e-3])), key=lambda t: t[0])[1]
        sim = np.array(best) + np.vstack([np.zeros(3), np.eye(3) * 0.7])
        r = minimize(f, list(best), method='Nelder-Mead', options={'initial_simplex': sim, 'xatol': 1e-2, 'fatol': 1e-2, 'maxiter': 300}); hp = tuple(np.exp(r.x))
    beta, Ai, reml, cf, X = gls(ug, yg, covm(ug, d, K, shape, *hp), trend)
    return dict(beta=beta, sd=np.sqrt(np.diag(Ai)), reml=reml, hp=hp, ug=ug, yg=yg, d=d)
def predict(F, y, lo, hi, K, shape, trend):
    # posterior predictive on (lo, hi] given the fit window data; z-scores of the window mean and pointwise rms z
    ug2, yg2 = grid(y, lo, hi, F['d']); ug = F['ug']; d = F['d']; ua = np.concatenate([ug, ug2])
    # covariance between arbitrary grid points: rebuild on the union grid (uniform step d, contiguous if windows adjoin)
    n1 = len(ug); g0 = ug[0]; idx = np.round((ua - g0) / d).astype(int); n = idx.max() + 1; ufull = g0 + np.arange(n) * d
    Kf = covm(ufull, d, K, shape, *F['hp']); K11 = Kf[np.ix_(idx[:n1], idx[:n1])]; K21 = Kf[np.ix_(idx[n1:], idx[:n1])]; K22 = Kf[np.ix_(idx[n1:], idx[n1:])]
    X1 = np.column_stack([np.ones(n1)] + ([1 / ug] if trend else [])); X2 = np.column_stack([np.ones(len(ug2))] + ([1 / ug2] if trend else []))
    cf = cho_factor(K11, lower=True); b = F['beta']; r1 = F['yg'] - X1 @ b
    mean2 = X2 @ b + K21 @ cho_solve(cf, r1)
    # predictive covariance incl. beta uncertainty
    Ai = np.diag(F['sd']) ** 2 if not trend else None
    R2 = X2 - K21 @ cho_solve(cf, X1); Ki_X1 = cho_solve(cf, X1); A = X1.T @ Ki_X1; Ain = np.linalg.inv(A)
    C2 = K22 - K21 @ cho_solve(cf, K21.T) + R2 @ Ain @ R2.T
    w = np.ones(len(ug2)) / len(ug2); mres = w @ (yg2 - mean2); msd = math.sqrt(w @ C2 @ w)
    z = (yg2 - mean2) / np.sqrt(np.diag(C2))
    return mres / msd, float(np.sqrt(np.mean(z ** 2))), float(np.sqrt(np.mean((yg2 - mean2) ** 2))), float(np.std(yg2))
def sd_window(hp, lo_u, hi_u, K, shape, trend, d=0.02):  # d = 0.02: Nyquist 157 > all line frequencies of zeros <= 50
    ug = np.arange(lo_u, hi_u + 1e-12, d); Km = covm(ug, d, K, shape, *hp)
    X = np.column_stack([np.ones_like(ug)] + ([1 / ug] if trend else [])); cf = cho_factor(Km, lower=True)
    return math.sqrt(np.linalg.inv(X.T @ cho_solve(cf, X))[0, 0])
mode = sys.argv[1] if len(sys.argv) > 1 else 'main'
if mode == 'main':
    print('target trend shape K | c_hat +- sd | C3 | REML | tau_p tau_s s0 | OOS: z(mean), rms z, rms err/sd(y)')
    for tgt, trend in (('C2', True), ('Q', False), ('Q', True)):
        for (lo, hi, plo, phi) in ((100, 640, None, None), (100, 4000, None, None)):
            for shape in ('A', 'B'):
                for K in (0, 5, 10, 20, 50):
                    if K == 0 and shape == 'B': continue
                    F = fit(Y[tgt], lo, hi, K, shape, trend)
                    c3 = f'{F["beta"][1]:+.4f}' if trend else '   -   '
                    print(f'{tgt:2s} {"C3/L" if trend else "none"} {shape} {K:2d} [{lo},{hi}] | {F["beta"][0]:+.5f} +- {F["sd"][0]:.5f} | {c3} | {F["reml"]:.1f} | {F["hp"][0]:.3g} {F["hp"][1]:.3g} {F["hp"][2]:.2e}', flush=True)

if mode == 'checks':
    print('window stability, out-of-sample, projection; K = 50 lines (zeros <= 50), both shapes')
    for tgt, trend in (('C2', True), ('Q', False)):
        for shape in (sys.argv[2] if len(sys.argv) > 2 else 'AB'):
            K = 50
            for (lo, hi) in ((100, 640), (100, 1280), (100, 2000), (100, 4000), (100, NMX), (300, NMX), (640, NMX), (1000, NMX)):
                if hi == NMX and NMX == 4000 and lo == 100 and (lo, hi) in seen: continue
                seen = locals().get('seen', set()) | {(lo, hi)}
                F = fit(Y[tgt], lo, hi, K, shape, trend)
                print(f'{tgt} {"C3/L" if trend else "none"} {shape} K={K} window [{lo},{hi}]: c = {F["beta"][0]:+.5f} +- {F["sd"][0]:.5f}' + (f', C3 = {F["beta"][1]:+.4f} +- {F["sd"][1]:.4f}' if trend else '') + f'; hp {F["hp"][0]:.3g} {F["hp"][1]:.3g} {F["hp"][2]:.2e}', flush=True)
                if (lo, hi) == (100, NMX): Ffull = F
            for (lo, mid, hi) in ((100, 400, 640), (100, 1280, NMX), (100, 2000, NMX), (100, 4000, NMX)):
                if mid >= hi: continue
                F = fit(Y[tgt], lo, mid, K, shape, trend)
                zm, zr, er, sy = predict(F, Y[tgt], mid + 1, hi, K, shape, trend)
                print(f'   OOS fit [{lo},{mid}] -> ({mid},{hi}]: z(window mean) {zm:+.2f}, rms pointwise z {zr:.2f}, rms error {er:.5f} vs sd(y) {sy:.5f}', flush=True)
            for uh in (math.log(4000), 10.0, 12.0, 15.0, 20.0):
                print(f'   projection (hp of [100,4000]): window [log 100, {uh:.2f}] (N <= {math.exp(uh):.2e}): sd(c) = {sd_window(Ffull["hp"], math.log(100), uh, K, shape, trend):.5f}', flush=True)
