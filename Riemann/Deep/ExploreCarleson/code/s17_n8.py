# explore42 s17 (N8 bridge): Carleson separation constant of the shifted zeros z = rho + eps in Re s > 1/2.
# Identity (RH zeros, simple; unconditional for eps > B - 1/2 with (2 beta + 2 eps - 1)): B(s) = prod b_{rho+eps}(s) = xi(s - eps)/xi(s + eps), so
#   delta_n(eps) = prod_{m != n} |z_n - z_m|/|z_n + conj(z_m) - 1| = 2 eps |xi'(rho_n)| / |xi(rho_n + 2 eps)|.
# Test: truncated product over 22,491 zeros (+ conjugates, + density tail) vs mpmath; then delta at eps = a/log T (unfolded) for zeta, CUE, Poisson.
import math, time
import numpy as np
import mpmath as mp
mp.mp.dps = 30
t0 = time.time()
zc = np.load('/Users/adriansager/playground/riemann/.scratch/explore42/zeros_cache.npz'); g = zc['g']
allg = np.concatenate([g, -g])
def prod_delta(n, eps):
    d = g[n] - allg
    d = d[np.abs(d) > 1e-9]
    s = -0.5 * np.sum(np.log1p(4 * eps * eps / (d * d)))
    T = g[-1]
    s += -2 * eps * eps * (math.log(T / (2 * math.pi)) + 1) / (math.pi * T)   # density tail beyond T (both signs)
    return math.exp(s)
def xi(s):
    return 0.5 * s * (s - 1) * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)
def xid(rho):
    return 0.5 * rho * (rho - 1) * mp.pi ** (-rho / 2) * mp.gamma(rho / 2) * mp.zeta(rho, derivative=1)
print('identity check')
for eps in [0.4, 0.1]:
    for n in [0, 1, 4, 99]:
        rho = mp.mpc(0.5, g[n])
        rhs = float(2 * eps * abs(xid(rho)) / abs(xi(rho + 2 * eps)))
        lhs = prod_delta(n, eps)
        print(f'eps={eps} n={n+1} gamma={g[n]:.4f}: product={lhs:.8f} formula={rhs:.8f} rel.diff={lhs/rhs-1:+.2e}', flush=True)
# fixed eps: delta over all zeros (product formula, cheap per zero: use local window + smooth remainder)
def delta_all(eps, idx):
    out = []
    for n in idx:
        out.append(prod_delta(n, eps))
    return np.array(out)
idx = np.arange(0, len(g), 50)
for eps in [0.4, 0.1]:
    dl = delta_all(eps, idx)
    T = g[idx]
    print(f'fixed eps={eps}: delta at T~100: {dl[T<200].mean():.4f}, T~1e4: {dl[(T>9000)&(T<11000)].mean():.4f}, T~2e4: {dl[T>18000].mean():.4f}; inf over sample {dl.min():.2e}; mean delta*(T/2pi)^eps: {np.mean(dl*(T/(2*np.pi))**eps):.3f}')
# unfolded scale eps = a/log(T/2pi): delta is a functional of the local configuration; compare zeta / CUE / Poisson (blocks of 400 around T0)
rng = np.random.default_rng(7)
def delta_local(x, eps, tailfac):
    # x: sorted points (one block), delta for the central half using the block, far zeros by density approx
    out = []
    n = len(x)
    for i in range(n // 4, 3 * n // 4):
        d = x[i] - np.delete(x, i)
        out.append(-0.5 * np.sum(np.log1p(4 * eps * eps / (d * d))))
    return np.exp(np.array(out) + tailfac)
for a in [0.5, 1.0, 2.0]:
    res = []
    for i0 in [5000, 15000]:
        x = g[i0:i0 + 400]; T0 = x[200]; Lt = math.log(T0 / (2 * math.pi)); eps = a / Lt
        sp = 2 * math.pi / Lt
        # far-field: density Lt/2pi outside the block half-width W ~ 200 sp: -0.5 * 2 * int_W^inf (Lt/2pi) * 4 eps^2/d^2 dd
        W = 100 * sp
        tail = -0.5 * 2 * (Lt / (2 * math.pi)) * 4 * eps * eps / W
        dz = delta_local(x, eps, tail)
        cue = []; poi = []
        for rep in range(6):
            Z = (rng.standard_normal((400, 400)) + 1j * rng.standard_normal((400, 400))) / np.sqrt(2)
            Q, R = np.linalg.qr(Z); Q = Q * (np.diag(R) / np.abs(np.diag(R)))[None, :]
            ang = np.sort(np.angle(np.linalg.eigvals(Q)) % (2 * np.pi))
            cue.append(delta_local(T0 + ang / (2 * np.pi) * 400 * sp, eps, tail))
            poi.append(delta_local(T0 + np.sort(rng.uniform(0, 1, 400)) * 400 * sp, eps, tail))
        cue = np.concatenate(cue); poi = np.concatenate(poi)
        res.append((T0, np.median(dz), dz.min(), np.median(cue), np.quantile(cue, 0.005), np.median(poi), np.quantile(poi, 0.005)))
        print(f'a={a} T0={T0:.0f}: zeta median delta={np.median(dz):.4f} min={dz.min():.4f} | CUE median={np.median(cue):.4f} q0.5%={np.quantile(cue,0.005):.4f} | Poisson median={np.median(poi):.4f} q0.5%={np.quantile(poi,0.005):.2e} | e^-a={math.exp(-a):.4f}', flush=True)
print(f'{time.time()-t0:.1f}s')
