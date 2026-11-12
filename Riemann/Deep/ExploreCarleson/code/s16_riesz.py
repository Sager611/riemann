# explore42 s16: Carleson / Shapiro-Shields test for Burnol's zero kernels (model-space form, explore38 main term):
#   Gram G(g,g') = (1 - e^{i(Phi(g)-Phi(g'))})/(i(g'-g)), G(g,g) = Phi'(g), Phi_N(t) = t log N + 2 theta(t) + 6 arctan 2t.
# (a) smallest case {rho1, conj}; (b) head (K pairs): two-sided estimate |c~|^2/lmax <= B <= |c~|^2/lmin vs B;
# (c) local blocks of 200 consecutive zeros at several heights: lmin, lmax of the normalized Gram vs min normalized gap
#     sigma = min Phi' dgamma/2pi and the box count nmax (zeros in a window 2pi/Phi'); zeta vs CUE vs Poisson (same density), DH (on-line zeros).
import sys, math, json, time, resource
import numpy as np
t0 = time.time()
rng = np.random.default_rng(42)
OUT = '/Users/adriansager/playground/riemann/.scratch/explore42/'
zc = np.load(OUT + 'zeros_cache.npz'); g = zc['g']
def theta(t):
    return t / 2 * np.log(t / (2 * np.pi)) - t / 2 - np.pi / 8 + 1 / (48 * t) + 7 / (5760 * t ** 3)
def dtheta(t):
    return 0.5 * np.log(t / (2 * np.pi)) - 1 / (48 * t ** 2) - 7 / (1920 * t ** 4)
def Phi(t, L, q=1.0):
    if q == 1.0:
        return t * L + 2 * theta(np.abs(t)) * np.sign(t) + 6 * np.arctan(2 * t)
    return t * L + np.sign(t) * (np.abs(t) * np.log(q * np.abs(t) / (2 * np.pi)) - np.abs(t))
def dPhi(t, L, q=1.0):
    if q == 1.0:
        return L + 2 * dtheta(np.abs(t)) + 3 / (0.25 + t ** 2)
    return L + np.log(q * np.abs(t) / (2 * np.pi))
def gram(x, L, q=1.0):
    P = Phi(x, L, q); dP = dPhi(x, L, q)
    d = x[None, :] - x[:, None]
    with np.errstate(divide='ignore', invalid='ignore'):
        G = (1 - np.exp(1j * (P[:, None] - P[None, :]))) / (1j * d)
    G[np.arange(len(x)), np.arange(len(x))] = dP
    return G, dP
def nstats(x, L, q=1.0):
    G, dP = gram(x, L, q)
    s = 1 / np.sqrt(dP)
    Gn = s[:, None] * G * s[None, :]
    ev = np.linalg.eigvalsh(Gn)
    xs = np.sort(x); gaps = np.diff(xs); mid = 0.5 * (xs[1:] + xs[:-1])
    sig = gaps * dPhi(mid, L, q) / (2 * np.pi)
    w = 2 * np.pi / dPhi(xs, L, q)
    cnt = np.searchsorted(xs, xs + w) - np.arange(len(xs))
    return ev[0], ev[-1], sig.min(), cnt.max()
Ns = [10, 1280, 10 ** 10]
print('(a) smallest case {rho1, conj(rho1)} and (b) head two-sided estimate')
for N in Ns:
    L = math.log(N)
    for K in [1, 50, 400]:
        x = np.concatenate([g[:K], -g[:K]])
        G, dP = gram(x, L)
        rho = 0.5 + 1j * x
        c = (rho - 1) / rho ** 2
        H = np.conj(G).T.conj()  # H_ij = (Y_j, Y_i): main-term Gram is Hermitian; use G^T = conj(G)
        H = G.T
        B = float(np.real(np.vdot(c, np.linalg.solve(H, c))))
        ct2 = float(np.sum(np.abs(c) ** 2 / dP))
        s = 1 / np.sqrt(dP); ev = np.linalg.eigvalsh(s[:, None] * H * s[None, :])
        print(f'N={N:.0e} K={K:4d}: lmin={ev[0]:.4f} lmax={ev[-1]:.4f} | |c~|^2 L={ct2*L:.6f}  B L={B*L:.6f}  bounds [{ct2/ev[-1]*L:.6f}, {ct2/ev[0]*L:.6f}]  B/|c~|^2={B/ct2:.5f}', flush=True)
print('(c) local blocks of 200 zeros: median over blocks/samples of lmin, lmax, sigma_min, nmax')
starts = [1000, 5000, 10000, 20000]
for N in Ns:
    L = math.log(N)
    for i0 in starts:
        x = g[i0:i0 + 200]
        a = nstats(x, L)
        T0 = x[0]; sp = 2 * np.pi / math.log(T0 / (2 * np.pi))
        cue = []; poi = []
        for rep in range(12):
            Z = (rng.standard_normal((200, 200)) + 1j * rng.standard_normal((200, 200))) / np.sqrt(2)
            Q, R = np.linalg.qr(Z); Q = Q * (np.diag(R) / np.abs(np.diag(R)))[None, :]
            ang = np.sort(np.angle(np.linalg.eigvals(Q)) % (2 * np.pi))
            cue.append(nstats(T0 + ang / (2 * np.pi) * 200 * sp, L))
            poi.append(nstats(T0 + np.sort(rng.uniform(0, 1, 200)) * 200 * sp, L))
        cue = np.median(np.array(cue), axis=0); poi = np.median(np.array(poi), axis=0)
        print(f'N={N:.0e} T~{T0:8.0f}: zeta lmin={a[0]:.4f} lmax={a[1]:.3f} sig={a[2]:.4f} nmax={a[3]} | CUE lmin={cue[0]:.4f} lmax={cue[1]:.3f} sig={cue[2]:.4f} nmax={cue[3]:.0f} | Poisson lmin={poi[0]:.2e} lmax={poi[1]:.3f} sig={poi[2]:.2e} nmax={poi[3]:.0f}', flush=True)
# whole-file trend: blocks of 200 across the file, record lmin and sigma_min at N = 1280
L = math.log(1280)
rows = []
for i0 in range(0, len(g) - 200, 200):
    a = nstats(g[i0:i0 + 200], L); rows.append((g[i0], *a))
rows = np.array(rows)
cc = np.corrcoef(np.log(rows[:, 1]), np.log(rows[:, 3]))[0, 1]
fit = np.polyfit(np.log(rows[:, 3]), np.log(rows[:, 1]), 1)
print(f'all {len(rows)} zeta blocks (N=1280): lmin range [{rows[:,1].min():.4f}, {rows[:,1].max():.4f}], sigma_min range [{rows[:,3].min():.4f}, {rows[:,3].max():.4f}]; corr(log lmin, log sigma_min)={cc:.3f}; lmin ~ sigma^{fit[0]:.2f} x {math.exp(fit[1]):.3f}; nmax max {rows[:,4].max():.0f}')
# DH on-line zeros (T < 420), q = 5, blocks of 100
d = json.load(open('/Users/adriansager/playground/riemann/Riemann/Deep/ExploreHelio/results/dh_zeros.json'))
on = np.array(sorted(d['online']))
for N in Ns:
    L = math.log(N)
    for i0 in [50, 150]:
        x = on[i0:i0 + 100]
        a = nstats(x, L, q=5.0)
        xz = g[(g > x[0]) & (g < x[-1])]
        b = nstats(xz, L) if len(xz) > 20 else (np.nan,) * 4
        print(f'N={N:.0e} DH on-line block [{x[0]:.0f},{x[-1]:.0f}] ({len(x)}): lmin={a[0]:.4f} lmax={a[1]:.3f} sig={a[2]:.4f} nmax={a[3]} | zeta same window ({len(xz)}): lmin={b[0]:.4f} lmax={b[1]:.3f} sig={b[2]:.4f}', flush=True)
print(f'{time.time()-t0:.1f}s maxrss {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6:.0f} MB')
