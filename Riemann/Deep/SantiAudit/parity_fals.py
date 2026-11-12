# Adversarial test of the finite inequalities used in Thm 5.1 on conjugation-invariant multisets:
#  (37) (Q-S)(N-O) >= 2(N-S)^2   [4, Prop 2.5];  Lamzouri S >= 2N-Q;  (33-weak) 3S >= 2N-Q+2O.
# K = FT of eta^2, eta^2 a probability density on (-lam/2, lam/2): indicator or cosine (optimal) density.
import numpy as np
rng = np.random.default_rng(12345)
lam = 0.527
def K_ind(xi):  # eta^2 = 1/lam on [-lam/2, lam/2]
    w = np.pi * lam * xi
    out = np.ones_like(w, dtype=complex)
    nz = np.abs(w) > 1e-12
    out[nz] = np.sin(w[nz]) / w[nz]
    return out
a = lam / 2; k = np.sqrt(2.0); Acos = k / (2 * np.sin(k * a))   # f(u) = Acos*cos(k u) on [-a,a], int f = 1
def K_cos(xi):
    # int_{-a}^{a} Acos cos(k u) e^{-2 pi i xi u} du = Acos * [sin((k-w)a)/(k-w) + sin((k+w)a)/(k+w)], w = 2 pi xi
    w = 2 * np.pi * xi
    def sc(b):
        out = np.full_like(b, a, dtype=complex); nz = np.abs(b) > 1e-12
        out[nz] = np.sin(b[nz] * a) / b[nz]; return out
    return Acos * (sc(k - w) + sc(k + w))
def stats(points, mults, K):
    z = np.array(points, dtype=complex); m = np.array(mults, dtype=float)
    D = z[:, None] - z[None, :]
    Q = np.sum(m[:, None] * m[None, :] * K(D) ** 2)
    N = m.sum(); real = np.abs(z.imag) < 1e-15
    S = np.sum(real & (m == 1)); O = np.sum(real & (m % 2 == 1))
    return Q, N, S, O
def random_set(spread):
    pts, ms = [], []
    nr = rng.integers(0, 8); nc = rng.integers(0, 4)
    for _ in range(nr):
        pts.append(rng.uniform(0, spread)); ms.append(int(rng.choice([1, 1, 1, 2, 3])))
    for _ in range(nc):
        x = rng.uniform(0, spread); y = rng.uniform(1e-4, 1.0) * rng.choice([1e-3, 1e-1, 1.0, 3.0])
        mm = int(rng.choice([1, 1, 2])); pts += [x + 1j * y, x - 1j * y]; ms += [mm, mm]
    if not pts: pts, ms = [0.0], [2]
    return pts, ms
worst = {}
for K, name in [(K_ind, "indicator"), (K_cos, "cosine")]:
    w37 = w_l = w33 = np.inf; maxim = 0.0
    for trial in range(60000):
        spread = rng.choice([0.5, 2.0, 5.0, 20.0])
        pts, ms = random_set(spread)
        Q, N, S, O = stats(pts, ms, K)
        maxim = max(maxim, abs(Q.imag)); Q = Q.real
        r37 = (Q - S) * (N - O) - 2 * (N - S) ** 2
        rl = S - (2 * N - Q); r33 = 3 * S - (2 * N - Q + 2 * O)
        # normalise by N^2 for comparability
        w37 = min(w37, r37 / N**2); w_l = min(w_l, rl / N); w33 = min(w33, r33 / N)
    worst[name] = (w37, w_l, w33, maxim)
    print(name, "min slack/N^2 (37):", w37, " min slack/N Lamzouri:", w_l, " min slack/N (33-weak):", w33, " max |Im Q|:", maxim)
# targeted family: one simple real point + a conjugate pair at distance a, height y (worst-case cross term)
for K, name in [(K_ind, "indicator"), (K_cos, "cosine")]:
    best = np.inf; arg = None
    for a0 in np.linspace(0.0, 6.0, 601):
        for y in [1e-3, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.2]:
            Q, N, S, O = stats([0.0, a0 + 1j * y, a0 - 1j * y], [1, 1, 1], K)
            r = (Q.real - S) * (N - O) - 2 * (N - S) ** 2
            if r < best: best, arg = r, (a0, y, Q.real)
    print(name, "targeted {x0, z, zbar}: min slack (37) =", best, "at (a, y, Q) =", arg)
