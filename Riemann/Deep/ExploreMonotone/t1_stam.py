# T1: Perelman-type entropy -> finite free Fisher information / Stam / EPI (Garza-Vargas-Srivastava-Stier, arXiv 2602.15822)
# Model: truncated Xi = prod_{j<=N} (1 - z^2/x_j^2), x_j = 2*gamma_j (de Bruijn-Newman units, H_0(z) = Xi(z/2)/8).
# DBN heat flow dH/dt = -H'': zeros move dz_j/dt = 2 sum_k 1/(z_j - z_k).  GSS heat time = 2 * DBN time.
import json, os, time, sys
import numpy as np, mpmath as mp
t0 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
SCR = os.path.join(HERE, "..", "..", "..", ".scratch", "explore7")
N = 50
zf = os.path.join(HERE, "zeros50.json")
if os.path.exists(zf):
    gam = json.load(open(zf))
else:
    mp.mp.dps = 20
    gam = [float(mp.zetazero(k).imag) for k in range(1, N + 1)]
    json.dump(gam, open(zf, "w"))
xs = 2 * np.array(gam)
base = np.concatenate([-xs[::-1], xs])          # degree n = 2N, real
n = len(base)

def score(z):
    D = z[:, None] - z[None, :]
    np.fill_diagonal(D, np.inf)
    return (1.0 / D).sum(1)

def fisher(x):       # GSS normalisation
    s = score(x).real
    return np.mean((2.0 / (n - 1) * s) ** 2)

def chi(x):
    D = np.abs(x[:, None] - x[None, :]); iu = np.triu_indices(n, 1)
    return 2.0 / (n * (n - 1)) * np.sum(np.log(D[iu]))

He = np.polynomial.hermite_e.hermeroots([0] * n + [1])
PhiHe, chiHe = fisher(He), chi(He)
print(f"n={n}  Phi(He_n)={PhiHe:.6e} (theory 1/(n-1)={1/(n-1):.6e})  chi(He_n)={chiHe:.6f}")

def rk4(z, t, t_end, cfl=0.002, ystop=None):
    # integrate dz/dt = 2*score(z) from t to t_end (either direction); stop near collisions
    sgn = 1.0 if t_end > t else -1.0
    hist = []
    while (t_end - t) * sgn > 1e-15:
        if np.iscomplexobj(z):
            D = np.abs(z[:, None] - z[None, :]); np.fill_diagonal(D, np.inf); g = D.min()
        else:
            zz = np.sort(z); g = np.min(np.diff(zz))
        dt = sgn * min(cfl * g * g, abs(t_end - t))
        if g < 1e-7:
            return z, t, hist, "collision"
        k1 = 2 * score(z); k2 = 2 * score(z + dt / 2 * k1); k3 = 2 * score(z + dt / 2 * k2); k4 = 2 * score(z + dt * k3)
        z = z + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4); t += dt
        if ystop is not None:
            ymax = np.max(np.abs(z.imag))
            hist.append((t, ymax, np.sum(z.imag[z.imag > 0] ** 2), np.sum(z.imag[z.imag > 0])))
            if ymax < ystop:
                return z, t, hist, "landed"
    return z, t, hist, "done"

# ---------- (a) true model: backward (forward-heat) collision time vs Stam/EPI lower bounds ----------
Phi0, chi0 = fisher(base), chi(base)
gaps = np.diff(base); gmin = gaps.min(); imin = int(np.argmin(gaps))
z, tc, _, st = rk4(base.copy(), 0.0, -5.0, cfl=0.001)
zz = np.sort(z); j = int(np.argmin(np.diff(zz)))
stam_lb = -1.0 / (2 * (n - 1) * Phi0)                                  # Lambda_p >= this (DBN time)
epi_lb = -np.exp(2 * chi0) / (2 * np.exp(2 * chiHe))                    # from N(p0) >= N(p_-tau) + N(h_{2tau})
print(f"(a) truncated Xi, N={N}: Phi_n(p0)={Phi0:.4e}, chi_n(p0)={chi0:.4f}, min gap {gmin:.4f} at x={base[imin]:.3f}")
print(f"    true first backward collision: t_c = {tc:.6f} ({st}) at x = {zz[j]:.3f}; isolated-pair estimate -g^2/8 = {-gmin**2/8:.6f}")
print(f"    Stam lower bound Lambda_p >= {stam_lb:.4e} ; EPI lower bound Lambda_p >= {epi_lb:.4e}")

# ---------- (b) Stam / EPI along the forward (de Bruijn) flow, true model ----------
z = base.copy(); t = 0.0; rows = []
for t1 in [0.02, 0.05, 0.108, 0.2]:
    z, t, _, st = rk4(z, t, t1)
    P1, c1 = fisher(z), chi(z)
    stam_slack = 1 / P1 - (1 / Phi0 + 2 * (n - 1) * t1)
    epi_slack = np.exp(2 * c1) - (np.exp(2 * chi0) + 2 * t1 * np.exp(2 * chiHe))
    rows.append((t1, P1, stam_slack * P1, epi_slack / np.exp(2 * c1)))
    print(f"(b) t={t1:5.3f}: Phi={P1:.6e} (monotone dec: {P1 < Phi0}); Stam rel. slack {stam_slack*P1:.3e}; EPI rel. slack {epi_slack/np.exp(2*c1):.3e}")

# ---------- (c) planted off-line pair: detection power of Stam/EPI from data at t1 > landing ----------
for (k, Y) in [(24, 0.30), (24, 0.75), (10, 0.75)]:
    a, b = xs[k], xs[k + 1]; m = 0.5 * (a + b)
    pos = [v for v in xs if v not in (a, b)]
    zc = np.array([complex(v) for v in pos] + [complex(-v) for v in pos] +
                  [complex(m, Y), complex(m, -Y), complex(-m, Y), complex(-m, -Y)])
    zc, tl, hist, st = rk4(zc, 0.0, 2.0, cfl=0.002, ystop=2e-4)
    V = np.array([h[2] for h in hist]); S = np.array([h[3] for h in hist])
    mono = bool(np.all(np.diff(V) <= 1e-12) and np.all(np.diff(S) <= 1e-12))
    tland = tl + (hist[-1][1] ** 2) / 2.0
    # replace each landed pair by its real splitting and continue on the real side
    zr = zc.real.copy()
    order = np.argsort(np.abs(zc.imag))[::-1][:4]
    xr = np.sort(np.array([v for i, v in enumerate(zr) if i not in order] + [m, m, -m, -m], dtype=float))
    t1 = max(0.108, tland + 0.02)
    dt0 = 1e-6
    for sgnm in (m, -m):
        idx = np.where(np.isclose(xr, sgnm))[0]
        xr[idx[0]] -= np.sqrt(2 * dt0); xr[idx[1]] += np.sqrt(2 * dt0)
    zr2, _, _, _ = rk4(xr, tland + dt0, t1)
    P1 = fisher(zr2); c1 = chi(zr2)
    thrPhi = 1.0 / (2 * (n - 1) * t1)            # Stam detects non-real-rootedness at t=0 iff Phi(p_t1) >= this
    thrN = 2 * t1 * np.exp(2 * chiHe)              # EPI detects iff N(p_t1) <= this
    print(f"(c) planted pair at m={m:.2f}, Y={Y}: V,S monotone={mono}; landing t*={tland:.5f} (de Bruijn Y^2/2={Y*Y/2:.5f}); at t1={t1:.3f}: Phi={P1:.3e} vs Stam threshold {thrPhi:.3e} (ratio {P1/thrPhi:.2e}); N={np.exp(2*c1):.3e} vs EPI threshold {thrN:.3e} (ratio {np.exp(2*c1)/thrN:.2e})")
print(f"time {time.time()-t0:.1f}s")
