"""explore51 analysis 2: (i) the periodic row against zeta's phase-uniform penalty at t > 0;
(ii) 'rigidified' irregular rows: CUE pair positions with their long-wavelength displacement removed (zeta's zeros are
prime-rigid at large scales), hitting time vs single-probe visibility at t = 0."""
import sys, math, json
import numpy as np, mpmath as mp
sys.path.insert(0, 'Riemann/Deep/ExploreRelax')
from row_exact import RHO, S, LOGT, Y0, t_row, hit_time, env_cue
PR = math.pi * RHO
A = 2 * S                     # pair spacing of a full row
mp.mp.dps = 20
# primes up to N0 = sqrt(x/4pi) = e^{L/2}
N0 = int(math.exp(LOGT / 2)) + 1
sv = np.ones(N0 + 1, dtype=bool); sv[:2] = False
for i in range(2, int(N0 ** 0.5) + 1):
    if sv[i]: sv[i * i::i] = False
pr = np.nonzero(sv)[0].astype(float)
lp = np.log(pr)
def pen_t(h, t):
    """lower bound for zeta's phase-uniform penalty of a height-h probe of H_t: primes p <= N0 only, exponent
    (1+h)/2 + (t/4) log(x/(4 pi p)) (Polymath 15 b_n^t factor), Lambda(p) = log p."""
    e = (1 + h) / 2 + (t / 4) * (LOGT - lp)
    return 0.5 * float(np.sum(lp * np.exp(-e * lp)))
def y_row(t, y0=Y0):
    v = math.cosh(PR * y0) * math.exp(-PR ** 2 * t)
    return math.acosh(v) / PR if v > 1 else 0.0
def deficit_periodic(h, y):
    return 1 - (math.tanh(PR * (h - y) / 2) + math.tanh(PR * (h + y) / 2)) / 2
print(f"(i) periodic row vs zeta penalty (primes <= N0 = {N0}); probes need sigma_eff = (1+h)/2 + tL/8 >= 1.02")
worst = 1e9
for t in [0.0, 0.0125, 0.025, 0.05, 0.075, 0.09]:
    y = y_row(t)
    hmin = max(2 * (1.02 - t * LOGT / 8) - 1, y + 0.02)
    rows = []
    for h in np.concatenate([np.linspace(hmin, hmin + 1, 21), np.linspace(hmin + 1, 6.5, 12)[1:]]):
        d = deficit_periodic(h, y); p = pen_t(h, t) / PR
        rows.append((h, d, p, p / d if d > 0 else float('inf')))
    m = min(rows, key=lambda r: r[3])
    worst = min(worst, m[3])
    print(f"   t={t:.4f} y_row={y:.4f} h>={hmin:.3f}: min pen/deficit = {m[3]:.2f} at h={m[0]:.3f} (deficit {m[1]:.3e}, pen/mean {m[2]:.3e})")
print(f"   overall min ratio {worst:.2f} (>1: the periodic row satisfies every single-probe source bound)")
# (ii) rigidified CUE rows
def J_min(zs, P, h, xs):
    tot = np.zeros_like(xs)
    for z in zs:
        c = h - z.imag
        tot += (np.pi / P) * np.sinh(2 * np.pi * c / P) / (np.cosh(2 * np.pi * c / P) - np.cos(2 * np.pi * (xs - z.real) / P))
    return tot.min()
M = 200; P = M * A
xs = np.linspace(0, P, 8001)
hs = [1.05, 1.2, 1.36, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0]
print("(ii) row of 200 pairs at height y0, positions = CUE(400) pair midpoints with displacement modes of")
print("     wavelength > lc pair spacings removed; hit time and min over h of pen/deficit (t = 0)")
res = {}
for seed in [1, 2, 3]:
    env = env_cue(2 * M, seed)
    mids = (env[0::2] + env[1::2]) / 2
    u = mids - (np.arange(M) + 0.5) * A - np.mean(mids - (np.arange(M) + 0.5) * A)
    U = np.fft.rfft(u)
    for lc in [None, 8, 4, 3, 2]:
        Uc = U.copy()
        if lc is not None:
            k = np.arange(len(U)); keep = (k > 0) & (M / np.maximum(k, 1) < lc)
            Uc[~keep] = 0
        uc = np.fft.irfft(Uc, n=M)
        pos = (np.arange(M) + 0.5) * A + uc
        zs = np.concatenate([pos + 1j * Y0, pos - 1j * Y0])
        th = hit_time(zs, P)
        ratios = []
        for h in hs:
            d = 1 - J_min(zs, P, h, xs) / PR
            p = pen_t(h, 0.0) / PR
            ratios.append(p / d if d > 0 else float('inf'))
        mi = int(np.argmin(ratios))
        mingap = np.min(np.diff(np.sort(pos))) / A; maxgap = np.max(np.diff(np.sort(pos))) / A
        print(f"   seed {seed} lc={str(lc):4s}: hit {th:.6f}  min pen/deficit {ratios[mi]:.2f} at h={hs[mi]}  gaps/a in [{mingap:.2f},{maxgap:.2f}]", flush=True)
        res[f'{seed}_{lc}'] = (th, ratios)
json.dump(res, open('Riemann/Deep/ExploreRelax/logs/analysis2.json', 'w'), indent=1)
print('done')
