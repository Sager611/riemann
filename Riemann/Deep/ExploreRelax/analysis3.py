"""explore51 analysis 3: single-probe consistency along the whole trajectory, J_t(x,h) = -Im F_t'/F_t(x+ih)
computed from the evolved Fourier coefficients (exact), for the periodic row and rigidified CUE rows (cutoffs lc)."""
import sys, math, json
import numpy as np
sys.path.insert(0, 'Riemann/Deep/ExploreRelax')
from row_exact import RHO, S, LOGT, Y0, coeffs, evolve, roots, hit_time, env_cue
PR = math.pi * RHO; A = 2 * S
N0 = int(math.exp(LOGT / 2)) + 1
sv = np.ones(N0 + 1, dtype=bool); sv[:2] = False
for i in range(2, int(N0 ** 0.5) + 1):
    if sv[i]: sv[i * i::i] = False
lp = np.log(np.nonzero(sv)[0].astype(float))
def pen_t(h, t):
    """as in analysis2.py: primes p <= N0, exponent (1+h)/2 + (t/4) log(x/(4 pi p))."""
    e = (1 + h) / 2 + (t / 4) * (LOGT - lp)
    return 0.5 * float(np.sum(lp * np.exp(-e * lp)))
M = 200; P = M * A; MG = 2 ** 13
def Jmin(m, ct, h):
    """min over x of -Im F'/F(x + ih) on MG grid points."""
    full = np.zeros(MG, dtype=complex); fd = np.zeros(MG, dtype=complex)
    damp = np.exp(-2 * np.pi * m * h / P)
    full[m % MG] = ct * damp
    fd[m % MG] = ct * damp * (2j * np.pi * m / P)
    F = np.fft.ifft(full); Fd = np.fft.ifft(fd)
    return float(np.min(-np.imag(Fd / F)))
def config(seed, lc):
    if seed is None:
        pos = (np.arange(M) + 0.5) * A
    else:
        env = env_cue(2 * M, seed); mids = (env[0::2] + env[1::2]) / 2
        base = (np.arange(M) + 0.5) * A
        u = mids - base; u -= u.mean(); U = np.fft.rfft(u); k = np.arange(len(U))
        U[~((k > 0) & (M / np.maximum(k, 1) < lc))] = 0
        pos = base + np.fft.irfft(U, n=M)
    return np.concatenate([pos + 1j * Y0, pos - 1j * Y0])
out = {}
for seed, lc in [(None, None), (1, 8), (2, 8), (3, 8), (1, 12), (2, 12), (3, 12), (1, 16), (2, 16), (1, 24)]:
    zs = config(seed, lc)
    th = hit_time(zs, P)
    m, c = coeffs(zs, P)
    worst = (1e9, None, None)
    for t in np.linspace(0, th, 9)[:-1]:
        ct = evolve(m, c, P, t)
        ymax = float(max(abs(roots(m, ct, P).imag)))
        hmin = max(2 * (1.02 - t * LOGT / 8) - 1, ymax + 0.03)
        for h in np.concatenate([np.linspace(hmin, hmin + 1.0, 11), [hmin + 1.5, hmin + 2.5, 4.0, 6.0]]):
            d = 1 - Jmin(m, ct, h) / PR
            r = pen_t(h, t) / PR / d if d > 0 else float('inf')
            if r < worst[0]:
                worst = (r, float(t), float(h))
    print(f"seed {seed} lc {lc}: hit {th:.6f}; min over (t,h) of pen/deficit = {worst[0]:.2f} at t={worst[1]:.4f}, h={worst[2]:.3f}", flush=True)
    out[f'{seed}_{lc}'] = (th, worst)
json.dump(out, open('Riemann/Deep/ExploreRelax/logs/analysis3.json', 'w'), indent=1)
print('done')
