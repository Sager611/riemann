"""explore51: the rigidified rows that pass Raymond's actual sources (raycheck.py) against his certified
comparison path (DBNB2 damping.txt, b0 = 0.87338196661320752): q_config(t) <= q_chain(t) must hold."""
import sys, math
import numpy as np
sys.path.insert(0, 'Riemann/Deep/ExploreRelax')
from row_exact import S, coeffs, evolve, roots, hit_time, env_cue
B0 = 0.87338196661320752; Y0B = 2 * B0 - 1
A = 2 * S; M = 200; P = M * A
rows = [l.split() for l in open('Riemann/Deep/DBNB2/evidence/damping.txt') if l.strip() and l.split()[0].isdigit()]
tq = [(float(r[1]), float(r[3])) for r in rows]
def rigid(seed, lc):
    env = env_cue(2 * M, seed); mids = (env[0::2] + env[1::2]) / 2
    base = (np.arange(M) + 0.5) * A
    u = mids - base; u -= u.mean(); U = np.fft.rfft(u); k = np.arange(len(U))
    U[~((k > 0) & (M / np.maximum(k, 1) < lc))] = 0
    pos = base + np.fft.irfft(U, n=M)
    return np.concatenate([pos + 1j * Y0B, pos - 1j * Y0B])
for seed, lc in [(1, 8), (1, 12), (2, 16)]:
    zs = rigid(seed, lc)
    m, c = coeffs(zs, P)
    qs = [float(max(abs(roots(m, evolve(m, c, P, t), P).imag)) ** 2) for t, _ in tq]
    gap = min(qc - q for (t, qc), q in zip(tq, qs))
    print(f"cue{seed} lc={lc}: hit {hit_time(zs, P):.6f}; q(t) " + ' '.join(f'{q:.5f}' for q in qs)
          + f"; min (q_chain - q) = {gap:+.5f}")
print("chain q(t) " + ' '.join(f'{q:.5f}' for _, q in tq) + " at t = " + ' '.join(f'{t:g}' for t, _ in tq))
