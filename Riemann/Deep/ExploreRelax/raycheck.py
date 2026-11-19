"""explore51: test model configurations against Raymond's ACTUAL source inequalities (Zenodo 23752103,
verification/arcs/ACTUAL_67_ROW_ARC_INPUTS.json: per row i a time interval, signed translated nodes (height h,
shift v, weight w), an optional jet overlay, and the certified affine source S_i(t) = a_i + b_i (t - start_i)).
For a configuration evolved by the exact periodic heat flow, J_i(x, t) = sum_l w_l G_t(x + v_l, h_l)
(+ jets), G_t(x, h) = -Im F_t'/F_t(x + i h) = sum over zeros of (h - b)/((x - a)^2 + (h - b)^2).
Report min over x and t of J_i - S_i (the configuration passes Raymond's inputs iff >= 0 everywhere) and the
value at the top zero.  Usage: python3 raycheck.py <rows67.json>"""
import sys, math, json
from fractions import Fraction as Fr
import numpy as np
sys.path.insert(0, 'Riemann/Deep/ExploreRelax')
from row_exact import RHO, S, Y0, coeffs, evolve, roots, hit_time, env_cue, insert_row
PR = math.pi * RHO; A = 2 * S; M = 200; P = M * A; MG = 2 ** 13
rows = []
for r in json.load(open(sys.argv[1])):
    f = lambda s: float(Fr(s))
    o = r['joint_jet_overlay']
    rows.append(dict(t0=f(r['time'][0]), t1=f(r['time'][1]), a=f(r['source_lower_intercept_at_row_start']),
                     b=f(r['source_lower_slope']), nodes=[(f(n['height']), f(n['shift']), f(n['weight'])) for n in r['nodes']],
                     ov=None if o is None else (f(o['alpha']), [(f(c['p']), [f(x) for x in c['b']]) for c in o['centres']])))
xs = np.arange(MG) * P / MG
def G(m, ct, h, v):
    """-Im F'/F(x + v + i h) on the grid xs."""
    damp = np.exp(-2 * np.pi * m * h / P) * np.exp(2j * np.pi * m * v / P)
    full = np.zeros(MG, dtype=complex); fd = np.zeros(MG, dtype=complex)
    full[m % MG] = ct * damp; fd[m % MG] = ct * damp * (2j * np.pi * m / P)
    return -np.imag(np.fft.ifft(fd) / np.fft.ifft(full))
def jets(zr, p, bs, x):
    """sum_m 2 b_m p^{m-1} (1/m!) d^m/dp^m G(x, p), via -Im[(-i)^m sum_j (z - z_j)^{-(m+1)}], images |n| <= 40."""
    z = x[:, None] + 1j * p
    tot = np.zeros(len(x))
    for m, b in enumerate(bs):
        s = np.zeros(len(x), dtype=complex)
        for n in range(-40, 41):
            s += np.sum((z - (zr[None, :] + n * P)) ** (-(m + 1)), axis=1)
        tot += 2 * b * p ** (m - 1) * (-np.imag(((-1j) ** m) * s)) if m > 0 else 2 * b * p ** (-1) * (-np.imag(s))
    return tot
def check(name, zs):
    th = hit_time(zs, P)
    m, c = coeffs(zs, P)
    worst, worst_top, info = 1e9, 1e9, None
    for r in rows:
        if r['t0'] >= th: break
        for t in (r['t0'], (r['t0'] + r['t1']) / 2, min(r['t1'], th) - 1e-9):
            ct = evolve(m, c, P, t)
            zr = roots(m, ct, P)
            itop = int(np.argmax(zr.imag)); xtop = zr[itop].real % P
            Sv = r['a'] + r['b'] * (t - r['t0'])
            for sgn in (1, -1):
                J = np.zeros(MG)
                for h, v, w in r['nodes']:
                    J += w * G(m, ct, h, sgn * v)
                if r['ov'] is not None:
                    alpha, cs = r['ov']; J *= alpha
                    sub = np.arange(0, MG, 16)
                    jj = np.zeros(len(sub))
                    for p, bs in cs:
                        jj += jets(zr, p, bs, xs[sub])
                    J = J[sub] + jj; xg = xs[sub]
                else:
                    xg = xs
                mm = float(J.min() - Sv)
                jt = float(np.interp(xtop, xg, J, period=P) - Sv)
                if mm < worst: worst, info = mm, (round(t, 5), sgn, round(float(xg[np.argmin(J)]), 2), round(Sv, 4))
                worst_top = min(worst_top, jt)
    print(f"{name:28s} hit {th:.6f}  min_(x,t) J-S = {worst:+.4f} (t, sign, x, S = {info})  at top zero: {worst_top:+.4f}", flush=True)
    return th, worst, worst_top
def rigid(seed, lc):
    env = env_cue(2 * M, seed); mids = (env[0::2] + env[1::2]) / 2
    base = (np.arange(M) + 0.5) * A
    if lc is None:
        pos = mids
    else:
        u = mids - base; u -= u.mean(); U = np.fft.rfft(u); k = np.arange(len(U))
        U[~((k > 0) & (M / np.maximum(k, 1) < lc))] = 0
        pos = base + np.fft.irfft(U, n=M)
    return np.concatenate([pos + 1j * Y0, pos - 1j * Y0])
out = {}
base = (np.arange(M) + 0.5) * A
# validation of G against the closed-form periodic Poisson sum (CUE row, t = 0)
zv = rigid(1, None); mv, cv = coeffs(zv, P)
for h, v in [(1.36, 0.0), (2.5, 0.7), (6.0, 1.95)]:
    cf = np.zeros(MG)
    for z in zv:
        cc = h - z.imag
        cf += (np.pi / P) * np.sinh(2 * np.pi * cc / P) / (np.cosh(2 * np.pi * cc / P) - np.cos(2 * np.pi * (xs + v - z.real) / P))
    print(f"validation h={h} v={v}: max |G - closed form| = {np.max(np.abs(G(mv, cv, h, v) - cf)):.2e} (mean {cf.mean():.4f})")
zr0 = roots(mv, cv, P)
jt = jets(zr0, 2.5, [1.0], xs[:64]); print(f"validation jets m=0: max |2G/p - jets| = {np.max(np.abs(2 * G(mv, cv, 2.5, 0.0)[:64] / 2.5 - jt)):.2e}")
out['periodic row'] = check('periodic row', np.concatenate([base + 1j * Y0, base - 1j * Y0]))
for seed, lc in [(1, 8), (1, 12), (2, 16), (1, None)]:
    out[f'cue{seed} row lc={lc}'] = check(f'cue{seed} row lc={lc}', rigid(seed, lc))
zz = json.load(open('Riemann/Deep/ExploreRelax/logs/zeta_unfolded.json'))
zs, _ = insert_row(np.array(zz['x']), 200, Y0); out['zeta-stat row'] = check('zeta-stat row', zs)
zs, _ = insert_row(env_cue(2 * M, 1), 1, Y0); out['one pair in CUE'] = check('one pair in CUE', zs)
json.dump(out, open('Riemann/Deep/ExploreRelax/logs/raycheck.json', 'w'), indent=1)
print('done')
