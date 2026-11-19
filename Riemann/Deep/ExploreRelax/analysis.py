"""explore51 analysis: (a) ceiling formula table, (b) Raymond's comparison path vs exact rows, (c) band rates,
(d) single-probe visibility of the periodic and the irregular rows against zeta's phase-uniform penalty."""
import sys, math, json
import numpy as np, mpmath as mp
sys.path.insert(0, 'Riemann/Deep/ExploreRelax')
from row_exact import (RHO, S, LOGT, Y0, t_row, coeffs, evolve, roots, hit_time, insert_row, env_cue,
                       env_lattice)
out = {}
PR = math.pi * RHO
B0 = 0.87338196661320752                   # strip of the DBNB2 path in damping.txt
Y0B = 2 * B0 - 1
print(f"rho={RHO:.9f} pi*rho=L/4={PR:.9f} L={LOGT:.9f}")
# (a) ceilings
print("(a) t_row(L, y0) = (16/L^2) log cosh(L y0/4) at T0 = X/2:")
for name, y0 in [('classical beta<1', 1.0), ('7/8', 0.75), ('b0=0.87338196661', Y0B), ('theta141', Y0)]:
    tb = t_row(RHO, y0)
    print(f"   {name:18s} y0={y0:.12f}  t_row={tb:.9f}  4y0/L-16log2/L^2={4*y0/LOGT-16*math.log(2)/LOGT**2:.9f}")
    out['ceiling_' + name] = tb
for target in [0.09, 0.08, 0.05, 0.02, 0.01]:
    lc = lambda u: u + math.log1p(math.exp(-2 * u)) - math.log(2)
    f = lambda L: (16 / L ** 2) * lc(L * Y0 / 4) - target
    lo, hi = 1.0, 1e5
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
    print(f"   theta141 ceiling = {target}: needs L = log(T0/2pi) = {hi:.3f}, T0 = 10^{(hi + math.log(2*math.pi))/math.log(10):.2f}")
# (b) Raymond path (DBNB2 damping.txt, b0 run) vs rows
rows = [l.split() for l in open('Riemann/Deep/DBNB2/evidence/damping.txt') if l.strip() and l.split()[0].isdigit()]
tq = [(float(r[1]), float(r[3])) for r in rows]
def y_row(t, y0):
    v = math.cosh(PR * y0) * math.exp(-PR ** 2 * t)
    return math.acosh(v) / PR if v > 1 else 0.0
N = 400; P = N * S
envs = {'cue1': env_cue(N, 1), 'cue2': env_cue(N, 2), 'lattice': env_lattice(N)}
zz = json.load(open('Riemann/Deep/ExploreRelax/logs/zeta_unfolded.json')); envs['zeta'] = np.array(zz['x'])
cfg = {'cue1 row': (envs['cue1'], 200), 'cue2 row': (envs['cue2'], 200), 'zeta row': (envs['zeta'], 200),
       'cue1 one pair': (envs['cue1'], 1), 'lattice one pair': (envs['lattice'], 1)}
curves = {}
for nm, (env, k) in cfg.items():
    zs, _ = insert_row(env, k, Y0B)
    m, c = coeffs(zs, P)
    curves[nm] = [float(max(abs(roots(m, evolve(m, c, P, t), P).imag)) ** 2) for t, _ in tq]
    curves[nm + ' hit'] = hit_time(zs, P)
print("(b) q(t): Raymond (DBNB2 path, b0) | periodic row | irregular rows | single pair")
print("   t        q_chain     q_row      cue1row    cue2row    zetarow    cue1pair   latpair")
for i, (t, q) in enumerate(tq):
    print(f"   {t:.5f}  {q:.7f}  {y_row(t, Y0B)**2:.7f}  " + '  '.join(f"{curves[nm][i]:.7f}" for nm in cfg))
print("   hit times: periodic row %.6f, " % t_row(RHO, Y0B) + ', '.join(f"{nm} {curves[nm+' hit']:.6f}" for nm in cfg))
print("   chain endpoint (DBNB2): 0.107559925317")
out['curves'] = curves; out['tq'] = tq
# (c) band rates (DBNAudit table, 7/8 run) vs periodic-row rate r_row(q) = (2 pi rho y coth(pi rho y) - 2)/q
def r_row(q):
    y = math.sqrt(q); x = PR * y
    return (2 * x / math.tanh(x) - 2) / q
bands = [(0.5625, 0.5, 12.56), (0.5, 0.4, 13.22), (0.4, 0.3, 14.16), (0.3, 0.2, 15.42), (0.2, 0.1, 17.38),
         (0.1, 0.05, 19.40), (0.05, 0.01, 21.81), (0.01, 1e-4, 22.68)]
print("(c) band: chain rate (max in band, 7/8 run) vs periodic-row rate at the band's upper q (its minimum)")
for qa, qb, rc in bands:
    print(f"   q {qa:.4g}-{qb:.4g}: chain {rc:.2f}  row {r_row(qa):.2f}  ratio {rc/r_row(qa):.3f}")
# (d) single-probe visibility at t = 0: row deficit vs zeta's phase-uniform penalty (1/2)(-zeta'/zeta((1+h)/2))/(pi rho)
mp.mp.dps = 20
def pen(sig):
    return float(-mp.zeta(sig, derivative=1) / mp.zeta(sig)) / 2
print("(d) single probe at height h (t=0, rows at y0): relative deficit of min_x J vs zeta's penalty/mean")
def Jper(xs, zs, P, h):
    """sum over zeros (periodic) of (h - b)/((x - a)^2 + (h - b)^2), closed form per zero."""
    tot = np.zeros_like(xs)
    for z in zs:
        c = h - z.imag
        tot += (np.pi / P) * np.sinh(2 * np.pi * c / P) / (np.cosh(2 * np.pi * c / P) - np.cos(2 * np.pi * (xs - z.real) / P))
    return tot
xs = np.linspace(0, P, 20001)
zs_cue, _ = insert_row(envs['cue1'], 200, Y0B)
zs_zeta, _ = insert_row(envs['zeta'], 200, Y0B)
real_cue = envs['cue1'].astype(complex)
print("   h      sigma   pen/mean   periodic row   cue1 row   zeta row   cue1 real")
vis = []
for h in [1.05, 1.1, 1.2, 1.36, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0]:
    sig = (1 + h) / 2
    pm = pen(sig) / PR
    dper = 1 - (math.tanh(PR * (h - Y0B) / 2) + math.tanh(PR * (h + Y0B) / 2)) / 2
    d = []
    for zs in [zs_cue, zs_zeta, real_cue]:
        J = Jper(xs, zs, P, h)
        d.append(1 - J.min() / PR)
    vis.append((h, pm, dper, *d))
    print(f"   {h:4.2f}  {sig:.3f}  {pm:.3e}   {dper:.3e}      {d[0]:.3e}  {d[1]:.3e}  {d[2]:.3e}")
out['visibility'] = vis
json.dump(out, open('Riemann/Deep/ExploreRelax/logs/analysis.json', 'w'), indent=1)
print('done')
