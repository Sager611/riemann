# N-E: M3 against a random planted configuration with density profile e^{-kappa a} (toy model of nc_toy.py, L = 30).
# Each background zero (rescaled zeta zeros, |gamma - t0| <= 75) becomes an off-line pair of depth (a+U)/L with
# probability e^{-kappa a} (a = 1..14, kappa = 8/7; total 0.47), plus forced pairs at a = 8, 11, 14.
# Moments int |M|^q over |t - t0| <= 60 (per unit length) vs background-only, the ratio to the Lemma Z majorant (A = 2),
# and the share of the moment carried by +-6/L windows around the forced deep pairs.
import numpy as np, sys, json, math
sys.path.insert(0, '/Users/adriansager/playground/riemann/Riemann/Deep/UPAudit')
from nc_toy import Toy, background
tiny = len(sys.argv) > 1
L, kap, t0 = 30.0, 8/7, 1000.0
m = Toy(L, du=0.02); rng = np.random.default_rng(7)
bg = background(L, t0, 90.0); g = bg.imag
pl = []
for gm in g[np.abs(g-t0) <= 75]:
    r = rng.random(); acc = 0.0
    for a in range(1, 15):
        acc += math.exp(-kap*a)
        if r < acc: pl.append((gm, (a+rng.random())/L)); break
for a, off in ((8, 20.3), (11, -33.7), (14, 47.1)): pl.append((t0+off, (a+0.5)/L))
onl = np.array([z for z in bg if not any(abs(z.imag-p[0]) < 1e-12 for p in pl)])
cfg = np.concatenate([onl] + [np.array([0.5+d+1j*gm, 0.5-d+1j*gm]) for gm, d in pl])
W = 4.0 if tiny else 60.0; step = (0.5 if tiny else 0.25)/L
ts = np.arange(t0-W, t0+W, step); qs = [1.0, 2.0, kap/m.th, 2.5, 3.0]
def win(zs, t, G=8.0): return zs[np.abs(zs.imag-t) <= G]
Mc = np.array([m.M(win(cfg, t), 0.5+1j*t) for t in ts]); Mb = np.array([m.M(win(bg, t), 0.5+1j*t) for t in ts])
def B2(zs, t, A=2.0):  # Lemma Z majorant at y = 0 with A = 2, constant 1
    kk = np.minimum(1.0, np.maximum(m.L*np.abs(zs.imag-t), 1e-300)**(-A)); d = zs.real-0.5
    return (1+np.sum(m.X**np.maximum(d, 0)*kk))*(1+kk.sum())**(m.k-1)
Bc = np.array([B2(win(cfg, t), t) for t in ts])
forced = {8: t0+20.3, 11: t0-33.7, 14: t0+47.1}
out = dict(n_bg=int(np.sum(np.abs(g-t0) <= 75)), n_planted=len(pl), deepest_random=sorted(round(d*L, 2) for _, d in pl[:-3])[-5:],
           max_ratio_to_majorant_A2=float(np.max(Mc/Bc)))
for q in qs:
    tot = float(np.sum(Mc**q)*step); row = dict(moment_per_unit=tot/(2*W), bg_per_unit=float(np.sum(Mb**q)*step)/(2*W),
           moment_over_majorant=tot/float(np.sum(Bc**q)*step))
    for a, gp in forced.items():
        sel = np.abs(ts-gp) <= 6/L
        if sel.any(): row[f'share_forced_a{a}'] = float(np.sum(Mc[sel]**q)*step)/tot
    out[f'q={q:.3f}'] = row
print(json.dumps(out, indent=0))
if not tiny: json.dump(out, open('/Users/adriansager/playground/riemann/Riemann/Deep/UPAudit/logs/ne_config.json', 'w'), indent=1)
