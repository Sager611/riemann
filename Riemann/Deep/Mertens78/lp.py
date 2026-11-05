# lp.py -- LP search for signed certificates lambda (float); rigorous check in cert.py
import numpy as np, time, json, sys, math
from scipy.optimize import linprog
from mpmath import zeta as mz
from opt import best, X_of, ramare_cross
B = 0.875
def P(x, y): return x/(np.pi*(x*x + y*y))
def lp(s1, lines, G, margin=2e-3, ygrid=None, eta_grid=None, mabs=0.0, Ymax=3.0, dy=0.1, nfine=601):
    x1 = s1 - B
    if ygrid is None:
        ygrid = np.unique(np.round(np.concatenate([np.arange(-Ymax, Ymax + 1e-9, dy), np.arange(-0.6, 0.6001, 0.02)]), 6))
    if eta_grid is None:
        gfar = np.geomspace(40, 1e5, 400)
        eta_grid = np.unique(np.concatenate([np.linspace(-40, 40, 801), np.linspace(-Ymax-1, Ymax+1, int(400*(Ymax+1))+1), np.linspace(-0.6, 0.6, nfine), np.linspace(-0.05, 0.05, 501), gfar, -gfar]))
    cols = []; cp = []; cq = []; meta = []
    for s0 in lines:
        x0 = s0 - B; K0 = math.log(float(mz(s0))/float(mz(2*s0))); K1 = math.log(float(mz(s0)))
        cols.append(P(x0, eta_grid[:, None] - ygrid[None, :])); cp += [K0 + G]*len(ygrid); cq += [K1 - G]*len(ygrid)
        meta += [(s0, y) for y in ygrid]
    A = np.hstack(cols); n = A.shape[1]
    Aub = np.hstack([-A, A]); bub = -(1 + margin)*P(x1, eta_grid) - mabs*P(5.0, eta_grid)
    # (C2) with margin: sum lam x_j >= (1+margin) x1
    xs = np.array([s0 - B for s0, _ in meta])
    Aub = np.vstack([Aub, np.concatenate([-xs, xs])]); bub = np.concatenate([bub, [-(1 + margin)*x1]])
    c = np.array(cp + cq)
    t0 = time.time()
    res = linprog(c, A_ub=Aub, b_ub=bub, bounds=[(0, None)]*(2*n), method='highs')
    lam = res.x[:n] - res.x[n:]
    masses = [(meta[k][0], float(meta[k][1]), float(lam[k])) for k in range(n) if abs(lam[k]) > 1e-10]
    return masses, res.status, time.time() - t0
if __name__ == "__main__":
    cfg = json.loads(sys.argv[1])
    out = []
    for s1 in cfg["s1"]:
        for G in cfg["G"]:
            masses, st, dt = lp(s1, cfg["lines"], G)
            Lam = sum(l for _,_,l in masses); L1 = sum(abs(l) for _,_,l in masses)
            line = f"s1={s1} G={G} st={st} {dt:.0f}s nm={len(masses)} Lam={Lam:.4f} |lam|={L1:.3f}"
            bestq = {}
            for eta in [0.005, 0.01, 0.02, 0.03, 0.05, 0.1]:
                r = best(s1, eta, masses)
                if not r: continue
                for tgt in [0.95, 0.94, 0.93, 0.92, 0.91]:
                    X = X_of(tgt, r)
                    if tgt not in bestq or X < bestq[tgt][0]: bestq[tgt] = (X, eta, r['theta'], r['log10C'])
                cr = ramare_cross(r)
                if 'cross' not in bestq or (cr and cr < bestq['cross'][0]): bestq['cross'] = (cr, eta, r['theta'], r['log10C'])
            print(line, {k: (round(v[0],1) if v[0] else None, v[1], round(v[2],4), round(v[3],2)) for k, v in bestq.items()}, flush=True)
            out.append(dict(s1=s1, G=G, lines=cfg["lines"], masses=masses))
    json.dump(out, open(cfg.get("save", "lpout.json"), "w"))

