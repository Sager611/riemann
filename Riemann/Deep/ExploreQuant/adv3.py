# explore29 adv3: Lloyd-Max uniqueness for k = 2..4 points (Lloyd iteration from many starts) on the merged interval of the
# off-line-pair model (cached |Z| grids from adv2.py) and on DH's merged interval at 85.70; compared with log-concavity.
# usage: uv run --no-project --with mpmath --with numpy python3 -u adv3.py <dl dir> <cache dir> <dh json>
import sys, gzip, math, time, json, os
import numpy as np, mpmath as mp
DL, CACHE, DHJ = sys.argv[1], sys.argv[2], sys.argv[3]
t0 = time.time(); TWO_PI = 2*math.pi
rng = np.random.default_rng(7)

def lloyd_limits(tg, p, k, starts=150, iters=4000):
    x = tg - tg[0]; L = x[-1]
    dt = np.diff(x)
    cp = np.concatenate(([0.0], np.cumsum((p[1:] + p[:-1])/2*dt)))
    xp = x*p
    cxp = np.concatenate(([0.0], np.cumsum((xp[1:] + xp[:-1])/2*dt)))
    def cum(F, b):  # F at real positions b by linear interpolation
        return np.interp(b, x, F)
    lims = []
    for s in range(starts):
        y = np.sort(rng.uniform(0, L, k))
        for it in range(iters):
            m = np.concatenate(([0.0], (y[1:] + y[:-1])/2, [L]))
            P = cum(cp, m[1:]) - cum(cp, m[:-1]); X = cum(cxp, m[1:]) - cum(cxp, m[:-1])
            ok = P > 1e-300
            ynew = np.where(ok, X/np.where(ok, P, 1.0), (m[1:] + m[:-1])/2)
            if np.max(np.abs(ynew - y)) < 1e-12*L:
                y = ynew; break
            y = ynew
        Dv = 0.0
        m = np.concatenate(([0.0], (y[1:] + y[:-1])/2, [L]))
        for i in range(k):
            msk = (x >= m[i]) & (x <= m[i+1])
            Dv += np.trapezoid(((x[msk] - y[i])**2)*p[msk], x[msk]) if msk.sum() > 1 else 0.0
        lims.append((tuple(np.round(y/L, 4)), Dv/cp[-1]/L**2))
    uniq = {}
    for y, Dv in lims:
        key = y
        if not any(max(abs(a - b) for a, b in zip(key, kk)) < 2e-3 for kk in uniq):
            uniq[key] = Dv
    return sorted(uniq.items(), key=lambda kv: kv[1])

def report(label, tg, lh):
    h = np.exp(lh - lh.max())
    d2 = (lh[2:] - 2*lh[1:-1] + lh[:-2])
    nmax = int(np.sum((h[1:-1] > h[:-2]) & (h[1:-1] > h[2:])))
    out = {"label": label, "log-concave": bool(np.all(d2[len(d2)//20:-len(d2)//20] <= 0)), "modes": nmax}
    for al in (1.0, 2.0):
        for k in (2, 3, 4):
            u = lloyd_limits(tg, h**al, k)
            out[f"a{al:g}k{k}"] = len(u)
    print(json.dumps(out), f"[{time.time()-t0:.0f}s]")

with gzip.open(DL + '/zeros1.gz', 'rt') as f:
    G = np.array([float(x) for x in f.read().split()])
mp.mp.dps = 25
for fn in sorted(os.listdir(CACHE)):
    if not fn.startswith('z_'): continue
    i0 = int(fn[2:-4]); c = np.load(f"{CACHE}/{fn}"); tg, lZ, lF = c['tg'], c['lZ'], c['lF']
    za, zb = mp.zetazero(i0 + 1).imag, mp.zetazero(i0 + 2).imag
    g0 = float((za + zb)/2); gbar = TWO_PI/math.log(g0/TWO_PI)
    llin = np.array([float(mp.log(abs((mp.mpf(t) - za)*(mp.mpf(t) - zb)))) for t in tg])
    lt = lF + lZ
    for (lo, hi) in [(tg[0], float(za)), (float(zb), tg[-1])]:
        msk = (tg > lo + 1e-6) & (tg < hi - 1e-6)
        report(f"zeta #{i0+1}: true Xi on ({lo:.4f},{hi:.4f})", tg[msk], lt[msk])
    for dd in [0.02, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0]:
        lh = lt + np.log((tg - g0)**2 + (dd*gbar)**2) - llin
        report(f"zeta #{i0+1},#{i0+2} -> pair g0 +- i {dd:g} gbar", tg, lh)

dh = json.load(open(DHJ)); on = np.array(dh['online'])
cf = f"{CACHE}/dh85.npz"
beta, gam = dh['offline'][0]; k = int(np.searchsorted(on, gam)); gl, gr = on[k-1], on[k]
if os.path.exists(cf):
    c = np.load(cf); tg, lh = c['tg'], c['lh']
else:
    mp.mp.dps = 20
    chi = [0, 1, 1j, -1j, -1]; chib = [0, 1, -1j, 1j, -1]
    kap = (mp.sqrt(10 - 2*mp.sqrt(5)) - 2)/(mp.sqrt(5) - 1)
    def lhDH(t):
        s = mp.mpf(0.5) + 1j*t
        f = ((1 - 1j*kap)*mp.dirichlet(s, chi) + (1 + 1j*kap)*mp.dirichlet(s, chib))/2
        return float(mp.log((mp.mpf(5)/mp.pi)**0.25*abs(mp.gamma(mp.mpf(0.75) + 1j*t/2))*abs(f)))
    tg = np.linspace(gl + 1e-7, gr - 1e-7, 901)
    lh = np.array([lhDH(mp.mpf(t)) for t in tg])
    np.savez(cf, tg=tg, lh=lh)
report(f"DH merged interval ({gl:.4f},{gr:.4f}) around 0.8085+85.6993i", tg, lh)
print(f"done [{time.time()-t0:.0f}s]")
