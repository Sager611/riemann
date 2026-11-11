# explore29 adversarial tests, zeta part (v2: log-space |Xi|, cached Z): does a quantizer criterion see a zero pair off the line?
# Model: Xi_delta(t) = Xi(t) * ((t-g0)^2 + delta^2) / ((t-ga)(t-gb)): consecutive zeros ga, gb replaced by the pair t = g0 +- i delta
# (s = 1/2 +- delta + i g0). Criteria: merged on-line gap; log-concavity of |Xi| (Fleischer's hypothesis = first Laguerre inequality);
# Lloyd-Max 2-point uniqueness for p = |Xi|^alpha on the merged interval and its sub-intervals; reverse water-filling R(D) of S = |Xi|^2.
# usage: uv run --no-project --with mpmath --with numpy python3 -u adv2.py <dl dir> <cache dir>
import sys, gzip, math, time, json, os
import numpy as np, mpmath as mp
mp.mp.dps = 25
DL, CACHE = sys.argv[1], sys.argv[2]
t0 = time.time()
TWO_PI = 2*math.pi
with gzip.open(DL + '/zeros1.gz', 'rt') as f:
    G = np.array([float(x) for x in f.read().split()])
thG = G/2*np.log(G/TWO_PI) - G/2 - math.pi/8 + 1/(48*G)
DZ = np.diff(thG/math.pi)

def stationary_2pt(tg, p):
    dt = np.diff(tg)
    cp = np.concatenate(([0.0], np.cumsum((p[1:] + p[:-1])/2*dt)))
    xp = (tg - tg[0])*p
    cxp = np.concatenate(([0.0], np.cumsum((xp[1:] + xp[:-1])/2*dt)))
    Tp, Txp = cp[-1], cxp[-1]
    i = np.arange(3, len(tg) - 3)
    cL = cxp[i]/cp[i]; cR = (Txp - cxp[i])/(Tp - cp[i])
    F = (cL + cR)/2 - (tg[i] - tg[0])
    ok = np.isfinite(F)
    F = F[ok]; ii = i[ok]
    sc = np.nonzero(np.sign(F[1:]) != np.sign(F[:-1]))[0]
    return len(sc), [round(float(tg[ii][k]), 4) for k in sc]

def analyse(label, tg, lh, gbar, sub=True):
    d2 = (lh[2:] - 2*lh[1:-1] + lh[:-2])/(tg[1] - tg[0])**2*gbar**2
    inner = slice(len(d2)//20, len(d2) - len(d2)//20)
    h = np.exp(lh - lh.max())
    nmax = int(np.sum((h[1:-1] > h[:-2]) & (h[1:-1] > h[2:])))
    conv = np.nonzero(d2[inner] > 0)[0]
    span = (float(tg[1:-1][inner][conv[0]]), float(tg[1:-1][inner][conv[-1]])) if len(conv) else None
    out = {"label": label, "max (log|Xi|)'' gbar^2": round(float(np.max(d2[inner])), 4),
           "log-convex span": None if span is None else [round(span[0], 4), round(span[1], 4)], "local maxima of |Xi|": nmax}
    for al in (1.0, 2.0):
        p = h**al
        n, ms = stationary_2pt(tg, p)
        rec = {"2-pt stationary (whole)": n, "thresholds": ms}
        if sub:
            nonu = 0; tot = 0; L = len(tg)
            for frac in np.linspace(0.3, 1.0, 36):
                for side in (0, 1):
                    k = int(frac*(L - 1))
                    tt, pp = (tg[:k+1], p[:k+1]) if side == 0 else (tg[L-1-k:], p[L-1-k:])
                    nn, _ = stationary_2pt(tt, pp); tot += 1
                    nonu += nn > 1
            rec["sub-intervals with >1 stationary"] = f"{nonu}/{tot}"
        out[f"alpha={al:g}"] = rec
    return out

def rd(lh, lh_ref):
    S = np.exp(2*(lh - lh_ref.max())); m = np.mean(np.exp(2*(lh_ref - lh_ref.max())))
    S = S/m
    res = []
    for th in [1e-1, 1e-2, 1e-3, 1e-4, 1e-6, 1e-8]:
        res.append([f"{th:.0e}", f"{float(np.mean(np.minimum(th, S))):.4e}", f"{float(np.mean(np.maximum(0.0, 0.5*np.log(S/th)))):.4f}"])
    return res

gi = np.diff(G[:10000]); ia = int(np.argmin(gi))
cands = [i for i in range(6600, 6700) if 0.97 < DZ[i] < 1.03]
ig = cands[0]
os.makedirs(CACHE, exist_ok=True)
for (name, i0) in [("Lehmer pair", ia), ("generic pair", ig)]:
    za, zb = mp.zetazero(i0 + 1).imag, mp.zetazero(i0 + 2).imag
    zl, zr = mp.zetazero(i0).imag, mp.zetazero(i0 + 3).imag
    g0 = (za + zb)/2
    gbar = TWO_PI/math.log(float(g0)/TWO_PI)
    M = 2401
    cf = f"{CACHE}/z_{i0}.npz"
    if os.path.exists(cf):
        c = np.load(cf); tg, lZ, lF = c['tg'], c['lZ'], c['lF']
    else:
        tg = np.linspace(float(zl), float(zr), M); tg[0] += 1e-9; tg[-1] -= 1e-9
        lZ = np.array([float(mp.log(abs(mp.siegelz(mp.mpf(t))))) for t in tg])
        lF = np.array([float(mp.log(mp.mpf(t)**2 + 0.25) + mp.re(mp.loggamma(mp.mpf(0.25) + 1j*mp.mpf(t)/2))) for t in tg])
        np.savez(cf, tg=tg, lZ=lZ, lF=lF)
    llin = np.array([float(mp.log(abs((mp.mpf(t) - za)*(mp.mpf(t) - zb)))) for t in tg])
    lt = lF + lZ
    merged = float(zr - zl)/gbar
    print(f"{name}: zeros {float(zl):.6f} {float(za):.6f} {float(zb):.6f} {float(zr):.6f}; mean gap {gbar:.4f}; merged gap {merged:.4f} mean gaps  [{time.time()-t0:.0f}s]")
    for (lo, hi) in [(zl, za), (za, zb), (zb, zr)]:
        msk = (tg > float(lo) + 1e-6) & (tg < float(hi) - 1e-6)
        if msk.sum() > 60:
            print("  ", json.dumps(analyse(f"{name}: true Xi on ({float(lo):.4f},{float(hi):.4f})", tg[msk], lt[msk], gbar, sub=False)))
    print(f"   true Xi R(D) (theta, D, R): {rd(lt, lt)}")
    for dd in [0.02, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0]:
        delta = dd*gbar
        lh = lt + np.log((tg - float(g0))**2 + delta**2) - llin
        r = analyse(f"{name}: pair t = g0 +- i {dd:g} gbar", tg, lh, gbar, sub=True)
        r["R(D)"] = rd(lh, lt)
        print("  ", json.dumps(r))
print(f"done [{time.time()-t0:.0f}s]")
