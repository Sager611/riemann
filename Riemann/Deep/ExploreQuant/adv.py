# explore29 adversarial tests: does a quantizer criterion see a zero pair off the line?
# (A) zeta with one pair of consecutive zeros replaced by an off-line pair (t = g0 +- i delta), (B) Davenport-Heilbronn.
# Criteria: merged on-line gap (Lloyd-Max/Zador see only gaps), log-concavity (Fleischer's hypothesis = Laguerre),
# Lloyd-Max 2-point uniqueness of |Xi|^alpha on the merged interval and on its sub-intervals, reverse water-filling R(D).
# usage: uv run --no-project --with mpmath --with numpy python3 adv.py <dl dir>
import sys, gzip, math, time, json
import numpy as np, mpmath as mp
mp.mp.dps = 25
DL = sys.argv[1]
t0 = time.time()
TWO_PI = 2*math.pi

with gzip.open(DL + '/zeros1.gz', 'rt') as f:
    G = np.array([float(x) for x in f.read().split()])
thG = G/2*np.log(G/TWO_PI) - G/2 - math.pi/8 + 1/(48*G)
DZ = np.diff(thG/math.pi)                     # unfolded gaps, first 1e5 zeros
V = []
for line in open(DL + '/zeros5'):
    try: V.append(float(line.strip()))
    except ValueError: pass
V = np.array(V); D22 = np.diff(V)*math.log(1.3709199099319953e21/TWO_PI)/TWO_PI

# sine-kernel tail P(s > x) = -E0'(x)
xg, wg = np.polynomial.legendre.leggauss(48)
def E0(s):
    x = (xg + 1)*s/2; w = wg*s/2
    Dm = x[:, None] - x[None, :]
    with np.errstate(invalid='ignore', divide='ignore'):
        K = np.where(np.abs(Dm) < 1e-15, 1.0, np.sin(math.pi*Dm)/(math.pi*Dm))
    A = np.sqrt(w)[:, None]*K*np.sqrt(w)[None, :]
    return float(np.linalg.det(np.eye(len(x)) - A))
def gue_tail(x, h=1e-4):
    return -(E0(x + h) - E0(x - h))/(2*h)

def tails(x):
    return (f"P_zeta(d>={x:.3f}): first 1e5 {np.mean(DZ >= x):.2e}, near 1e22 {np.mean(D22 >= x):.2e}; sine kernel {gue_tail(x):.2e}")

def stationary_2pt(tg, p):
    # zeros of F(m) = (c_L(m) + c_R(m))/2 - m on the grid (2-point Lloyd-Max stationary thresholds)
    dt = np.diff(tg)
    cp = np.concatenate(([0.0], np.cumsum((p[1:] + p[:-1])/2*dt)))
    xp = tg*p
    cxp = np.concatenate(([0.0], np.cumsum((xp[1:] + xp[:-1])/2*dt)))
    Tp, Txp = cp[-1], cxp[-1]
    i = np.arange(2, len(tg) - 2)
    cL = cxp[i]/cp[i]; cR = (Txp - cxp[i])/(Tp - cp[i])
    F = (cL + cR)/2 - tg[i]
    sc = np.nonzero(np.sign(F[1:]) != np.sign(F[:-1]))[0]
    return len(sc), [float(tg[i][k]) for k in sc]

def analyse(label, tg, h, gbar, sub=True, want_rd=None):
    lh = np.log(h)
    d2 = (lh[2:] - 2*lh[1:-1] + lh[:-2])/(tg[1] - tg[0])**2*gbar**2
    inner = slice(len(d2)//20, len(d2) - len(d2)//20)
    nmax = int(np.sum((h[1:-1] > h[:-2]) & (h[1:-1] > h[2:])))
    out = {"label": label, "max (log h)'' * gbar^2": float(np.max(d2[inner])), "local maxima of h": nmax}
    for al in (1.0, 2.0):
        p = h**al
        n, ms = stationary_2pt(tg, p)
        rec = {"stationary 2-pt on whole interval": n}
        if sub:
            nonu = 0; tot = 0
            L = len(tg)
            for frac in np.linspace(0.5, 1.0, 26):
                for side in (0, 1):
                    k = int(frac*(L - 1))
                    if side == 0: tt, pp = tg[:k+1], p[:k+1]
                    else: tt, pp = tg[L-1-k:], p[L-1-k:]
                    if len(tt) < 50: continue
                    nn, _ = stationary_2pt(tt, pp); tot += 1
                    if nn > 1: nonu += 1
            rec["sub-intervals with >1 stationary 2-pt quantizer"] = f"{nonu}/{tot}"
        out[f"alpha={al:g}"] = rec
    if want_rd is not None:
        S = h**2; S = S/np.mean(want_rd**2)
        rd = []
        for th in [1e-1, 1e-2, 1e-3, 1e-4, 1e-6, 1e-8]:
            Dd = float(np.mean(np.minimum(th, S))); R = float(np.mean(np.maximum(0.0, 0.5*np.log(S/th))))
            rd.append((th, Dd, R))
        out["R(D) water-filling (theta, D, R)"] = [(f"{a:.0e}", f"{b:.3e}", f"{c:.4f}") for a, b, c in rd]
    return out

# ---------- (A) zeta
gi = np.diff(G[:10000]); ia = int(np.argmin(gi))            # Lehmer pair (0-based index ia, ia+1)
cands = [i for i in range(6600, 6700) if 0.97 < DZ[i] < 1.03]
ig = cands[0]
print(f"Lehmer pair: zeros #{ia+1},#{ia+2} at {G[ia]:.6f}, {G[ia+1]:.6f} (gap {gi[ia]:.6f}, normalised {DZ[ia]:.4f}); generic pair #{ig+1},#{ig+2} (normalised gap {DZ[ig]:.4f})")
def Fx(t):   # Xi = -F Z, F = (t^2+1/4)/2 * pi^(-1/4) |Gamma(1/4+it/2)|
    return (t*t + 0.25)/2*mp.pi**(-0.25)*abs(mp.gamma(mp.mpf(0.25) + 1j*t/2))
results = []
for (name, i0) in [("Lehmer pair", ia), ("generic pair", ig)]:
    za, zb = mp.zetazero(i0 + 1).imag, mp.zetazero(i0 + 2).imag
    zl, zr = mp.zetazero(i0).imag, mp.zetazero(i0 + 3).imag
    g0 = (za + zb)/2
    gbar = TWO_PI/math.log(float(g0)/TWO_PI)
    M = 2401
    tg = np.linspace(float(zl), float(zr), M)
    tg[0] += 1e-9; tg[-1] -= 1e-9
    absZ = np.array([float(abs(mp.siegelz(mp.mpf(t)))) for t in tg])
    Fv = np.array([float(Fx(mp.mpf(t))) for t in tg])
    lin = np.array([float(abs((mp.mpf(t) - za)*(mp.mpf(t) - zb))) for t in tg])
    merged = float(zr - zl)/gbar
    print(f"{name}: zeros {float(zl):.6f} {float(za):.6f} {float(zb):.6f} {float(zr):.6f}; mean gap {gbar:.4f}; merged gap {merged:.4f} mean gaps; {tails(merged)}  [{time.time()-t0:.0f}s]")
    htrue = Fv*absZ
    # true Xi on each consecutive-zero sub-interval
    for (lo, hi) in [(zl, za), (za, zb), (zb, zr)]:
        msk = (tg > float(lo)) & (tg < float(hi))
        if msk.sum() > 60:
            r = analyse(f"{name}: true Xi on ({float(lo):.4f},{float(hi):.4f})", tg[msk], htrue[msk], gbar, sub=False)
            print("  ", json.dumps(r))
    for dd in [0.02, 0.1, 0.2, 0.35, 0.5]:
        delta = dd*gbar
        quad = (tg - float(g0))**2 + delta**2
        h = Fv*absZ*quad/lin
        r = analyse(f"{name}: off-line pair t = g0 +- i {dd:g} gbar", tg, h, gbar, sub=True, want_rd=htrue)
        print("  ", json.dumps(r)); results.append(r)
    S0 = htrue**2/np.mean(htrue**2)
    rd0 = [(f"{th:.0e}", f"{float(np.mean(np.minimum(th, S0))):.3e}", f"{float(np.mean(np.maximum(0.0, 0.5*np.log(S0/th)))):.4f}") for th in [1e-1, 1e-2, 1e-3, 1e-4, 1e-6, 1e-8]]
    print(f"   true Xi R(D) water-filling (theta, D, R): {rd0}")
print(f"[{time.time()-t0:.0f}s]")

# ---------- (B) Davenport-Heilbronn
dh = json.load(open(sys.argv[2]))
on = np.array(dh['online'])
uDH = np.array([ (t/TWO_PI)*math.log(5*t/(TWO_PI*math.e)) for t in on])
dDH = np.diff(uDH)
print(f"DH: {len(on)} on-line zeros up to {dh['TMAX']}, unfolded mean gap {np.mean(dDH[20:]):.4f}")
merged_list = []
for (beta, gam) in dh['offline']:
    k = int(np.searchsorted(on, gam))
    gl, gr = on[k-1], on[k]
    gbar = TWO_PI/math.log(5*gam/TWO_PI)
    mg = (gr - gl)/gbar
    merged_list.append(mg)
    print(f"  off-line zero {beta:.6f}+{gam:.6f}i: flanking on-line zeros {gl:.6f}, {gr:.6f}, merged gap {mg:.4f} mean gaps; {tails(mg)}")
others = np.array([d for j, d in enumerate(np.diff(on)/np.array([TWO_PI/math.log(5*t/TWO_PI) for t in on[:-1]])) if on[j] > 20 and not any(on[j] < g < on[j+1] for _, g in dh['offline'])])
print(f"  DH ordinary on-line gaps (t>20): n={len(others)}, max {others.max():.4f}, 99th pct {np.percentile(others, 99):.4f}; merged gaps min {min(merged_list):.4f}, max {max(merged_list):.4f}")
# Laguerre + Lloyd-Max for DH around the first off-line zero
mp.mp.dps = 20
chi = [0, 1, 1j, -1j, -1]; chib = [0, 1, -1j, 1j, -1]
kap = (mp.sqrt(10 - 2*mp.sqrt(5)) - 2)/(mp.sqrt(5) - 1)
def hDH(t):
    s = mp.mpf(0.5) + 1j*t
    f = ((1 - 1j*kap)*mp.dirichlet(s, chi) + (1 + 1j*kap)*mp.dirichlet(s, chib))/2
    return float((mp.mpf(5)/mp.pi)**0.25*abs(mp.gamma(mp.mpf(0.75) + 1j*t/2))*abs(f))
beta, gam = dh['offline'][0]
k = int(np.searchsorted(on, gam)); gl, gr = on[k-1], on[k]
gbar = TWO_PI/math.log(5*gam/TWO_PI)
tg = np.linspace(gl + 1e-7, gr - 1e-7, 901)
hv = np.array([hDH(mp.mpf(t)) for t in tg])
print("  ", json.dumps(analyse(f"DH merged interval ({gl:.4f},{gr:.4f})", tg, hv, gbar, sub=True)))
gl2, gr2 = on[k-3], on[k-2]
tg2 = np.linspace(gl2 + 1e-7, gr2 - 1e-7, 401)
hv2 = np.array([hDH(mp.mpf(t)) for t in tg2])
print("  ", json.dumps(analyse(f"DH ordinary interval ({gl2:.4f},{gr2:.4f})", tg2, hv2, gbar, sub=False)))
print(f"done [{time.time()-t0:.0f}s]")
