# explore29: zeros of zeta as a 1-D quantizer (Zador / Fejes Toth efficiency, Lloyd-Max necessary conditions).
# usage: uv run --no-project --with numpy --with scipy python3 quant_stats.py <dl dir>
import sys, gzip, math, time
import numpy as np
from scipy.optimize import brentq

DL = sys.argv[1]
TWO_PI = 2*math.pi

def theta(t):
    t = np.asarray(t, dtype=float)
    return t/2*np.log(t/TWO_PI) - t/2 - math.pi/8 + 1/(48*t) + 7/(5760*t**3)

def load_first():
    with gzip.open(DL + '/zeros1.gz', 'rt') as f:
        return np.array([float(x) for x in f.read().split()])

def load_high(name):
    vals = []
    for line in open(DL + '/' + name):
        s = line.strip()
        try:
            vals.append(float(s))
        except ValueError:
            pass
    return np.array(vals)

BASES = {'zeros3': 267653395647.0, 'zeros4': 1.44176897509546973e20, 'zeros5': 1.3709199099319953e21}

def peaks(g):
    # strict interior peaks g[j-1] < g[j] > g[j+1]
    return int(np.sum((g[1:-1] > g[:-2]) & (g[1:-1] > g[2:])))

def tilt(a, b):
    # beta with int_{-a}^{b} x e^{beta x} dx = 0 (exponential density on the cell [-a, b] with centroid 0)
    if abs(a - b) < 1e-15*(a+b):
        return 0.0
    def F(beta):
        if abs(beta) < 1e-9:
            return (b*b - a*a)/2 + beta*(b**3 + a**3)/3
        # int x e^{bx} = e^{bx}(x/beta - 1/beta^2)
        return (math.exp(beta*b)*(b/beta - 1/beta**2) - math.exp(-beta*a)*(-a/beta - 1/beta**2))
    lo, hi = -1.0, 1.0
    while F(lo)*F(hi) > 0:
        lo *= 2; hi *= 2
        if hi > 1e6:
            return float('nan')
    return brentq(F, lo, hi, xtol=1e-14, rtol=1e-12)

def stats(name, gam, delta, raw_gaps, meangap_local):
    d = delta
    N = len(d)
    k1 = np.mean(d**2)/np.mean(d)**2
    k2 = np.mean(d**3)/np.mean(d)**3
    k3 = np.mean(d**4)/np.mean(d)**4
    var = np.var(d/np.mean(d))
    c1 = np.corrcoef(d[:-1], d[1:])[0, 1]
    P = peaks(raw_gaps)
    Pn = peaks(d)
    # Lloyd residual for p = zero density (slowly varying): cell midpoint offset (g_n - g_{n-1})/4, normalised
    res = (d[1:] - d[:-1])/4.0
    # exponential tilt needed per interior cell, in units of 1/mean gap
    m = min(N-1, 20000)
    betas = np.array([tilt(d[i]/2, d[i+1]/2) for i in range(m)])
    sgn = np.sign(betas[np.isfinite(betas)])
    signchg = np.mean(sgn[1:] != sgn[:-1])
    tv = np.nansum(np.abs(betas)*(d[:m] + d[1:m+1])/2)/m   # total variation of log p per cell
    print(f"{name}: N={N+1} zeros, mean normalised gap={np.mean(d):.6f}")
    print(f"  Zador/Fejes Toth inefficiency kappa_r=<d^(r+1)>/<d>^(r+1): r=1 {k1:.5f}  r=2 {k2:.5f} ({10*math.log10(k2):.3f} dB)  r=3 {k3:.5f}")
    print(f"  var={var:.5f}  corr(d_n,d_n+1)={c1:+.4f}  max d={d.max()/np.mean(d):.4f}  min d={d.min()/np.mean(d):.5f}  <log d>={np.mean(np.log(d/np.mean(d))):+.5f}")
    print(f"  strict gap peaks: raw {P} ({P/(N-2):.4f} per interior gap), normalised {Pn} ({Pn/(N-2):.4f})")
    print(f"  Lloyd residual for the zero density: rms {np.sqrt(np.mean(res**2)):.4f} mean gaps; tilt |beta| median {np.nanmedian(np.abs(betas)):.4f}, "
          f"mean {np.nanmean(np.abs(betas)):.4f} per mean gap; sign changes of tilt {signchg:.4f}; TV(log p) per cell {tv:.4f}")
    return dict(N=N+1, k1=k1, k2=k2, k3=k3, var=var, c1=c1, peakfrac=P/(N-2))

t0 = time.time()
g1 = load_first()
print(f"loaded first {len(g1)} zeros, last {g1[-1]:.6f} [{time.time()-t0:.1f}s]")

# smallest case
print("smallest case: first 8 zeros and gaps")
gs = np.diff(g1[:8])
print("  zeros:", ' '.join(f"{x:.6f}" for x in g1[:8]))
print("  gaps :", ' '.join(f"{x:.6f}" for x in gs))
for j in range(1, len(gs)-1):
    if gs[j] > gs[j-1] and gs[j] > gs[j+1]:
        print(f"  strict peak at gap gamma_{j+2}-gamma_{j+1} = {gs[j]:.6f} (> {gs[j-1]:.6f}, > {gs[j+1]:.6f}): zeros {j}..{j+3} need a valley of p in V_{j+1} u V_{j+2}")

x = theta(g1)/math.pi + 1
d_all = np.diff(x)
gaps_all = np.diff(g1)
for (a, b) in [(0, 1000), (1000, 10000), (10000, 100000), (0, 100000)]:
    stats(f"zeros {a+1}..{b}", g1[a:b], d_all[a:b-1], gaps_all[a:b-1], None)

for nm in ['zeros3', 'zeros4', 'zeros5']:
    v = load_high(nm)
    base = BASES[nm]
    gam_approx = base + v
    gaps = np.diff(v)
    dens = math.log(base/TWO_PI)/TWO_PI
    stats(f"{nm} (height {base:.4g})", gam_approx, gaps*dens, gaps, None)
print(f"done [{time.time()-t0:.1f}s]")
