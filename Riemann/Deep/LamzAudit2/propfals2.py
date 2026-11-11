# propfals2.py -- lamz2audit: adversarial search on Prop 2.1 (2.1)-(2.3) of 2609.02882 with off-line pairs kept at Im z >= y0
# (no collapse to double real zeros) and with multiplicities; best cases rechecked in Arb.
# Run: uv run --no-project --with numpy --with python-flint python3 propfals2.py > logs/propfals2.log
import numpy as np, sys, time
from flint import acb, arb, ctx
rng = np.random.default_rng(7)
def sinc(z):
    z = np.asarray(z, dtype=complex); sm = np.abs(z) < 1e-6
    return np.where(sm, 1 - z*z/6 + z**4/120, np.sin(z)/np.where(sm, 1, z))
def Kcos(w):
    N = np.sinc(w/(2*np.pi)) if w > 0 else 1.0
    return lambda xi: (sinc((w - 2*np.pi*xi)/2) + sinc((w + 2*np.pi*xi)/2))/(2*N)
def acb_sinc(z): return z.sinc()
def Kcos_acb(w):
    N = 2*(w/2).sin()/w; pi = arb.pi()
    return lambda xi: (acb_sinc((w - 2*pi*xi)/2) + acb_sinc((w + 2*pi*xi)/2))/(2*N)
K = Kcos(np.sqrt(2)); ctx.prec = 256; Ka = Kcos_acb(arb(2).sqrt())
def S_of(pts, m, KK):
    D = pts[:, None] - pts[None, :]; return np.sum(m[:, None]*m[None, :]*KK(D)**2)
def lhs(xr, mr, zc, mc):
    Z = mr.sum() + 2*mc.sum(); nsr = np.sum(mr == 1); dist = len(xr) + 2*len(zc)
    simple = np.sum(mr == 1) + 2*np.sum(mc == 1); real_mult = mr.sum()
    return Z, nsr, dist, simple + real_mult
def slack(xr, mr, zc, mc, which):
    pts = np.concatenate([xr.astype(complex), zc, np.conj(zc)]); m = np.concatenate([mr, mc, mc]).astype(float)
    S = S_of(pts, m, K).real; Z, nsr, dist, sr = lhs(xr, mr, zc, mc)
    return [nsr - (2*Z - S), dist - (1.5*Z - 0.5*S), sr - (3*Z - S)][which]
t0 = time.time()
print("Kernel: Montgomery-Taylor f_0. Minimise slack (LHS - RHS, unnormalised) by random local search.")
for which, name in enumerate(["(2.1)", "(2.2)", "(2.3)"]):
    for y0 in [0.02, 0.1, 0.3]:
        best = (np.inf, None)
        for trial in range(16):
            n = int(rng.integers(0, 8)); r = int(rng.integers(0, 3)); k = int(rng.integers(1, 4))
            xr = rng.uniform(0, 4, n + r); mr = np.concatenate([np.ones(n, int), rng.integers(2, 4, r)])
            zc = rng.uniform(0, 4, k) + 1j*(y0 + rng.uniform(0, 1, k)); mc = rng.choice([1, 1, 2], k)
            cur = slack(xr, mr, zc, mc, which); step = 0.4
            for it in range(1200):
                xr2 = xr + rng.normal(0, step, len(xr))*(rng.random(len(xr)) < 0.5)
                zc2 = zc + (rng.normal(0, step, k) + 1j*rng.normal(0, step, k))*(rng.random(k) < 0.5)
                zc2 = zc2.real + 1j*np.maximum(np.abs(zc2.imag), y0)
                v = slack(xr2, mr, zc2, mc, which)
                if v < cur: xr, zc, cur = xr2, zc2, v
                if it % 300 == 299: step *= 0.5
            if cur < best[0]: best = (cur, (xr.copy(), mr.copy(), zc.copy(), mc.copy()))
        xr, mr, zc, mc = best[1]
        # Arb recheck
        pts = [acb(float(x)) for x in xr] + [acb(float(z.real), float(z.imag)) for z in zc] + [acb(float(z.real), -float(z.imag)) for z in zc]
        m = [int(v) for v in mr] + [int(v) for v in mc]*2
        S = acb(0)
        for i in range(len(pts)):
            for j in range(len(pts)): S += m[i]*m[j]*Ka(pts[i] - pts[j])**2
        Z, nsr, dist, sr = lhs(xr, mr, zc, mc)
        sl = [nsr - (2*int(Z) - S.real), dist - (arb(3)/2*int(Z) - S.real/2), int(sr) - (3*int(Z) - S.real)][which]
        print(" %s y0=%-5g min slack float %.6e | Arb %s | n_real=%d (mult %s) pairs=%d (mult %s, Im %s) | %s" % (
            name, y0, best[0], sl.str(8), len(xr), list(mr), len(zc), list(mc), np.round(zc.imag, 3), "OK" if (sl > 0 or sl.contains(0)) else "VIOLATION"))
        sys.stdout.flush()
print("t = %.1f s" % (time.time() - t0))
