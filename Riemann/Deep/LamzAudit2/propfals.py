# propfals.py -- lamz2audit: falsification tests of Lamzouri 2609.02882 Prop 2.1 (Pro:FourierHilbert, TeX v2 l.303-324)
# on random and adversarial finite conjugation-invariant multisets; tightest cases rechecked in Arb (python-flint).
# Run: uv run --no-project --with numpy --with python-flint python3 propfals.py > logs/propfals.log
import numpy as np, time, sys
from flint import acb, arb, ctx
rng = np.random.default_rng(20261111)
t0 = time.time()
def P(*a): print(*a); sys.stdout.flush()
def sinc(z):
    z = np.asarray(z, dtype=complex); sm = np.abs(z) < 1e-6
    return np.where(sm, 1 - z*z/6 + z**4/120, np.sin(z)/np.where(sm, 1, z))
# kernels K = FT of f = eta^2 (f >= 0, even, int f = 1, supp in [-lam, lam]); Prop 2.1 needs only eta in L^2.
def Kcos(w):   # f = cos(w x)/N on [-1/2,1/2], 0 <= w < pi ; w = sqrt2: Montgomery-Taylor f_0 ; w = 0: indicator
    N = np.sinc(w/(2*np.pi)) if w > 0 else 1.0     # = 2 sin(w/2)/w
    return lambda xi: (sinc((w - 2*np.pi*xi)/2) + sinc((w + 2*np.pi*xi)/2))/(2*N)
def Kfejer(lam):  # f = (lam-|x|)/lam^2 on [-lam,lam]
    return lambda xi: sinc(np.pi*lam*xi)**2
def Kbox(lam):    # f = 1/(2 lam) on [-lam,lam]
    return lambda xi: sinc(2*np.pi*lam*xi)
KERNELS = {"MT f0 (w=sqrt2)": Kcos(np.sqrt(2)), "cos w=2.5": Kcos(2.5), "box lam=1/2": Kbox(0.5),
           "Fejer lam=1/2": Kfejer(0.5), "box lam=2": Kbox(2.0), "Fejer lam=3": Kfejer(3.0)}
# arb versions for recheck
def acb_sinc(z):
    return z.sinc() if hasattr(z, "sinc") else (z.sin()/z if not z.contains(0) else acb(1))
def Kcos_acb(w):
    w = arb(w) if not isinstance(w, arb) else w
    N = 2*(w/2).sin()/w if w != 0 else arb(1)
    pi = arb.pi()
    return lambda xi: (acb_sinc((w - 2*pi*xi)/2) + acb_sinc((w + 2*pi*xi)/2))/(2*N)
def Kfejer_acb(lam):
    lam = arb(lam); return lambda xi: acb_sinc(arb.pi()*lam*xi)**2
def Kbox_acb(lam):
    lam = arb(lam); return lambda xi: acb_sinc(2*arb.pi()*lam*xi)
KACB = {"MT f0 (w=sqrt2)": Kcos_acb(arb(2).sqrt()), "cos w=2.5": Kcos_acb(arb("2.5")), "box lam=1/2": Kbox_acb("0.5"),
        "Fejer lam=1/2": Kfejer_acb("0.5"), "box lam=2": Kbox_acb(2), "Fejer lam=3": Kfejer_acb(3)}
# configuration: reals (positions xr, mult mr), pairs (positions zc with Im>0, mult mc); multiset = xr^mr U zc^mc U conj(zc)^mc
def slacks(K, xr, mr, zc, mc):
    pts = np.concatenate([xr.astype(complex), zc, np.conj(zc)]); m = np.concatenate([mr, mc, mc]).astype(float)
    D = pts[:, None] - pts[None, :]
    S = np.sum(m[:, None]*m[None, :]*K(D)**2)
    Z = m.sum(); nsr = np.sum(mr == 1); distinct = len(pts)
    simple = np.sum(mr == 1) + 2*np.sum(mc == 1); real_mult = mr.sum()
    sor = real_mult + 2*np.sum(mc == 1)
    s1 = nsr - (2*Z - S.real); s2 = distinct - (1.5*Z - 0.5*S.real); s3 = simple + real_mult - (3*Z - S.real)
    A = S.real/Z; s4 = np.nan
    if 1 <= A < 2: s4 = sor - (5 + 2*np.sqrt(2) - 2*A)/(3 + 2*np.sqrt(2))*Z
    return np.array([s1/Z, s2/Z, s3/Z, s4/Z, S.real/Z - 1]), abs(S.imag)/max(1, abs(S.real))
def rand_conf(kind):
    n = rng.integers(0, 12); r = rng.integers(0, 4); k = rng.integers(0, 5)
    if n + r + k == 0: n = 1
    W = {"spaced": n + r + k + 1.0, "cluster": rng.uniform(0.05, 1.0), "wide": 5.0*(n + r + k + 1)}[kind]
    xr = np.sort(rng.uniform(0, W, n + r)); mr = np.concatenate([np.ones(n, int), rng.integers(2, 4, r)]); rng.shuffle(mr)
    ymax = rng.choice([1e-6, 1e-3, 0.05, 0.3, 1.0, 3.0])
    zc = rng.uniform(0, W, k) + 1j*rng.uniform(1e-9, ymax, k); mc = rng.choice([1, 1, 1, 2, 3], k)
    return xr, mr, zc, mc
names = ["(2.1) simple-real", "(2.2) distinct", "(2.3) simple+real", "(2.4) simple-or-real", "Rem2.2 S/|Z|-1"]
worst = {}
P("== random scan (float64): min slack/|Z| per inequality (negative = violation)")
for kn, K in KERNELS.items():
    mins = np.full(5, np.inf); arg = [None]*5; cnt = 0; imax = 0.0
    for kind in ["spaced", "cluster", "wide"]:
        for it in range(1500):
            c = rand_conf(kind); s, im = slacks(K, *c); cnt += 1; imax = max(imax, im)
            for i in range(5):
                if not np.isnan(s[i]) and s[i] < mins[i]: mins[i] = s[i]; arg[i] = c
    worst[kn] = arg
    P(" %-18s %d configs, max |Im S|/|S| = %.1e" % (kn, cnt, imax), " | ".join("%s %.3e" % (names[i].split()[0], mins[i]) for i in range(5)))
# adversarial: n simple real + k simple pairs, minimise S - (n + 4k) (= slack of (2.1)*|Z|) by random local search
P("== adversarial local search on (2.1): minimise S - (n+4k) over positions (n simple real, k simple pairs)")
adv = {}
for kn, K in KERNELS.items():
    best_overall = (np.inf, None)
    for trial in range(24):
        n = int(rng.integers(1, 9)); k = int(rng.integers(1, 4))
        xr = np.sort(rng.uniform(0, n + k, n)); zc = rng.uniform(0, n + k, k) + 1j*rng.uniform(0.01, 1.5, k)
        mr = np.ones(n, int); mc = np.ones(k, int)
        def obj(xr, zc):
            s, _ = slacks(K, xr, mr, zc, mc); return s[0]*(n + 2*k)
        cur = obj(xr, zc); step = 0.5
        for it in range(1500):
            xr2 = xr + rng.normal(0, step, n)*(rng.random(n) < 0.5)
            zc2 = zc + (rng.normal(0, step, k) + 1j*rng.normal(0, step, k))*(rng.random(k) < 0.5)
            zc2 = zc2.real + 1j*np.clip(np.abs(zc2.imag), 1e-9, 5)
            v = obj(xr2, zc2)
            if v < cur: xr, zc, cur = xr2, zc2, v
            if it % 300 == 299: step *= 0.5
        if cur < best_overall[0]: best_overall = (cur, (xr.copy(), mr.copy(), zc.copy(), mc.copy()))
    adv[kn] = best_overall
    xr, mr, zc, mc = best_overall[1]
    P(" %-18s min S-(n+4k) = %.6e  (n=%d, k=%d, min Im z = %.2e)" % (kn, best_overall[0], len(xr), len(zc), np.min(zc.imag)))
# Arb recheck of the tightest cases
P("== Arb recheck (256 bits) of the tightest configurations: slack balls")
ctx.prec = 256
def slacks_acb(Ka, xr, mr, zc, mc):
    pts = [acb(float(x)) for x in xr] + [acb(float(z.real), float(z.imag)) for z in zc] + [acb(float(z.real), -float(z.imag)) for z in zc]
    m = [int(v) for v in mr] + [int(v) for v in mc]*2
    S = acb(0)
    for i in range(len(pts)):
        for j in range(len(pts)):
            S += m[i]*m[j]*Ka(pts[i] - pts[j])**2
    Z = sum(m); nsr = int(np.sum(mr == 1)); distinct = len(pts)
    simple = int(np.sum(mr == 1)) + 2*int(np.sum(mc == 1)); real_mult = int(mr.sum())
    return [nsr - (2*Z - S.real), distinct - (arb(3)/2*Z - S.real/2), simple + real_mult - (3*Z - S.real), S.real - Z], S.imag
for kn in KERNELS:
    Ka = KACB[kn]
    cases = [("adv(2.1)", adv[kn][1])] + [("rand" + names[i].split()[0], worst[kn][i]) for i in range(4) if worst[kn][i] is not None]
    for lab, c in cases:
        sl, im = slacks_acb(Ka, *c)
        ok = all((x >= 0) or x.contains(0) for x in sl)
        P(" %-18s %-14s slacks(2.1),(2.2),(2.3),S-|Z|: %s  Im S: %s  %s" % (kn, lab, ", ".join(x.str(6) for x in sl), im.str(3), "OK" if ok else "VIOLATION"))
P("t = %.1f s" % (time.time() - t0))
