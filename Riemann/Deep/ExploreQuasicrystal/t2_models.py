# t2_models.py -- explore5 T-B: finite-X Euler-product models with the true Gamma factor ("Kronecker-torus" models),
# zeta vs Davenport-Heilbronn: where do their zeros sit?
#   zeta model (analytic, = Gonek's "wrong AFE"):  G_X(s) = E(s) + E(1-s),  E(s) = Gamma_R(s) P_X(s),  P_X = exp(sum_{n<=X} Lambda(n) n^-s / log n)
#   DH model:  H_X(s) = Phi_X(s) + conj(Phi_X(1-conj s)),  Phi_X = (5/pi)^{s/2} Gamma((s+1)/2) [c1 P_X(s,chi) + c2 P_X(s,chibar)]
# On Re s = 1/2 both are 2 Re E (resp. 2 Re Phi_X); their real zeros coincide with Gonek's non-analytic model (0704.3448).
# Counts: real zeros on [T,T+W] (phase crossings / sign changes) vs all zeros in the strip 0.02<=sigma<=0.98 (argument principle).
# run: OMP_NUM_THREADS=1 uv run --no-project --with numpy --with scipy python3 t2_models.py
import numpy as np, math, sys, time
from scipy.special import loggamma
ZER = '/Users/adriansager/playground/riemann/.scratch/explore5/dl/zeros1'
zz = np.loadtxt(ZER)
def prime_powers(X):
    out = []
    for p in range(2, X + 1):
        if all(p % q for q in range(2, int(p ** 0.5) + 1)):
            q = p
            while q <= X: out.append((q, math.log(p), p)); q *= p
    return out
KAP = (math.sqrt(10 - 2 * math.sqrt(5)) - 2) / (math.sqrt(5) - 1)
C1 = (1 - 1j * KAP) / 2; C2 = (1 + 1j * KAP) / 2
CHI = {0: 0, 1: 1, 2: 1j, 3: -1j, 4: -1}
def dirpoly(s, pp, chi=None):
    # sum Lambda(n) chi(n) n^-s / log n
    acc = np.zeros_like(s, dtype=complex)
    for (n, L, p) in pp:
        c = 1.0 if chi is None else CHI[n % 5]
        if c == 0: continue
        acc += (c * L / math.log(n)) * np.exp(-s * math.log(n))
    return acc
def logGR(s):  return -(s / 2) * math.log(math.pi) + loggamma(s / 2)
def logGDH(s): return (s / 2) * math.log(5 / math.pi) + loggamma((s + 1) / 2)
def logE_zeta(s, pp):  return logGR(s) + dirpoly(s, pp)
def logEs_zeta(s, pp): return logGR(1 - s) + dirpoly(1 - s, pp)
def logsum2(a, b):
    # log(e^a + e^b) up to 2 pi i, returned as complex value of the scaled sum and the real scale
    m = np.maximum(a.real, b.real)
    return np.exp(a - m) + np.exp(b - m)
def G_scaled(s, model, pp):
    if model == 'zeta':
        return logsum2(logE_zeta(s, pp), logEs_zeta(s, pp))
    # DH: Phi(s) = exp(logGDH(s)) [C1 e^{D_chi(s)} + C2 e^{D_chibar(s)}];  Phi#(s) = exp(logGDH(1-s)) [C2 e^{D_chibar(1-s)} + C1 e^{D_chi(1-s)}]
    Dc = dirpoly(s, pp, 'c'); Dcb = np.conj(dirpoly(np.conj(s), pp, 'c'))   # D_chibar(s) = conj(D_chi(conj s))
    Dc1 = dirpoly(1 - s, pp, 'c'); Dcb1 = np.conj(dirpoly(np.conj(1 - s), pp, 'c'))
    a1 = logGDH(s) + Dc; a2 = logGDH(s) + Dcb; b1 = logGDH(1 - s) + Dcb1; b2 = logGDH(1 - s) + Dc1
    m = np.maximum.reduce([a1.real, a2.real, b1.real, b2.real])
    return C1 * np.exp(a1 - m) + C2 * np.exp(a2 - m) + C2 * np.exp(b1 - m) + C1 * np.exp(b2 - m)
def winding(model, pp, s1, s2, T1, T2, dt=None, ds=0.0005):
    if dt is None:
        rate = math.log(T2) + sum(L * n ** (-min(s1, 1 - s2)) for (n, L, p) in pp)   # bound on |d/dt arg| of the dominant half
        dt = min(0.004, 0.25 / rate)
    bot = s1 + (np.arange(0, s2 - s1, ds)) + 1j * T1
    rig = s2 + 1j * (T1 + np.arange(0, T2 - T1, dt))
    top = s2 - np.arange(0, s2 - s1, ds) + 1j * T2
    lef = s1 + 1j * (T2 - np.arange(0, T2 - T1, dt))
    path = np.concatenate([bot, rig, top, lef, [bot[0]]])
    g = G_scaled(path, model, pp)
    ph = np.unwrap(np.angle(g))
    jumps = np.max(np.abs(np.diff(ph)))
    return (ph[-1] - ph[0]) / (2 * math.pi), jumps, np.min(np.abs(g))
def real_zeros(model, pp, T1, T2, dt=0.002):
    t = T1 + np.arange(0, T2 - T1, dt); s = 0.5 + 1j * t
    if model == 'zeta':
        phi = (logGR(s)).imag - np.sum([ (L / math.log(n)) * n ** -0.5 * np.sin(t * math.log(n)) for (n, L, p) in pp ], axis=0) if pp else (logGR(s)).imag
        k = np.floor((phi - math.pi / 2) / math.pi)
        idx = np.nonzero(np.diff(k))[0]
        z = []
        for i in idx:
            lev = math.pi / 2 + math.pi * max(k[i], k[i + 1])
            if k[i + 1] < k[i]: lev = math.pi / 2 + math.pi * k[i]
            f0, f1 = phi[i] - lev, phi[i + 1] - lev
            z.append(t[i] + dt * f0 / (f0 - f1) if f0 != f1 else t[i])
        dphi = np.gradient(phi, dt)
        return np.array(z), dphi.min()
    g = G_scaled(s, model, pp).real
    idx = np.nonzero(np.sign(g[:-1]) != np.sign(g[1:]))[0]
    z = np.array([t[i] + dt * g[i] / (g[i] - g[i + 1]) for i in idx])
    return z, float('nan')
def unfold(z, model):
    s = 0.5 + 1j * z
    th = (logGR(s)).imag if model == 'zeta' else (logGDH(s)).imag
    return th / math.pi
def stats(z, model):
    u = unfold(z, model); g = np.diff(u)
    return g.var() / g.mean() ** 2 if len(g) > 2 else float('nan'), (g < 0.25 * g.mean()).mean() if len(g) > 2 else float('nan')
t0 = time.time()
W = 200.0
print("model  T      X    #pp  real  strip(arg.pr)  offline(0.02-0.98) offline(0.3-0.7)  min phi'   gapvar  small<.25  match_err(zeta)  [winding jumps, min|G|]")
for T in [1000.0, 10000.0, 50000.0, 70000.0]:
    zt = zz[(zz >= T) & (zz <= T + W)]
    vz, sz = stats(zt, 'zeta')
    print("zeta zeros on [%g,%g]: %d ; gap variance %.3f ; frac gaps<0.25 mean: %.3f" % (T, T + W, len(zt), vz, sz))
    for model in ['zeta', 'dh']:
        for X in [2, 5, 10, 20, 30, 50, 100, 300, 1000]:
            pp = prime_powers(X)
            z, mphi = real_zeros(model, pp, T, T + W)
            wnd, jmp, ming = winding(model, pp, 0.02, 0.98, T, T + W)
            wnd2, jmp2, _ = winding(model, pp, 0.30, 0.70, T, T + W)
            off = wnd - len(z)
            v, sm = stats(z, model)
            err = float('nan')
            if model == 'zeta' and len(z) > 0:
                d = np.abs(zt[:, None] - z[None, :]).min(axis=1)
                sp = 2 * math.pi / math.log(T / (2 * math.pi))
                err = np.median(d) / sp
            print("%-5s %7g %5d %4d %5d %10.3f %8.2f %8.2f %9.3f %8.3f %8.3f %10.3f   [%.2f %.2f, %.1e]" % (model, T, X, len(pp), len(z), wnd, off, wnd2 - len(z), mphi, v, sm, err, jmp, jmp2, ming))
            sys.stdout.flush()
print("elapsed %.0f s" % (time.time() - t0))
