# validation tests for the dbn5 certificate machinery (small, quick)
import numpy as np, math, random, time
from fractions import Fraction as Fr
import lib5 as L, cert5 as C
from scipy.optimize import minimize
rng = np.random.default_rng(7)
# ---------- T1: interval arithmetic vs exact rationals ----------
bad = 0
for trial in range(2000):
    x = rng.normal(size=2)*10.0**rng.integers(-5, 5, size=2); y = rng.normal(size=2)*10.0**rng.integers(-5, 5, size=2)
    X = C.IA(min(x), max(x)); Y = C.IA(min(y), max(y))
    for op in ('+', '-', '*'):
        Z = X + Y if op == '+' else (X - Y if op == '-' else X*Y)
        ex = []
        for a in (min(x), max(x)):
            for b in (min(y), max(y)):
                A = Fr(a); B = Fr(b); ex.append(A + B if op == '+' else (A - B if op == '-' else A*B))
        if not (Fr(float(Z.lo)) <= min(ex) and Fr(float(Z.hi)) >= max(ex)): bad += 1
    s = abs(x[0]); Q = C.IA(s, s).sqrt()
    if not (Fr(float(Q.lo))**2 <= Fr(s) <= Fr(float(Q.hi))**2): bad += 1
v = rng.random(10000)*1e-3
if not Fr(C.sum_up(v)) >= sum(Fr(float(t)) for t in v): bad += 1
print('T1 interval arithmetic failures:', bad, flush=True)
# ---------- T3: toy paired certificate vs brute force over free label phases ----------
class Syn: pass
def toy(n, aw, seed):
    r = np.random.default_rng(seed)
    lam = np.sort(r.uniform(1.9, 14, n)); a = r.uniform(0.2, 1, n)*np.exp(-0.55*lam)*aw
    return a, lam
def synth(a, lam, ash):
    G = Syn(); n = len(a)
    G.M = a.copy(); G.a = a.copy(); G.lam_lo = lam.copy(); G.lam_hi = lam.copy(); G.nb = n
    G.Dg = float(np.sum(a*a*lam)); G.Un = 0.0
    I, J = np.triu_indices(n, 1)
    G.Slo = lam[I] + lam[J]; G.Shi = G.Slo.copy(); W = a[I]*a[J]; G.W0 = W
    th = ash/2*(lam[I] - lam[J])
    S = {'Clo': np.cos(ash*lam), 'Chi': np.cos(ash*lam), 'Wc': W*np.abs(np.cos(th)), 'Ws': W*np.abs(np.sin(th))}
    return G, S
def brute(a, lam, ash, starts=40, seed=0):
    r = np.random.default_rng(seed); best = -1e9
    def obj(ph):
        tot = 0.0; g = np.zeros_like(ph)
        for e in (+1, -1):
            T = a*np.exp(1j*(ph + e*ash/2*lam)); Q = 1 + T.sum(); V = (lam*T).sum()
            tot += 0.5*(V/Q).real; g += 0.5*np.real(1j*T*(lam/Q - V/Q**2))
        return -tot, -g
    for s in range(starts):
        x0 = r.uniform(-np.pi, np.pi, len(a))
        res = minimize(obj, x0, jac=True, method='L-BFGS-B')
        best = max(best, -res.fun)
    return best
worst = 1e9
for seed in range(2, 14):
    n = 8 + seed % 5; aw = [1.0, 2.0, 3.5][seed % 3]; ash = [0.5, 0.9, 1.3][seed % 3]
    a, lam = toy(n, aw, seed)
    if a.sum() > 0.35: a *= 0.35/a.sum()
    G, S = synth(a, lam, ash)
    Ct = C.one_point_fast(G, S); res = C.paired_fast(G, S, Ct, K=32)
    bf = brute(a, lam, ash)
    # one-point brute force
    def obj1(ph):
        T = a*np.exp(1j*ph); Q = 1 + T.sum(); V = (lam*T).sum(); return -(V/Q).real
    b1 = max(-minimize(obj1, np.random.default_rng(s).uniform(-np.pi, np.pi, len(a))).fun for s in range(30))
    worst = min(worst, res[0] - bf, Ct - b1)
    print(f'T3 seed {seed}: n {n} A {a.sum():.3f} B {(a*lam).sum():.3f}: one-point cert {Ct:.5f} >= brute {b1:.5f};  paired cert {res[0]:.5f} >= brute {bf:.5f}', flush=True)
print('T3 min margin (must be >= 0):', worst, flush=True)
# ---------- T4: Euler paired bound vs fine float grid ----------
from flint import arb  # pyright: ignore[reportMissingImports]
for (rho_, ash) in ((0.3, 0.5), (0.45, 1.1), (0.15, 0.8)):
    rmax = {p: arb(rho_)**(math.log(p)/math.log(2)) for p in (2, 3, 5)}
    Ep = C.euler_pair(rmax, rmax, Fr(ash))
    th = np.linspace(0, 2*np.pi, 200001); tot = 0.0
    for p in (2, 3, 5):
        r_ = float(rmax[p].mid()); d = ash/2*math.log(p)
        g = lambda u: (r_*np.cos(u) - r_*r_)/(1 - 2*r_*np.cos(u) + r_*r_)
        tot += math.log(p)/4*np.max(g(th + d) + g(th - d))
    print(f'T4 rho {rho_} a {ash}: bound {float(Ep.upper()):.10f} grid {tot:.10f} diff {float(Ep.upper()) - tot:.2e}', flush=True)
