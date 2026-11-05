#!/usr/bin/env python3
# densityb2 explore2.py: numeric cross-check that the half-plane inputs do not change the zero-density
# exponent A(sigma) below sigma* = 0.856092 (and do change it above, reproducing DensityB).
# Own code. Reuses our DensityB pipeline (../DensityB/repro/explore.py: literature large-value estimates,
# Halasz-Montgomery with best beta, subdivision, raising to powers, Corollary zero-large-cor-0) and adds:
#   - fix: Ivic's zeta-from-exp bound is a MAX of its two terms (DensityB's explore.py used MIN);
#   - the fourth-moment zeta large-value bound LV_z <= tau - 4(sigma - 1/2);
#   - Ivic's large-value lemma (11.74) (ANTEDB Lemma a-ivt-1) for m = 2..6 with known exponent pairs;
#   - new mode: the zeta pair (B - 1/2, 1/2) in (11.74) and in zeta-from-exp (optimistic: admissibility for a
#     zeta-only pair is not checked; adding it can only lower the new bound), plus DensityB's inputs
#     (Montgomery range from Lemma 1 via beta_new, LV_zeta = -inf above B).
# Usage: python3 explore2.py [sigma ...]   (standard library only; about 10 s per sigma)
import sys, os, math, importlib.util
sys.dont_write_bytecode = True  # do not leave __pycache__ in ../DensityB/repro
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('ex1', os.path.join(HERE, '..', 'DensityB', 'repro', 'explore.py'))
ex1 = importlib.util.module_from_spec(spec); spec.loader.exec_module(ex1)
B = ex1.B
OLD_PAIRS = [(1/6, 2/3), (2/7, 4/7), (97/251, 132/251), (13/84 + 1e-9, 55/84 + 1e-9)]
NEW_PAIRS = OLD_PAIRS + [(B - 0.5, 0.5)]

def ivic_1174(s, t, m, k, l):  # ANTEDB Lemma a-ivt-1 (Ivic (11.74)), 3/4 < s <= 1, m >= 2
    if s <= 0.75: return math.inf
    X = 2*t/3 + 4*m*(3 - 4*s)/3
    Y = max(t + 3*m*(3 - 4*s), ((k + l)*t + m*(1 + 2*k + 2*l)*(3 - 4*s))/k)
    return max(2 - 2*s, m*(2 - 4*s) + m*t, min(X, Y))

def L_base2(s, t, new, mc):
    v = min(f(s, t) for _, f in ex1.LIT)
    v = min(v, max(2 - 2*s, 2 - 2*s + t - mc))          # HM with beta (old or new)
    for (k, l) in (NEW_PAIRS if new else OLD_PAIRS):
        for m in range(2, 7):
            v = min(v, ivic_1174(s, t, m, k, l))
    return v

def make_Ls(s, new):
    mc = ex1.mc_range_hm(s, ex1.g_new if new else ex1.g_old)
    lvl = 2 - 2*s; t = 0.0; step = 0.0005
    while t < 6 and L_base2(s, t + step, new, mc) <= lvl + 1e-12: t += step
    tM = t
    cache = {}
    def L1(t):
        key = round(t, 6)
        if key not in cache: cache[key] = min(L_base2(s, t, new, mc), lvl + max(0.0, t - tM))
        return cache[key]
    def Ls(t): return min(k*L1(t/k) for k in range(1, 9))
    return Ls, tM, mc

def Z2(s, t, Ls, new):  # zeta large values, t >= 2
    if s > t*ex1.beta_low(1/t): return -math.inf
    if new and s > B: return -math.inf
    v = min(Ls(t), 2*t - 12*s + 6, t - 4*(s - 0.5))      # LV_z <= LV; 12th and 4th moments
    for (k, l) in (NEW_PAIRS if new else OLD_PAIRS):     # zeta-from-exp (Ivic Thm 8.2): MAX of the two terms
        v = min(v, max(t - 6*(s - 0.5), ((k + l)*t - 2*(1 + 2*k + 2*l)*(s - 0.5))/k))
    return v

def A2(s, new):
    Ls, tM, mc = make_Ls(s, new)
    ts = [1.0 + 0.004*i for i in range(0, 1500)]
    Lr = {t: Ls(t)/t for t in ts}
    Zr = {t: Z2(s, t, Ls, new)/t for t in ts if t >= 2}
    best = math.inf; bt = None
    for t0 in [1.0 + 0.004*i for i in range(0, 700)]:
        sup1 = max(Lr[t] for t in ts if t0 - 1e-9 <= t <= 2*t0 + 1e-9)
        sup2 = max([Zr[t] for t in Zr if t < t0 - 1e-9] or [-math.inf])
        v = max(sup1, sup2)
        if v < best: best, bt = v, t0
    return best/(1 - s), bt, tM, mc

def record(s):  # ANTEDB Table "Current best upper bound on A(sigma)" (main, read 2026-11-05)
    if s < 0.76: return 15/(3 + 5*s)
    if s < 127/167: return 9/(8*s - 2)
    if s < 13/17: return 15/(13*s - 3)
    if s < 17/22: return 6/(5*s - 1)
    if s < 41/53: return 2/(9*s - 6)
    if s < 7/9: return 9/(7*s - 1)
    if s < 1867/2347: return 9/(8*(2*s - 1))
    if s < 7/8: return 3/(2*s)
    return float('nan')

if __name__ == '__main__':
    sig = [float(x) for x in sys.argv[1:]] or [0.75, 0.77, 0.79, 0.80, 0.82, 0.84, 0.85, 0.855, 0.86, 0.87]
    print('sigma   record    pipeline-old  pipeline-new  new-old    tauMC(old,new)   tauB')
    for s in sig:
        a0 = A2(s, False); a1 = A2(s, True)
        print(f'{s:.4f}  {record(s):.5f}  {a0[0]:.5f}      {a1[0]:.5f}      {a1[0]-a0[0]:+.5f}   {a0[2]:.4f} {a1[2]:.4f}   {ex1.tauB(s):.4f}', flush=True)
