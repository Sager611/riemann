#!/usr/bin/env python3
"""serrenumK follow-up certificate (python-flint Arb): Lemma 2 for every number field k (Momose Type 2 isogeny primes).
Run: cd Riemann/Deep/SerreNumK && uv run --no-project --with python-flint python3 cert2.py > logs/cert2_output.txt
For k with Galois closure M (degree n_M, discriminant Delta_M) and y = log p:
  R_F(y) = 2 log Delta_M + n_M (y + log 4)              >= log Delta of F = M(sqrt(-p))
  Q1     = 32 log G(R_F)                                >= log q1   (chebexp over Q, class {1} of F)
  R_2    = 2 R_F + 2 n_M (log 4 + Q1)                   >= log Delta of F(sqrt q1)   ([F:Q] <= 2 n_M)
  Q2     = 32 log G(R_2)                                >= log q2   (q2 != q1, chebexp, class {1})
Certified: Q1 + log 4 < y and Q2 + log 4 < y for every y >= y*(k) (intervals, then the closed-form tail of Lemma 2').
G = envelope of chebexp's piecewise constants (as in cert.py). Comparisons are certain Arb comparisons.
"""
import time, math
from fractions import Fraction as Fr
from flint import arb, ctx
ctx.prec = 200
T0 = time.time()
L3, L4 = arb(3).log(), arb(4).log()
C = {k: arb(v) for k, v in [('c0', '0.484'), ('c1', '0.29'), ('c2', '0.151'), ('c3', '0.108'), ('c4', '0.094')]}
OK = []
def check(name, cond):
    OK.append(bool(cond)); print(('PASS ' if cond else 'FAIL ') + name)
def up(x): return x.mid() + x.rad()
def amax(x, y): return (x + y + abs(x - y)) / 2
def amin(x, y): return (x + y - abs(x - y)) / 2
def G_up(R):
    Ru = up(R); m = C['c0'] * (amin(Ru, arb(10)) + L3)
    for start, cap, key in ((10, 100, 'c1'), (100, 1000, 'c2'), (1000, 10**4, 'c3')):
        if not (Ru < start): m = amax(m, C[key] * (amin(Ru, arb(cap)) + L3))
    if not (Ru < 10**4): m = amax(m, C['c4'] * (Ru + L3))
    return m
check('envelope caps (as cert.py): 0.484(10+log3) < 5.374 < 152.3, 0.29(100+log3) < 29.32, 0.151(1000+log3) < 151.17, 0.108(1e4+log3) < 1080.12',
      C['c0'] * (10 + L3) < arb('5.374') and C['c1'] * (100 + L3) < arb('29.32') and C['c2'] * (1000 + L3) < arb('151.17') and C['c3'] * (10**4 + L3) < arb('1080.12'))
check('ghat(R) = 0.108R + 152.3 >= G(R) (R >= 0) and ghat(R) <= R for R >= 171', arb('152.3') > arb('151.17') and C['c3'] * L3 < arb('152.3') and C['c4'] * L3 < arb('152.3') and arb('1080.12') < C['c3'] * 10**4 + arb('152.3') and arb('0.892') * 171 > arb('152.3'))

def interval_ok(ya, yb, nM, lDM):
    RF = 2 * lDM + nM * (yb + L4)
    Q1 = 32 * G_up(RF).log()
    R2 = 2 * RF + 2 * nM * (L4 + Q1)
    Q2 = 32 * G_up(R2).log()
    return bool((Q1 + L4 < ya) and (Q2 + L4 < ya))

def certify(name, nM, DM):
    lDM = arb(DM).log(); n = arb(nM)
    Y0 = max(3000 * nM, math.ceil(float(up(2 * lDM + n * L4))) + 1)
    y0 = arb(Y0)
    # closed-form tail (Lemma 2'): valid for all y >= Y0 by monotonicity (DERIVATION2 section 2)
    t = (y0 > 2 * lDM + n * L4) and (2 * L4 + 64 * (2 * n * y0).log() < y0) and (32 * (5 * n * y0).log() + L4 < y0) and (n * y0 > 171)
    # interval chain from 100 to Y0: width 1/4 below 2000, width 4 above
    grid = [Fr(100) + Fr(k, 4) for k in range(0, 4 * 1900 + 1)] + [Fr(2000 + 4 * k) for k in range(1, (Y0 - 2000) // 4 + 2)]
    grid = [g for g in grid if g <= Y0 + 4]
    res = [interval_ok(arb(a.numerator) / a.denominator, arb(b.numerator) / b.denominator, nM, lDM) for a, b in zip(grid[:-1], grid[1:])]
    lf = max([i for i, r in enumerate(res) if not r], default=-1)
    ystar = grid[lf + 1]
    good = t and lf + 1 < len(res) and all(res[lf + 1:]) and grid[-1] >= Y0
    OK.append(good)
    return float(ystar), Y0, good, len(res) - lf - 1

def lv_route(d, lDk):
    """float comparison only: LV's route with chebexp (over k) on E' = k(sqrt(+-l), sqrt(+-3)): l <= (1 + G(R_E)^16)^4,
    R_E = 4 log Delta_k + 2d log(12 l) (BD's bound for log Delta_E'); returns the least y on a 1/4 grid beyond which it fails"""
    def Gf(R):
        v = [0.484 * (min(R, 10) + math.log(3))]
        for s, cap, c in ((10, 100, 0.29), (100, 1000, 0.151), (1000, 1e4, 0.108)):
            if R >= s: v.append(c * (min(R, cap) + math.log(3)))
        if R >= 1e4: v.append(0.094 * (R + math.log(3)))
        return max(v)
    last = 0
    for k in range(4 * 100, 4 * 5000):
        y = k / 4; RE = 4 * lDk + 2 * d * (math.log(12) + y + 0.25)
        if 4 * math.log1p(Gf(RE) ** 16) >= y: last = y + 0.25
    return last

# BD Table 'tab:output_summary_type_2' (Galois, class number one, no imaginary quadratic HCF): d, Delta_K, label, GRH Type 2 bound
BD = [(2, 5, '2.2.5.1', '5.65e10'), (3, 49, '3.3.49.1', '4.09e11'), (4, 125, '4.0.125.1', '1.46e12'), (5, 14641, '5.5.14641.1', '4.75e12'),
      (6, 300125, '6.6.300125.1', '1.12e13'), (7, 594823321, '7.7.594823321.1', '2.65e13'), (8, 64000000, '8.0.64000000.2', '4.16e13'),
      (9, 16983563041, '9.9.16983563041.1', '7.60e13'), (10, 572981288913, '10.10.572981288913.1', '1.24e14')]
# non-Galois cubics: Galois closure degree 6, |Delta_M| from PARI nfdisc(nfsplitting(.)) (logs/galclosure.txt)
NG = [('Q(2^(1/3))', 3, 108, 6, 34992), ('x^3-x^2-2x-20 (BD Thm 1.8)', 3, 2891, 6, 493114979)]
print(f'\n{"field":>28} {"d":>3} {"n_M":>4} {"log Delta_M":>11} {"y*(k)":>7} {"p >= 10^":>9} {"chain":>6} {"Y0":>6} {"GRH (BD)":>9} {"LV route (float)":>16}')
rows = [(lab, d, d, Dk, Dk, g) for d, Dk, lab, g in BD] + [(nm, d, nM, Dk, DM, '-') for nm, d, Dk, nM, DM in NG]
for lab, d, nM, Dk, DM, g in rows:
    ys, Y0, good, nint = certify(lab, nM, DM)
    lv = lv_route(d, math.log(Dk))
    print(f'{lab:>28} {d:>3} {nM:>4} {math.log(DM):>11.3f} {ys:>7.2f} {ys / math.log(10):>9.3f} {"PASS" if good else "FAIL":>6} {Y0:>6} {g:>9} {"y=" + format(lv, ".2f") + " (10^" + format(lv / math.log(10), ".1f") + ")":>16}')
print('\nLemma 2 (closed form, every k): no Momose Type 2 isogeny prime with log p >= max(2 log Delta_M + n_M log 4, 3000 n_M) (proof: DERIVATION2 section 2)')
print(f'\n{"ALL CERTIFIED" if all(OK) else "SOME CHECK FAILED"} ({time.time() - T0:.1f} s)')
