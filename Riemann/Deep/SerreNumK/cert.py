#!/usr/bin/env python3
"""serrenumK certificate (python-flint Arb). Run:
  cd Riemann/Deep/SerreNumK && uv run --no-project --with python-flint python3 cert.py > logs/cert_output.txt
Checks (DERIVATION.md section numbers):
  S2  exact coefficients of R(K, rho) = paper's bound for log d_{L~} (Fractions + integer identity for the logs)
  S3  envelope G(R) = sup_{Lambda<=R} c(Lambda)(Lambda+log 3) <= 0.108(R+log 3) (R>=1400), <= 0.094(R+log 3) (R>=11490)
  S3  Theorem A constants a, b(K) (rounded up) for K = Q and K = Q(sqrt D), D squarefree, 2 <= D <= 47; table
  S5  Momose Type 2 threshold y* (Condition CC + chebexp twice): no Type 2 isogeny prime p with log p >= y*
All decimal constants enter as Arb balls from strings; every comparison is a strict Arb comparison (True only if certain).
"""
import time, math
from fractions import Fraction as Fr
from flint import arb, ctx
ctx.prec = 200
T0 = time.time()
L3 = arb(3).log()
C = {k: arb(v) for k, v in [('c0', '0.484'), ('c1', '0.29'), ('c2', '0.151'), ('c3', '0.108'), ('c4', '0.094')]}
OK = []
def check(name, cond):
    OK.append(bool(cond)); print(('PASS ' if cond else 'FAIL ') + name)
def up(x): return x.mid() + x.rad()
def lo(x): return x.mid() - x.rad()

# ---------- S2: R(K, rho) ----------
# basis: rho, n*rho, ld, n*ld, n*log n, n*log192, n*log4   (n = n_K, ld = log d_K)
def add(*ds):
    out = {}
    for d in ds:
        for k, v in d.items(): out[k] = out.get(k, 0) + v
    return out
def scale(c, d): return {k: c * v for k, v in d.items()}
logdL = {'ld': 192, 'rho': 191, 'nlog192': 192}                         # (4.4), [L:K] <= 192
logNalpha = {'nlogn': 2, 'nld': Fr(1, 2), 'nrho': 2}                      # (4.2), Lemma 'norm bound'
R_terms = add(scale(2, logdL), scale(192, add({'nlog4': 1}, logNalpha)))  # (4.6) with [L:K] <= 192
print('R terms:', {k: str(v) for k, v in R_terms.items()})
check('R: rho coeff 382, n*rho coeff 384', R_terms['rho'] == 382 and R_terms['nrho'] == 384)
check('R: log d_K coeff 384 + 96 n_K', R_terms['ld'] == 384 and R_terms['nld'] == 96)
check('R: 384 n log192 + 192 n log4 = 384 n log384 (384^2 = 4*192^2) and 384 n log n', R_terms['nlog192'] == 384 and R_terms['nlog4'] == 192 and 384**2 == 4 * 192**2 and R_terms['nlogn'] == 384)
def R_of(n, dK, rho):  # arb
    n = arb(n); return (382 + 384 * n) * rho + (384 + 96 * n) * arb(dK).log() + 384 * n * (384 * n).log()

# ---------- S3: envelope ----------
p0 = C['c0'] * (10 + L3); p1 = C['c1'] * (100 + L3); p2 = C['c2'] * (1000 + L3); p3 = C['c3'] * (10**4 + L3)
print('pieces (upper):', [round(float(up(x)), 4) for x in (p0, p1, p2, p3)])
check('piece bounds 5.374, 29.32, 151.17, 1080.12', p0 < arb('5.374') and p1 < arb('29.32') and p2 < arb('151.17') and p3 < arb('1080.12'))
check('0.108(1400 + log 3) > 151.17 (so G(R) <= 0.108(R+log3) for R >= 1400)', C['c3'] * (1400 + L3) > arb('151.17'))
check('0.094(11490 + log 3) > 1080.12 (so G(R) <= 0.094(R+log3) for R >= 11490)', C['c4'] * (11490 + L3) > arb('1080.12'))
check('ghat(R) = 0.108R + 152.3 >= G(R) for R >= 0: 152.3 > 151.17 (capped pieces), 0.108 log3 < 152.3, 0.094 log3 < 152.3, 1080.12 < 0.108e4 + 152.3', arb('152.3') > arb('151.17') and C['c3'] * L3 < arb('152.3') and C['c4'] * L3 < arb('152.3') and arb('1080.12') < C['c3'] * 10**4 + arb('152.3'))
def amax(x, y): return (x + y + abs(x - y)) / 2          # ball enclosure of max
def amin(x, y): return (x + y - abs(x - y)) / 2          # ball enclosure of min
def G_up(R):
    """enclosure of an upper bound for G(R') for every R' <= R (G nondecreasing); a piece is kept unless R < its start
    is certain (keeping an inactive piece only enlarges the bound)"""
    Ru = up(R); m = C['c0'] * (amin(Ru, arb(10)) + L3)
    for start, cap, key in ((10, 100, 'c1'), (100, 1000, 'c2'), (1000, 10**4, 'c3')):
        if not (Ru < start): m = amax(m, C[key] * (amin(Ru, arb(cap)) + L3))
    if not (Ru < 10**4): m = amax(m, C['c4'] * (Ru + L3))
    return m

# ---------- S3: Theorem A constants ----------
def ceilto(x, d=3):
    s = 10**d; return math.ceil(float(up(x)) * s - 1e-12 + 1e-9) / s
def sqfree(n): return all(n % (q*q) for q in range(2, int(n**0.5) + 1))
fields = [('Q', 1, 1)] + [(f'Q(sqrt {D})', 2, D if D % 4 == 1 else 4 * D) for D in range(2, 50) if sqfree(D)]
print('\nTheorem A: N p <= (a rho + b)^32 [0.108], and (a1 rho + b1)^32 [0.094] once rho >= rho1;  rho = log N rad(2 O_K N1 N2)')
print(f'{"K":>12} {"d_K":>4} {"a":>9} {"b":>9} {"a1":>9} {"b1":>9} {"rho1":>7}')
CONST = {}
LOG2 = arb(2).log()
for name, n, dK in fields:
    na = arb(n); ld = arb(dK).log()
    a_ex = C['c3'] * (382 + 384 * na); b_ex = C['c3'] * ((384 + 96 * na) * ld + 384 * na * (384 * na).log() + L3)
    a1_ex = C['c4'] * (382 + 384 * na); b1_ex = C['c4'] * ((384 + 96 * na) * ld + 384 * na * (384 * na).log() + L3)
    a, b, a1, b1 = (ceilto(x) for x in (a_ex, b_ex, a1_ex, b1_ex))
    A, B, A1, B1 = (arb(f'{v:.3f}') for v in (a, b, a1, b1))
    rho1 = (arb(11490) - (384 + 96 * na) * ld - 384 * na * (384 * na).log()) / (382 + 384 * na)
    r1 = max(0.0, math.ceil(float(up(rho1)) * 1000) / 1000)
    Rmin = R_of(n, dK, LOG2)
    good = (Rmin > 1400) and (A > a_ex) and (B > b_ex) and (A1 > a1_ex) and (B1 > b1_ex)
    # identity: 0.108(R + log 3) = a_ex rho + b_ex, checked at rho = 1 and 2 (affine)
    for rho in (arb(1), arb(2)):
        good = good and abs(C['c3'] * (R_of(n, dK, rho) + L3) - (a_ex * rho + b_ex)) < arb('1e-40')
    OK.append(good)
    CONST[name] = (n, dK, a, b, a1, b1, r1)
    print(f'{name:>12} {dK:>4} {a:>9.3f} {b:>9.3f} {a1:>9.3f} {b1:>9.3f} {r1:>7.3f} {"PASS" if good else "FAIL"}')

check('Theorem A constants, all fields', all(OK[-len(fields):]))

# ---------- tables (Theorems A, B, Corollary) ----------
def C12(n, dK):  # GRH constants of 2611.19556, C2 as proved (l.1094-1100)
    n = arb(n); ld = arb(dK).log()
    return arb('510.352') + arb('513.024') * n, (arb('513.024') + arb('128.256') * n) * ld + arb('513.024') * n * (384 * n).log() + arb('1324.8') * n + arb('45.7')
LN10 = arb(10).log()
def l10(x): return float(up(x.log() / LN10))

# ---------- S5: Momose Type 2 threshold ----------
def type2_interval(ya, yb, dK):
    RF = 2 * arb(yb) + 2 * arb(dK).log() + arb(16).log()          # log Delta_F <= 2 log p + 2 log d_K + log 16
    Q1 = 32 * G_up(RF).log()                                       # log q1
    R2 = 2 * RF + 4 * arb(4).log() + 4 * Q1                        # log Delta_{F(sqrt q1)} <= 2 log Delta_F + 4 log(4 q1)
    Q2 = 32 * G_up(R2).log()                                       # log q2
    return bool((Q1 + arb(4).log() < ya) and (Q2 + arb(4).log() < ya)), float(up(Q2))
dKmax = max(f[2] for f in fields)
step = Fr(1, 4); Y0 = 2000
grid = [Fr(150) + k * step for k in range(int((Y0 - 150) / step) + 1)]
res = [type2_interval(arb(g.numerator) / g.denominator, arb((g + step).numerator) / (g + step).denominator, dKmax)[0] for g in grid[:-1]]
last_fail = max([i for i, r in enumerate(res) if not r], default=-1)
ystar = grid[last_fail + 1]
print(f'\nType 2 chain: d_K <= {dKmax}, intervals of width {step} on [150, {Y0}]: last failing interval ends at y* = {float(ystar)}')
check(f'Type 2 chain passes on [{float(ystar)}, {Y0}] ({len(res) - last_fail - 1} intervals)', last_fail + 1 < len(res) and all(res[last_fail + 1:]))
y0 = arb(Y0)
check('tail t1: 2 log d_K + log 16 <= Y0 (so R_F <= 3y)', 2 * arb(dKmax).log() + arb(16).log() < y0)
check('tail t2: 0.108R + 152.3 <= R for R >= 171', arb('0.892') * 171 > arb('152.3'))
check('tail t3: 6y + 4 log 4 + 128 log(3y) <= 7y at y = Y0 (difference increasing for y > 128)', 6 * y0 + 4 * arb(4).log() + 128 * (3 * y0).log() < 7 * y0)
check('tail t4: 32 log(7y) + log 4 < y at y = Y0 (difference increasing for y > 32)', 32 * (7 * y0).log() + arb(4).log() < y0)
print(f'=> no Momose Type 2 isogeny prime p over any real quadratic K with d_K <= {dKmax} has log p >= {float(ystar)} (p >= 10^{float(ystar)/math.log(10):.3f})')

print('\nTable: log10 of the bounds (rho = log N rad(2 O_K N...)); A = N p (Thm A); B = 2 sqrt(A) (Thm B, l <= B); GRH = paper (C1 rho + C2)^2 and 2(C1 rho + C2)')
print(f'{"K":>12} {"N rad":>7} {"A":>9} {"A GRH":>7} {"B":>8} {"B GRH":>6} {"Cor":>8}')
for name in ('Q', 'Q(sqrt 2)', 'Q(sqrt 5)', 'Q(sqrt 47)'):
    n, dK = CONST[name][0], CONST[name][1]
    for e in (1, 6, 100):
        rho = e * LN10
        Gv = G_up(R_of(n, dK, rho)); A = 32 * l10(Gv); B = float(up((2 * Gv**16).log() / LN10))
        c1, c2 = C12(n, dK); g = c1 * rho + c2
        corv = max(B, float(ystar) / math.log(10)) if name != 'Q' else B
        print(f'{name:>12} {"10^"+str(e):>7} {A:>9.3f} {2*l10(g):>7.3f} {B:>8.3f} {l10(2*g):>6.3f} {corv:>8.3f}')
# informational (float): N rad above which Theorem B's bound exceeds the Type 2 threshold e^{y*}
for name in ('Q(sqrt 5)', 'Q(sqrt 47)'):
    n, dK = CONST[name][0], CONST[name][1]
    Rc = math.exp((float(ystar) - math.log(2)) / 16) / 0.094 - math.log(3)
    rc = (Rc - (384 + 96 * n) * math.log(dK) - 384 * n * math.log(384 * n)) / (382 + 384 * n)
    print(f'info {name}: 2G(R)^16 > e^y* once rho > {rc:.1f} (N rad > 10^{rc / math.log(10):.1f}) [float]')
print(f'\n{"ALL CERTIFIED" if all(OK) else "SOME CHECK FAILED"} ({time.time() - T0:.1f} s)')
