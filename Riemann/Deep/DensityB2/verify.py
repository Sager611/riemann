#!/usr/bin/env python3
# densityb2 verify.py: exact rational checks behind REPORT.md / PROOF.md (standard library only).
# Claim certified: on [3/4, sigma*) the half-plane inputs (H) do not lower any entry of the A(sigma) table through
# the channels of PROOF.md; the binding inputs are (H)-independent. Run: python3 verify.py
from fractions import Fraction as F
import time
t_start = time.time()
B = F(87338196661320752, 10**17); h = F(1, 2)
ok = True
def check(name, cond):
    global ok
    print(('PASS ' if cond else 'FAIL ') + name); ok = ok and bool(cond)
def tauB(s): return (2*s - 2 + B)/(B - h)
sstar = (5 - 4*B)/(7 - 6*B)
check('sigma* = (5-4B)/(7-6B) solves tauB = 3 sigma - 1', tauB(sstar) == 3*sstar - 1)
print('   sigma* =', float(sstar))

# ---------- Part A: the table and its binding inputs ----------
# Ivic's Lemma gzd (ANTEDB): A <= 3/tau0_m, tau0_m = ((3m-2)s + 2 - m)/m; MC cap tau_m = tau0_m + s - 1.
def tau0m(m, s): return ((3*m - 2)*s + 2 - m)/F(m)
def taum(m, s): return tau0m(m, s) + s - 1
table_A = {2: lambda s: 3/(2*s), 3: lambda s: 9/(7*s - 1), 4: lambda s: 6/(5*s - 1), 5: lambda s: 15/(13*s - 3), 6: lambda s: 9/(8*s - 2)}
for m, f in table_A.items():
    check(f'A1 m={m}: 3/tau0_m equals the table formula (two points)', all(3/tau0m(m, s) == f(s) for s in (F(4, 5), F(17, 20))))
    # A2: tau_m is where the diagonal term m(tau + 2 - 4s) of the m-th power estimate reaches 2 - 2s
    check(f'A2 m={m}: m(tau_m + 2 - 4s) = 2 - 2s', all(m*(taum(m, s) + 2 - 4*s) == 2 - 2*s for s in (F(4, 5), F(17, 20))))
def c_m(m): return F(6*m*m - 5*m + 2, 8*m*m - 7*m + 2)
def a_m(m): return F(9*m*m - 4*m + 2, 12*m*m - 6*m + 2)
def g_m(m, k, l):
    return (3*m*m*(1 + 2*k + 2*l) - (4*k + 2*l)*m + 2*k + 2*l)/(4*m*m*(1 + 2*k + 2*l) - (6*k + 4*l)*m + 2*k + 2*l)
known_c = {2: F(4, 5), 3: F(41, 53), 4: F(13, 17), 5: F(127, 167), 6: F(47, 62)}
for m in range(2, 7):
    # A3: c_m is where Heath-Brown's double-zeta off-diagonal term X = 2tau/3 + 4m(3-4s)/3 equals 2-2s at tau_m
    s = c_m(m); X = 2*taum(m, s)/3 + 4*m*(3 - 4*s)/3
    check(f'A3 m={m}: c_m = {c_m(m)} (table/Ivic start) and X(tau_m) = 2 - 2s there', c_m(m) == known_c[m] and X == 2 - 2*s)
check('A3 sanity: Ivic pair (97/251,132/251) gives 3831/4791', g_m(2, F(97, 251), F(132, 251)) == F(3831, 4791))
check('A3 sanity: Lindelof pair (0,1/2) gives 23/29 (Kerr\'s unpublished start)', g_m(2, F(0), h) == F(23, 29))
# A4: the zeta pair (B-1/2, 1/2) in Lemma gzd
kz, lz = B - h, h
sI = g_m(2, kz, lz)
print('   m=2 with zeta pair: 3/(2 sigma) from sigma >=', float(sI))
check('A4 m=2: zeta pair start > Bourgain 1867/2347 (no record) and < Ivic 3831/4791', F(1867, 2347) < sI < F(3831, 4791))
for m in range(3, 7):
    check(f'A4 m={m}: zeta pair leaves the range start c_m unchanged (c_m <= max(a_m, g_m))', c_m(m) <= max(a_m(m), g_m(m, kz, lz)))
# A5: DensityB's channel tau_B gives A <= 3/(tau_B + 1 - s); compare with each table entry's tau0 on its range.
# All tau0 here are linear in s, so endpoint checks suffice.
rows = [  # (lo, hi, tau0_record(s), name)
    (F(3, 4), F(19, 25), lambda s: (3 + 5*s)/5, '15/(3+5s) Guth-Maynard'),
    (F(19, 25), F(127, 167), lambda s: tau0m(6, s), '9/(8s-2)'),
    (F(127, 167), F(13, 17), lambda s: tau0m(5, s), '15/(13s-3)'),
    (F(13, 17), F(17, 22), lambda s: tau0m(4, s), '6/(5s-1)'),
    (F(17, 22), F(41, 53), lambda s: 3*(9*s - 6)/2, '2/(9s-6) TTY-Bourgain'),
    (F(41, 53), F(7, 9), lambda s: tau0m(3, s), '9/(7s-1)'),
    (F(7, 9), F(1867, 2347), lambda s: 8*(2*s - 1)/3, '9/(8(2s-1)) TTY-Bourgain'),
    (F(1867, 2347), sstar, lambda s: 2*s, '3/(2s) Bourgain/Ivic'),
]
for lo, hi, t0, name in rows:
    okr = all(tauB(s) + 1 - s <= t0(s) for s in (lo, hi)) and tauB(lo) + 1 - lo < t0(lo)
    check(f'A5 {name} on [{float(lo):.4f},{float(hi):.4f}]: tau_B + 1 - s <= tau0_record (half-plane HM channel gives nothing)', okr)

# ---------- Part B: every sup-norm channel of (H) reproduces tau_B ----------
# B1 mu-lv: MC for tau <= (2s-1-th)/mu(th). New hull point (1-B, B-1/2); between it and (1/2, 13/84) mu is <= linear;
# for th in [0, 1-B] mu(th) = 1/2 - th. Ratio of linear functions -> endpoint; check the endpoint values.
for s in (F(3, 4) + F(1, 1000), F(4, 5), F(17, 20), sstar):
    at_1mB = (2*s - 1 - (1 - B))/(B - h); at_half = (2*s - 1 - h)/F(13, 84)
    check(f'B1 mu-lv at s={float(s):.4f}: best endpoint th=1-B gives tau_B ({float(at_1mB):.4f} >= {float(at_half):.4f})', at_1mB == tauB(s) and at_1mB >= at_half)
# B2 reflection principle, k=1, (H) off-diagonal 2rho + 2B alpha (alpha = tau - 1 <= tau/2): contradiction iff tau < tau_B
for s in (F(4, 5), F(17, 20)):
    check(f'B2 k=1 reflect + (0,B): threshold 1 + (4s-3)/(2B-1) = tau_B (s={float(s)})', 1 + (4*s - 3)/(2*B - 1) == tauB(s))
# B3 k=1 reflect + Heath-Brown double zeta (5rho/4 + tau/2 term): MC iff tau <= 5s - 3; diagonal: tau <= 2s
for s in (F(4, 5), F(17, 20)):
    rho = 2 - 2*s; tau = 5*s - 3
    check(f'B3 k=1 reflect + HB: rho = 2-2s at tau = 5s-3 (s={float(s)})', 4 - F(16, 3)*s + F(2, 3)*tau == rho)
# B4 Heath-Brown's double zeta bound is diagonal-optimal (s <= rho + 2) for tau <= 3/2, 0 <= rho <= 1 (vertex check)
check('B4 max(rho+1, 2rho, 5rho/4+tau/2) = rho+1 on tau<=3/2, 0<=rho<=1',
      all(max(r + 1, 2*r, F(5, 4)*r + t/2) == r + 1 for r in (F(0), F(1)) for t in (F(0), F(3, 2))))
check('B4 at the binding points (tau = 3s-1 <= 3sigma*-1, rho = 2-2s): 5rho/4 + tau/2 <= rho + 1',
      all(F(5, 4)*(2 - 2*s) + (3*s - 1)/2 <= (2 - 2*s) + 1 for s in (F(3, 4), sstar)))
# B5 combined thresholds below sigma*: max(tau_B, 5s-3) < 3s-1 (k=2 diagonal cap) on [4/5, sigma*)
for s in (F(4, 5), F(17, 20), sstar - F(1, 10**6)):
    check(f'B5 s={float(s):.6f}: max(tau_B, 5s-3) < 3s-1', max(tauB(s), 5*s - 3) < 3*s - 1)
# B6 the (H) off-diagonal condition at tau = 3s-1 fails exactly below sigma*: 3 - 4s + (2B-1)(3s-2) > 0
for s in (F(4, 5), F(17, 20)):
    check(f'B6 s={float(s)}: (H) pointwise off-diagonal too weak at tau = 3s-1', 3 - 4*s + (2*B - 1)*(3*s - 2) > 0)
check('B6 equality exactly at sigma*', 3 - 4*sstar + (2*B - 1)*(3*sstar - 2) == 0)
# B7 reflection principle (Ivic 11.58-11.59): corrected form rho <= k(3-4s)/2 + (s' - alpha)/2 reproduces the
# known terms; the ANTEDB blueprint's printed form rho <= k(3-4s)/2 + s' - 1 does not (typo there).
for k in (1, 2, 3):
    s = F(17, 20); tau = F(3, 2); al = k*(tau - 1)
    # diagonal s' = rho + 2 alpha: rho <= k(3-4s) + alpha = k(tau + 2 - 4s)
    rho = k*(3 - 4*s) + al
    check(f'B7 k={k}: corrected reflect, diagonal case gives k(tau+2-4s)', rho == k*(tau + 2 - 4*s) and rho == F(k, 2)*(3 - 4*s) + ((rho + 2*al) - al)/2)
    # HB off-diagonal s' = 5rho/4 + tau/2 + alpha gives the X-term 2tau/3 + 4k(3-4s)/3
    X = 2*tau/3 + F(4*k, 3)*(3 - 4*s)
    check(f'B7 k={k}: corrected reflect, HB term gives X', X == F(k, 2)*(3 - 4*s) + ((F(5, 4)*X + tau/2 + al) - al)/2)
    # printed form with the diagonal: rho <= k(3-4s)/2 + (rho + 2alpha) - 1, i.e. RHS - rho does not depend on rho
    gap = lambda r: F(k, 2)*(3 - 4*s) + (r + 2*al) - 1 - r
    check(f'B7 k={k}: printed form is vacuous in the diagonal case (RHS - rho independent of rho)', gap(F(0)) == gap(F(1)))

# ---------- Part C: the zeta large-value side has slack (Type II not binding) ----------
# For 3/(2s) (tau0 = 2s): need LV_z(s,tau) <= (3-3s)tau/(2s) for 2 <= tau < 8s/3; twelfth moment 2tau - 12(s - 1/2)
# has larger slope, so check tau = 8s/3: margin (4-4s) - (6 - 20s/3) = 8s/3 - 2 > 0.
for s in (F(1867, 2347), F(4, 5), sstar):
    marg = (4 - 4*s) - (2*F(8, 3)*s - 12*(s - h))
    check(f'C1 s={float(s):.4f}: twelfth-moment margin at tau=8s/3 is {float(marg):.4f} > 0', marg > 0 and marg == F(8, 3)*s - 2)
# For the m-entries (tau0_m, 4tau0_m/3 vs 2) and the TTY-Bourgain entries, the same twelfth-moment check:
for lo, hi, t0, name in rows[:-1]:
    for s in (lo, hi):
        T0 = t0(s); top = 4*T0/3
        if top <= 2:
            good = True  # (lvoz) vacuous
        else:
            good = (2*top - 12*(s - h)) <= (3 - 3*s)*top/T0 and (4 - 12*(s - h)) <= (3 - 3*s)*2/T0
        check(f'C2 {name} at s={float(s):.4f}: twelfth moment covers (lvoz)', good)

# ---------- values ----------
print()
print('sigma     record     DensityB-channel 3/(tauB+1-s)   tauB     3s-1    5s-3')
for lo, hi, t0, name in rows:
    for s in (lo, (lo + hi)/2):
        print(f'{float(s):.4f}   {float(3/t0(s)):.5f}    {float(3/(tauB(s) + 1 - s)):.5f}                    {float(tauB(s)):.4f}   {float(3*s - 1):.4f}  {float(5*s - 3):.4f}')
print(f'\nruntime {time.time() - t_start:.3f} s')
print('ALL PASS' if ok else 'SOME FAIL')
