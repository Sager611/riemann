#!/usr/bin/env python3
"""densityb verify.py: exact rational checks for Riemann/Deep/DensityB/PROOF.md.
Run: uv run --no-project python3 verify.py   (standard library only)."""
from fractions import Fraction as F
import time
t_start = time.time()
FAILS = []
def check(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond: FAILS.append(name)
def dec(x, d=12):
    # decimal string of a Fraction, truncated toward zero, d digits
    s = '-' if x < 0 else ''; x = abs(x); ip = x.numerator // x.denominator
    fr = x - ip; digs = ''
    for _ in range(d):
        fr *= 10; dgt = fr.numerator // fr.denominator; digs += str(dgt); fr -= dgt
    return f"{s}{ip}.{digs}"

B = F(87338196661320752, 10**17)          # Lean-checked half-plane Re s > B (Bebel, 31 Oct 2026)
half = F(1, 2)
mu_half = F(13, 84)                        # Bourgain (2017): mu(1/2) <= 13/84

# --- consistency of B with Liu's contact B2 = (777 - sqrt(1081))/852 (B must lie above it)
# B > B2  <=>  777 - 852 B < sqrt(1081)  <=>  (777-852B)^2 < 1081 (777-852B > 0)
lhs = 777 - 852 * B
check("B above Liu's contact B2 = (777-sqrt1081)/852 (0 < 777-852B and (777-852B)^2 < 1081)", lhs > 0 and lhs * lhs < 1081)
check("B < 7/8", B < F(7, 8))

# --- (C1') reflected exponent: mu(1-B) = B - 1/2 (functional equation); zeta-sum exponent line
def gB(tau):            # sum_{n in I} n^{-it} << N^{gB(tau)+eps}, |t| = N^tau, 1 <= tau <= 2
    return 1 - B + (B - half) * tau
check("gB(1) = 1/2", gB(F(1)) == half)
check("gB(2) = B", gB(F(2)) == B)

def tauB(s):  return (2 * s - 2 + B) / (B - half)
def tau0(s):  return tauB(s) + 1 - s
def Anew(s):  return 3 / tau0(s)
def Aclosed(s): return 3 * (B - half) / ((F(5, 2) - B) * s - (F(5, 2) - 2 * B))
sstar = (5 - 4 * B) / (7 - 6 * B)
print("B      =", dec(B, 17))
print("sigma* =", dec(sstar, 15), " (5-4B)/(7-6B)")
check("3/4 < sigma* < B", F(3, 4) < sstar < B)
check("B < (1+B)/2 (HM range for Theorem 1 contains [sigma*, B])", B < (1 + B) / 2)

# grid of exact rationals in [sigma*, B] including endpoints
grid = [sstar + (B - sstar) * F(i, 400) for i in range(401)]
ok = {k: True for k in "abcdefghij"}
for s in grid:
    tb, t0 = tauB(s), tau0(s)
    ok['a'] &= (F(1) <= tb <= 2)                                   # tauB in [1,2]
    ok['b'] &= (gB(tb) == 2 * s - 1)                               # HM threshold exactly at tauB
    ok['c'] &= (t0 > 3 - 3 * s)                                    # Cor zero-large-cor3 needs tau0 >= 3-3s
    # (lvo) on [2 tau0/3, tau0]: max(2-2s, 2-2s+tau-tauB) <= (3-3s) tau/tau0; convex - linear => check endpoints
    for tau in (2 * t0 / 3, t0):
        lv = max(2 - 2 * s, 2 - 2 * s + tau - tb)
        ok['d'] &= (lv <= (3 - 3 * s) * tau / t0)
    # (lvoz) for 2 <= tau < 4 tau0/3: LV_zeta = -inf when s > 1/2 + tau*13/84 (lvz-mu with Bourgain)
    ok['e'] &= (4 * t0 / 3 <= (s - half) / mu_half)
    # alternative (lvoz) via Heath-Brown twelfth moment: 2tau-12(s-1/2) <= (3-3s)tau/tau0 on [2,4tau0/3] iff tau0 <= 3s-3/4
    ok['f'] &= (t0 <= 3 * s - F(3, 4))
    ok['g'] &= (Anew(s) == Aclosed(s))
    ok['h'] &= (Anew(s) <= F(3, 2) / s) and ((s == sstar) == (Anew(s) == F(3, 2) / s))
    # beats the routine HM bound with the B-image (178/641, 365/641) of TTY's pair: MC for tau < (641 s - 414)/89
    ok['i'] &= (tb > (641 * s - 414) / 89)
    ok['j'] &= (4 * t0 / 3 > 2)                                    # so (lvoz) is not vacuous; it is checked in e/f
check("(a) 1 <= tauB(s) <= 2 on [sigma*, B]", ok['a'])
check("(b) gB(tauB(s)) = 2s-1 exactly (HM threshold)", ok['b'])
check("(c) tau0 > 3-3s", ok['c'])
check("(d) (lvo): LV <= (3-3s) tau/tau0 on [2tau0/3, tau0]", ok['d'])
check("(e) (lvoz) via mu(1/2)<=13/84: 4tau0/3 <= (s-1/2)*84/13", ok['e'])
check("(f) (lvoz) alternatively via HB twelfth moment: tau0 <= 3s-3/4", ok['f'])
check("(g) 3/tau0 equals closed form 3(B-1/2)/((5/2-B)s-(5/2-2B))", ok['g'])
check("(h) Anew <= 3/(2s), equality only at sigma*", ok['h'])
check("(i) tauB(s) > (641s-414)/89 (beats routine HM with (178/641,365/641))", ok['i'])
check("(j) 4tau0/3 > 2", ok['j'])
check("tauB(sigma*) = 3 sigma* - 1 (crossover with Jutila k=2)", tauB(sstar) == 3 * sstar - 1)
check("tauB(s) > 3s-1 for s in (sigma*, B]", all(tauB(s) > 3 * s - 1 for s in grid[1:]))

# --- old best zeta-sum exponents from ANTEDB best beta bounds (compute_best_beta_bounds, 2026-11-05 main):
# beta(a) <= c + m a on [x0,x1] (a <= 1/2); for tau in [1,2]: g_old(tau) = 1 - tau/2 + tau*beta(1 - 1/tau)
pieces = [(F(423,1295), F(227,601), F(89,1282), F(454,641)),
          (F(227,601), F(12,31), F(29,280), F(173,280)),
          (F(12,31), F(1508,3825), F(1,32), F(103,128)),
          (F(1508,3825), F(62831,155153), F(18,199), F(521,796))]
def g_line(c, m):  # g = (1-m) + tau (c + m - 1/2)
    return (1 - m, c + m - half)
lo_tau, hi_tau = tauB(sstar), tauB(B)
cov_lo = 1 / (1 - pieces[0][0]); cov_hi = 1 / (1 - pieces[-1][1])
check("beta pieces cover tau in [tauB(sigma*), tauB(B)]", cov_lo <= lo_tau and hi_tau <= cov_hi)
okg = True
for x0, x1, c, m in pieces:
    t_a, t_b = 1 / (1 - x0), 1 / (1 - x1)
    a0, a1 = g_line(c, m)
    for t in (t_a, t_b):
        okg &= (a0 + a1 * t > gB(t))
check("g_old(tau) > gB(tau) on [1.485, 1.68] (new zeta-sum bound strictly better there)", okg)

# --- (C2) and (C3) sample values
def mu_line(a): return mu_half * (B - a) / (B - half)
print("C2 mu(0.8) <=", dec(mu_line(F(4,5)), 6), " mu(0.85) <=", dec(mu_line(F(17,20)), 6))
m34 = F(12408) / (4537 - 4890 * F(3, 4))                      # Ivic Thm 8.4 at sigma = 3/4
m78_new = m34 * (B - F(3, 4)) / (B - F(78, 100))
m78_ivic = F(12408) / (4537 - 4890 * F(78, 100))
print("C3 m(3/4) =", dec(m34, 4), " m(0.78) >=", dec(m78_new, 4), " Ivic m(0.78) =", dec(m78_ivic, 4))
check("C3 example m(0.78) >= 18.8 > Ivic 17.1", m78_new > F(188, 10) and m78_ivic < F(172, 10))
# no-gain of (C3) for LV_zeta below B: m(s)(sig - s) <= m1 (sig - s1) for s1 < s < B, sig < B
okng = True
for s1 in (half, F(3,5), F(3,4)):
    for s in (s1 + (B - s1) * F(k, 10) for k in range(1, 10)):
        for sig in (s1 + (B - s1) * F(k, 20) for k in range(1, 20)):
            if sig <= s: continue
            okng &= ((B - s1) / (B - s) * (sig - s) <= sig - s1)
check("C3 no-gain: (B-s1)/(B-s)*(sig-s) <= sig-s1 for s1<s<sig<B", okng)
# (C2) no gain through lvz-mu: on the segment mu=c(B-a), (sig-a)/(c(B-a)) is increasing as a decreases (sig<B)
okc2 = True
c = mu_half / (B - half)
for sig in (F(80,100), F(85,100), F(87,100)):
    vals = [(sig - a) / (c * (B - a)) for a in (half + (sig - half) * F(k, 50) for k in range(0, 50))]
    okc2 &= all(vals[i] > vals[i + 1] for i in range(len(vals) - 1))
check("C2 no-gain: (sig-a)/mu_line(a) maximal at the left end a=1/2", okc2)

# --- table of values
print()
print("sigma        old 3/(2s)     routine-HM 267/(552s-325)   new 3/tau0       tauB")
for s in [sstar, F(86,100), F(865,1000), F(325,374), F(87,100), F(8725,10000), B]:
    old = F(3, 2) / s
    hm = (F(267) / (552 * s - 325)) if s > F(325, 374) else None
    print(f"{dec(s,6):>10}   {dec(old,5):>10}   {dec(hm,5) if hm else '   (n/a)':>12}            {dec(Anew(s),5):>10}   {dec(tauB(s),5)}")
AB = Anew(B)
print("A(B-) <=", dec(AB, 6), " = 3(B-1/2)/(9B/2 - B^2 - 5/2);  3/(2B) =", dec(F(3,2)/B, 6), "; A(B-)(1-B) =", dec(AB * (1 - B), 6))
check("A(B-) closed form", AB == 3 * (B - half) / (F(9, 2) * B - B * B - F(5, 2)))
print()
print("runtime %.3f s" % (time.time() - t_start))
print("ALL PASS" if not FAILS else "FAILURES: " + ", ".join(FAILS))
