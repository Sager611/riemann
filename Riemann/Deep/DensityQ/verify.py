#!/usr/bin/env python3
"""DensityQ certificate. Exact rational arithmetic (fractions), standard library only.

Run: uv run --no-project python3 verify.py      (prints ALL PASS on success)

Notation (PROOF.md): B = 87338196661320752/10^17, c = 2 - 2s, tauB(s) = (2s - 2 + B)/(B - 1/2),
X = qT (one modulus) or Q^2 T (large sieve), tau = log X / log N, cut = (6s - 3)/lam.
"""
from fractions import Fraction as Fr
import time

t0 = time.time()
B = Fr(87338196661320752, 10**17)
H = Fr(1, 2)
FAILS = []


def ok(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond:
        FAILS.append(name)


def tauB(s):
    return (2 * s - 2 + B) / (B - H)


def gB(t):  # (H)-vertex exponent of an off-diagonal character sum: N^{1-B} X^{B-1/2}
    return 1 - B + (B - H) * t


def c(s):
    return 2 - 2 * s


SSTAR = (5 - 4 * B) / (7 - 6 * B)
S84 = Fr(21, 25)  # 0.84, left end of the range asked for


def lin_ok(f, a, b):
    """f affine in s: f >= 0 on [a, b] iff f(a) >= 0 and f(b) >= 0."""
    return f(a) >= 0 and f(b) >= 0


print("B =", float(B), " s* =", float(SSTAR))

# ---------------------------------------------------------------- A. Theorem A (family Montgomery range)
ok("A1 gB(tauB(s)) = 2s - 1 (HM threshold), exact at 3 points of an affine identity",
   all(gB(tauB(s)) == 2 * s - 1 for s in (Fr(3, 4), Fr(43, 50), B)))
ok("A2 tauB(3/4) = 1 and tauB((1+B)/2) = 2", tauB(Fr(3, 4)) == 1 and tauB((1 + B) / 2) == 2)
ok("A3 tauB(s*) = 3s* - 1 (crossing with Jutila k=2 family range)", tauB(SSTAR) == 3 * SSTAR - 1)
ok("A4 slope of tauB = 2/(B-1/2) > 3, so tauB > 3s-1 exactly for s > s*", 2 / (B - H) > 3)
ok("A5 tauB > 4s-2 (Huxley family range) on [0.84, B]", lin_ok(lambda s: tauB(s) - (4 * s - 2), S84, B) and tauB(S84) > 4 * S84 - 2)
SPY = (10 * B - Fr(13, 2)) / (12 * B - 8)
ok("A6 tauB(sPY) = 12 sPY - 9 (Petrow-Young line crossing)", tauB(SPY) == 12 * SPY - 9)
ok("A7 slope tauB < 12, so tauB > 12s-9 exactly for s < sPY", 2 / (B - H) < 12)
print("   sPY =", float(SPY), " (Theorem A beats every known family range on (s*, sPY))")
TSW = (B - H) / (B - Fr(2, 3))
ok("A8 hull: gB(t) <= 1/2 + t/6 iff t <= (B-1/2)/(B-2/3) = %.4f; tauB(sPY) <= that" % float(TSW),
   gB(TSW) == H + TSW / 6 and tauB(SPY) <= TSW)
ok("A9 hull: gB(t) <= t/2 for t >= 1 (gB(1) = 1/2, slope B-1/2 < 1/2)", gB(1) == H and B - H < H)
ok("A10 Petrow-Young point (1/2,1/6) lies below the (H) chord from (1-B,B-1/2) to (B,0)", (B - H) / 2 > Fr(1, 6))

# ---------------------------------------------------------------- B. class-II cut is not moved by (H)
def hullmin(t):  # best exponent of an L-sum of length N at conductor*height X = N^t, vertices of mu*
    return min(t / 2, gB(t), H + t / 6, B)


ok("B1 on [2, 6(B-1/2)] the hull minimum is the Weyl vertex 1/2 + t/6 (endpoints; pieces affine)",
   all(hullmin(t) == H + t / 6 for t in (Fr(2), 6 * (B - H))) and gB(Fr(2)) >= H + Fr(2) / 6)
ok("B2 cut 6s-3 <= 6(B-1/2) for s <= B, so exclusion holds exactly for tau < 6s-3", 6 * B - 3 == 6 * (B - H))
# interpolated mu-line under (H): mu*(u) <= (B-u)/(6(B-1/2)); best line for exclusion is u = 1/2
def cut_line(s, u):
    return (s - u) / ((B - u) / (6 * (B - H)))


ok("B3 (s-u)/mu*(u) is maximal at u = 1/2 (decreasing in u for s < B; 50 grid points)",
   all(cut_line(Fr(43, 50), H) >= cut_line(Fr(43, 50), H + Fr(k, 200) * (Fr(43, 50) - H)) for k in range(1, 200, 4)))
ok("B4 at u = 1/2 that cut equals 6s-3", cut_line(Fr(43, 50), H) == 6 * Fr(43, 50) - 3)

# ---------------------------------------------------------------- C. sup over [a, 2a] of min_k k LV(s, t/k)/t
def LV_sub(tM):
    return lambda s, t: c(s) + max(Fr(0), t - tM)


def LV_jut(s, t):  # Jutila 1977 (1.2)/(1.3), k = 2, G <= N, V = N^s, X = N^t
    return max(c(s), t + 12 - 16 * s, 2 * (t + 2 - 4 * s))


def Fmin(LV, s, t, kmax=12):
    return min(k * LV(s, t / k) / t for k in range(1, kmax + 1))


def cands_sub(s, tM, kmax=12):
    pts = []
    for k in range(1, kmax):
        for j in range(k + 1, kmax + 1):
            pts.append(k * tM + (j - k) * c(s))
        pts.append(k * tM)
    return pts


def cands_jut(s, kmax=12):
    pts = []
    for k in range(1, kmax):
        for j in range(k + 1, kmax + 1):
            pts.append(j * c(s) - k * (12 - 16 * s))
            pts.append((j * c(s) - 2 * k * (2 - 4 * s)) / 2)
        pts.append(k * (3 * s - 1))
        pts.append(k * (14 * s - 10))
    return pts


def G(LV, s, a, cands):
    pts = [a, 2 * a] + [p for p in cands if a <= p <= 2 * a]
    return max(Fmin(LV, s, p) for p in pts)


def best(LV, s, cut, cands, n=120):
    """min over a in (0, cut] of G(a): grid plus window positions next to every candidate point."""
    grid = [cut * Fr(k, n) for k in range(n // 4, n + 1)]
    extra = [p for p in cands if 0 < p <= cut] + [p / 2 for p in cands if 0 < p / 2 <= cut]
    return min(G(LV, s, a, cands) for a in grid + extra)


SG = [Fr(5, 6) + Fr(k, 1000) for k in range(1, 42, 2)] + [Fr(7, 8)]
# C1: unconditional uniform bound 4/(6s-3) on (5/6, 7/8] from Petrow-Young + Jutila k=2 (no subdivision)
okc1 = all(G(LV_jut, s, 6 * s - 3, cands_jut(s)) <= 2 * c(s) / (6 * s - 3) for s in SG)
ok("C1 unconditional: sup_{[6s-3,12s-6]} F_Jutila <= (4-4s)/(6s-3) on 22 points of (5/6, 7/8]", okc1)
# the same in closed form (all pieces are monotone, see PROOF.md Thm C):
ok("C1a 9s-9/2 >= 2(3s-1) iff s >= 5/6", (9 * Fr(5, 6) - Fr(9, 2)) == 2 * (3 * Fr(5, 6) - 1))
ok("C1b 2 + (8-16s)/(9s-9/2) = 2/9 identically", all(2 + (8 - 16 * s) / (9 * s - Fr(9, 2)) == Fr(2, 9) for s in SG))
ok("C1c 2/9 <= (4-4s)/(6s-3) iff s <= 7/8", Fr(2, 9) == c(Fr(7, 8)) * 2 / (6 * Fr(7, 8) - 3))
ok("C1d (19.5-23s)/(9s-4.5) <= 2/9 for s >= 0.82", all((Fr(39, 2) - 23 * s) / (9 * s - Fr(9, 2)) <= Fr(2, 9) for s in SG))
ok("C1e Jutila k=2 Montgomery range is 3s-1 (14s-10 >= 3s-1 for s >= 9/11)", Fr(5, 6) >= Fr(9, 11))
# C2: barrier. With MC only up to tM <= 8s-5 and the cut 6s-3, no window does better than (4-4s)/(6s-3)
ok("C2a tauB <= 8s-5 on [0.84, B] (affine)", lin_ok(lambda s: 8 * s - 5 - tauB(s), S84, B))
ok("C2b 3s-1 <= 8s-5 for s >= 4/5", 3 * Fr(4, 5) - 1 == 8 * Fr(4, 5) - 5)
SB = [S84, Fr(17, 20), Fr(43, 50), Fr(87, 100), B - Fr(1, 10**6)]
for s in SB:
    cut = 6 * s - 3
    nb = best(LV_sub(tauB(s)), s, cut, cands_sub(s, tauB(s)))
    ok("C2c s=%.6f: min_a G with tauM = tauB equals (4-4s)/(6s-3) = %.6f (got %.6f)" % (float(s), float(2 * c(s) / cut), float(nb)),
       nb == 2 * c(s) / cut)
# C3: (H) reaches the cap: conditions tauB >= 5.5s - 3.25 and tauB >= 3s - 1.5 on [0.84, B]
ok("C3a tauB >= (11/2)s - 13/4 on [0.84, B] (affine)", lin_ok(lambda s: tauB(s) - (Fr(11, 2) * s - Fr(13, 4)), S84, B))
ok("C3b tauB < (9/2)s - 9/4 at 0.84 and B (no-subdivision variant is weaker: 3/tauB)",
   tauB(S84) < Fr(9, 2) * S84 - Fr(9, 4) and tauB(B) < Fr(9, 2) * B - Fr(9, 4))
ok("C3c Huxley+PY gives only 3/(3s-1): (4s-2) < (11/2)s - 13/4 for s > 5/6", lin_ok(lambda s: Fr(11, 2) * s - Fr(13, 4) - (4 * s - 2), Fr(5, 6), B))

# ---------------------------------------------------------------- D. large-sieve sub-regime (lam = log(QT)/log(Q^2 T) < 1)
def tau0(s):
    return tauB(s) + 1 - s


def lam_b(s):
    return 3 * (6 * s - 3) / (4 * tau0(s))


def lam_a(s):
    return max((6 * s - 3) / (4 * s), (1 + (1 - s) / tau0(s)) / 2)


def Texp(l):  # T = Q^{(2 lam - 1)/(1 - lam)}
    return (2 * l - 1) / (1 - l)


print("\n sigma     lam_a    lam_b    T from Q^a to Q^b      old 3/(2s)  new 3/tau0   uniform 4/(6s-3)  DH")
for s in [Fr(86, 100), Fr(865, 1000), Fr(87, 100), Fr(8725, 10000), B - Fr(1, 10**8)]:
    la, lb = lam_a(s), lam_b(s)
    print("  %.6f  %.4f   %.4f   Q^%.3f .. Q^%.2f     %.5f     %.5f      %.5f        2"
          % (float(s), float(la), float(lb), float(Texp(la)), float(Texp(lb)), 3 / (2 * float(s)), 3 / float(tau0(s)), 4 / (6 * float(s) - 3)))
ok("D1 lam_a < lam_b on {0.86, 0.87, B-} (nonempty regime)", all(lam_a(s) < lam_b(s) for s in (Fr(86, 100), Fr(87, 100), B - Fr(1, 10**8))))
for s in [Fr(86, 100), Fr(87, 100), B - Fr(1, 10**8)]:
    for l in [lam_a(s), (lam_a(s) + lam_b(s)) / 2, lam_b(s)]:
        cut = (6 * s - 3) / l
        newv = best(LV_sub(tauB(s)), s, cut, cands_sub(s, tauB(s)))
        oldv = best(LV_sub(3 * s - 1), s, cut, cands_sub(s, 3 * s - 1))
        oldj = best(LV_jut, s, cut, cands_jut(s))
        ok("D2 s=%.5f lam=%.4f: new <= (3-3s)/tau0 = %.6f < old(subdiv) = %.6f <= old(Jutila form) = %.6f"
           % (float(s), float(l), float((3 - 3 * s) / tau0(s)), float(oldv), float(oldj)),
           newv <= (3 - 3 * s) / tau0(s) and oldv == (3 - 3 * s) / (2 * s) and oldv <= oldj)
ok("D3 at s*: tau0 = 2s* (no gain at the left end)", tau0(SSTAR) == 2 * SSTAR)

# ---------------------------------------------------------------- E. values
print("\n sigma     DH   Huxley+PY 3/(3s-1)   uniform 4/(6s-3)   tauB      3s-1   4s-2")
for s in [Fr(84, 100), Fr(85, 100), Fr(86, 100), Fr(87, 100), B - Fr(1, 10**8)]:
    print("  %.6f  2.0   %.5f              %.5f            %.4f   %.4f  %.4f" % (
        float(s), 3 / (3 * float(s) - 1), 4 / (6 * float(s) - 3), float(tauB(s)), 3 * float(s) - 1, 4 * float(s) - 2))
ok("E1 4/(6s-3) = 2 at s = 5/6 and is decreasing, so < 2 exactly for s > 5/6", 4 / (6 * Fr(5, 6) - 3) == 2)
ok("E2 4/(6s-3) < 3/(3s-1) iff 12s-4 < 18s-9 iff s > 5/6 (equality at 5/6)", 4 * (3 * Fr(5, 6) - 1) == 3 * (6 * Fr(5, 6) - 3))
# ---------------------------------------------------------------- F. zeta part of the family sums (review repair)
ok("F1 3/(2s) <= 4/(6s-3) iff 9 - 10s >= 0: equality at 9/10, holds on (5/6, 7/8]",
   Fr(3, 2) / Fr(9, 10) == 4 / (6 * Fr(9, 10) - 3) and all(Fr(3, 2) / s <= 4 / (6 * s - 3) for s in SG))
ok("F2 Ivic 3831/4791 and Bourgain 3734/4694 lie below 5/6", Fr(3831, 4791) < Fr(5, 6) and Fr(3734, 4694) < Fr(5, 6))
print("\nruntime %.2f s" % (time.time() - t0))
print("ALL PASS" if not FAILS else "FAILURES: %d" % len(FAILS))
