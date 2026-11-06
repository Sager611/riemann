#!/usr/bin/env python3
"""exp48 certificate (python-flint arb, 256-bit balls). Run:
  uv run --no-project --with python-flint --with mpmath python3 cert.py
Every inequality is decided on arb balls (rigorous); a FAIL aborts with exit code 1.
Notation: d = |D|, D < -4 a fundamental discriminant; h = h(D); omega = number of prime factors of D;
  L78(d)   = 0.0350 sqrt(d)/max(3, log log d)   (Littlewood78, Cor 2: h >= L78(d), every D < -4)
  Chen1(d) = sqrt(d)/(7.732 pi log d)            (EKN Lemma Chen, first term of (b_ch))
  Sieg(d)  = 0.655/(pi e) sqrt(d)/log d          (EKN Lemma Siegel)
  d_k      = product of the first k primes       (EKN Lemma lem7: d >= d_omega)
"""
import sys
from flint import arb, ctx
ctx.prec = 256
OK = True
def check(name, cond):
    global OK
    print(("PASS " if cond else "FAIL ") + name)
    if not cond: OK = False

def primes(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]
P = primes(20000)
PI, E1, LOG2 = arb.pi(), arb(1).exp(), arb(2).log()
C78 = arb("0.0350")
CCHEN = 1 / (arb("7.732") * PI)
CSIEG = arb("0.655") / (PI * E1)
def M(lx):  # max(3, log log x) given lx = log x (arb); decided rigorously
    ll = lx.log()
    if ll > 3: return ll
    if ll < 3: return arb(3)
    raise ValueError("undecided max at log log x = 3")

print("== (0) constants")
print("  C78 = 0.0350; 0.11/pi =", (arb("0.11") / PI).str(12), "; 1/(7.732 pi) =", CCHEN.str(12), "; 0.655/(pi e) =", CSIEG.str(12))
check("(0a) 0.11/pi >= 0.0350 (Littlewood78 Thm 1(1) -> Cor 2 constant)", arb("0.11") / PI > C78)
check("(0b) 1/(7.732 pi) < 0.0412 (brief's rounded Chen constant is an upper rounding)", CCHEN < arb("0.0412"))
check("(0c) 0.0767 < 0.655/(pi e) < 0.0768 (the brief's 0.0767 is EKN's constant rounded down; comparisons below use the exact ball)", arb("0.0767") < CSIEG < arb("0.0768"))

print("== (1) L78(d) > Chen1(d) for every d >= 35  [so (b_ch) holds for every D with |D| >= 10^6]")
# (1a) 35 <= d <= e^(e^3): max = 3; need 0.0350 log d > 3 CCHEN; LHS increasing in d -> check at d = 35
check("(1a) 0.0350 log 35 > 3/(7.732 pi)", C78 * arb(35).log() > 3 * CCHEN)
print("     3/(7.732 pi) =", (3 * CCHEN).str(6), "; 0.0350 log 35 =", (C78 * arb(35).log()).str(6), "; 0.0350 e^3/3 =", (C78 * arb(3).exp() / 3).str(4))
print("     margin at d=35:", (C78 * arb(35).log() - 3 * CCHEN).str(6), "; at d = 34:", (C78 * arb(34).log() - 3 * CCHEN).str(6), "(negative: 35 is the first integer)")
# (1b) u = log d >= e^3: need 0.0350 u/log u > CCHEN; u/log u increasing for u >= e -> check at u = e^3
u = arb(3).exp()
check("(1b) 0.0350 e^3/3 > 1/(7.732 pi)", C78 * u / 3 > CCHEN)
# (1c) the brief's form with 0.0412, same split
check("(1c) 0.0350 log 35 > 3*0.0412 and 0.0350 e^3/3 > 0.0412", C78 * arb(35).log() > 3 * arb("0.0412") and C78 * u / 3 > arb("0.0412"))
# ratio L78/Chen1 at d = 10^6 and its minimum over d >= 10^6 (attained at d = 10^6 on the max=3 branch, increasing after)
r6 = C78 * arb(10**6).log() / (3 * CCHEN)
check("(1d) L78/Chen1 >= 3.9 for every d >= 10^6 (min at d = 10^6: %s; on log log d >= 3 the ratio is 0.0350 u/(log u * CCHEN) >= %s)" % (r6.str(5), (C78 * u / (3 * CCHEN)).str(5)), r6 > arb("3.9") and C78 * u / (3 * CCHEN) > arb("3.9"))

print("== (2) L78(d) > Sieg(d) for every d >= e^6.6  [so EKN's faster 'no Siegel zero' run is also complete]")
check("(2a) 0.0350 * 6.6 > 3 * 0.655/(pi e)", C78 * arb("6.6") > 3 * CSIEG)
check("(2b) 0.0350 e^3/3 > 0.655/(pi e)", C78 * u / 3 > CSIEG)
check("(2c) e^6.6 < 10^6 and e^11.2 < 10^6 (EKN's ranges start above both)", arb("6.6").exp() < 10**6 and arb("11.2").exp() < 10**6)
print("     0.0350*6.6 =", (C78 * arb("6.6")).str(4), "; 3*0.655/(pi e) =", (3 * CSIEG).str(6))
B = arb("0.87338196661320752")
check("(2d) half-plane gives Lemma Siegel's hypothesis: 1 - 1/(4 log|D|) > B iff log|D| > 1/(4(1-B)) = %s < 11.2 (7/8: 2)" % (1 / (4 * (1 - B))).str(5), 1 / (4 * (1 - B)) < arb("1.975") and 1 / (4 * (1 - arb(7) / 8)) <= 2)

print("== (3) number of prime factors: L78(d_k) > c^(k-1) for all k >= K_c, so E(D) | c forces omega(D) <= K_c - 1")
ld = [arb(0)]
for p in P[:400]: ld.append(ld[-1] + arb(p).log())
def logL78(lx): return C78.log() + lx / 2 - M(lx).log()
RAT = 1 + LOG2 / 3   # bound for M(d_{k+1})/M(d_k), since p_{k+1} < d_k for k >= 2
Kc = {}
for c in (2, 4, 8):
    lc = arb(c).log()
    ok = [logL78(ld[k]) > (k - 1) * lc for k in range(1, 400)]
    K = next(k for k in range(1, 399) if all(ok[k - 1:]))
    Kc[c] = K
    fails_before = not ok[K - 2]
    check("(3) c=%d: K_c = %d; L78(d_k) > c^(k-1) for %d <= k <= 399 (exact d_k), fails at k = %d" % (c, K, K, K - 1), all(ok[K - 1:]) and fails_before)
    # induction for k >= 399: L78(d_{k+1})/L78(d_k) >= sqrt(p_{k+1})/RAT > c
    check("(3) c=%d: induction step sqrt(p_{k+1})/(1 + log2/3) > c for all k >= %d (p_{k+1} >= %d)" % (c, K, P[K]), arb(P[K]).sqrt() / RAT > c)
print("  K_2, K_4, K_8 =", Kc[2], Kc[4], Kc[8], " (EKN's N_2, N_4, N_8 = 11, 24, 58)")
check("(3) K_8 <= 58 and K_4 <= 24 and K_2 <= 11: every omega left open by L78 lies inside EKN's search range (k <= N_c)", Kc[8] <= 58 and Kc[4] <= 24 and Kc[2] <= 11)

print("== (3') side facts used in (3)")
dk = 1; okE = True
for k in range(1, 399):
    dk *= P[k - 1]
    if k >= 2 and not P[k] < dk: okE = False
check("(3') p_{k+1} < d_k for 2 <= k <= 398 (Euclid: d_k - 1 >= 2 has a prime factor > p_k, so p_{k+1} <= d_k - 1 for all k >= 2)", okE)
# L78 increasing in x: branch log log x < 3 is 0.0350 sqrt(x)/3; branch u = log x > e^3: d/du (u/2 - log log u) = 1/2 - 1/(u log u) > 0
check("(3') 1/2 - 1/(u log u) > 0 at u = e^3 (and u log u increases), so L78 is increasing on [e^(e^3), oo)", arb(1) / 2 - 1 / (u * 3) > 0)

print("== (4) EKN Theorem upp_bound recomputed (three conditions; constant 0.655/(pi e)) and with Chen's constant")
def ekn_N(r, cst):
    for n in range(1, 400):
        lx = ld[n]
        c1 = lx > arb("11.2")
        c2 = (arb(1) / 2 - 1 / lx) * arb(P[n - 1]).log() > r * LOG2
        c3 = cst.log() + lx / 2 - lx.log() > r * (n - 1) * LOG2
        if c1 and c2 and c3: return n
for r, Nexp, dmax in ((1, 11, "2.01e11"), (2, 24, "2.38e34"), (3, 58, "3.17e110")):
    N = ekn_N(r, CSIEG)
    dN = 1
    for p in P[:N]: dN *= p
    check("(4) c=%d: N = %d (EKN: %d), d_N = %.4e <= %s" % (2**r, N, Nexp, float(dN), dmax), N == Nexp and dN <= float(dmax))
for r in (1, 2, 3):
    print("     with Chen's constant 1/(7.732 pi) instead: N_%d = %d" % (2**r, ekn_N(r, CCHEN)))

print("== (5) fallback: OpenAI's 7/8 half-plane instead of B (Littlewood78 cert_output_7_8: |L(1,chi)| >= 0.105/max(3, log log q))")
C7 = arb("0.0334")
check("(5a) 0.105/pi >= 0.0334", arb("0.105") / PI > C7)
check("(5b) 0.0334 log 41 > 3/(7.732 pi) and 0.0334 e^3/3 > 1/(7.732 pi): L78-fallback > Chen1 for d >= 41", C7 * arb(41).log() > 3 * CCHEN and C7 * u / 3 > CCHEN)
def logL7(lx): return C7.log() + lx / 2 - M(lx).log()
for c in (2, 4, 8):
    lc = arb(c).log()
    ok = [logL7(ld[k]) > (k - 1) * lc for k in range(1, 400)]
    K = next(k for k in range(1, 399) if all(ok[k - 1:]))
    check("(5c) c=%d: fallback K_c = %d <= EKN's N_c, induction holds (sqrt(p_{K+1})/(1+log2/3) > c)" % (c, K), K <= {2: 11, 4: 24, 8: 58}[c] and arb(P[K]).sqrt() / RAT > c)

print("ALL PASS" if OK else "SOME FAIL")
sys.exit(0 if OK else 1)
