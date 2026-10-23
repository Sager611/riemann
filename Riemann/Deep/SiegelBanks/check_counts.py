#!/usr/bin/env python3
"""Exact checks for Banks, arXiv 2610.31402 (Lemma count, Prop main bookkeeping,
parameter arithmetic of Thm quadratic / Thm main / Cor siegel, toy run of Lemma pigeon).
Pure Python, < 50 MB, < 30 s.  Run: python3 check_counts.py > check_counts_output.txt"""
from fractions import Fraction as Fr
from math import comb, log, ceil, isqrt, factorial

def sums(s, H, D, ps=()):
    """Exact A, W, sum alpha_1, E_p over idx = {alpha in N^s : alpha_1 + H|alpha*| <= D}."""
    A = W1 = Ws = 0
    E = {p: 0 for p in ps}
    for j in range(D // H + 1):
        m = 1 if s == 1 else comb(j + s - 2, s - 2)
        if s == 1 and j > 0:
            break
        L = D - H * j
        A += m * (L + 1)
        W1 += m * L * (L + 1) // 2
        Ws += m * j * (L + 1)
        for p in ps:
            q, r = divmod(L, p)
            E[p] += m * (p * q * (q - 1) // 2 + q * (r + 1))
    return A, W1, Ws, E

def I_of(s, H, D):
    return Fr(D ** (s + 1), factorial(s + 1) * H ** (s - 1))

print("== 1. Lemma count: exact sums vs simplex integrals")
worst = 0.0
for (s, H) in [(1, 2), (2, 2), (3, 8), (3, 16), (4, 5), (5, 3)]:
    for D in [10**3, 10**4, 10**5]:
        ps = sorted({H, 2 * H + 1, isqrt(D), D // 3, D // 2, D})
        A, W1, Ws, E = sums(s, H, D, ps)
        I = I_of(s, H, D)
        eA = float((A - (s + 1) * I / D) / Fr(D ** (s - 1)))
        eW = float((W1 + Ws - (1 + Fr(s - 1, H)) * I) / Fr(D ** s))
        e1 = float((W1 - I) / Fr(D ** s))
        eE = max(abs(float((E[p] - I / p) / Fr(D ** s))) for p in ps)
        # the inequality actually used in Prop main: E_p >= (sum alpha_1)/p - A
        ok = all(E[p] * p >= W1 - A * p for p in ps)
        worst = max(worst, abs(eA), abs(eW), eE)
        print(f"s={s} H={H} D={D:>6}: (A-(s+1)I/D)/D^(s-1)={eA:+.4f}  (W-(1+(s-1)/H)I)/D^s={eW:+.4f}"
              f"  (sum a1 - I)/D^s={e1:+.4f}  max|E_p-I/p|/D^s={eE:.4f}  E_p>=sum/p-A:{ok}")
print("max normalized error", round(worst, 4), "(bounded => O(D^(s-1)), O(D^s) as claimed)")

print("\n== 2. Kerr special case s=3, H=8 (Kerr: r=D^3/384, W=5D^4/6144, E_p=D^4/(1536p))")
D = 10**5
A, W1, Ws, E = sums(3, 8, D, [101, 10007])
print(" A/(D^3/384) =", float(Fr(A) / Fr(D**3, 384)), " W/(5D^4/6144) =", float(Fr(W1 + Ws) / Fr(5 * D**4, 6144)),
      " E_101/(D^4/(1536*101)) =", float(Fr(E[101]) / Fr(D**4, 1536 * 101)))
print(" I = D^4/1536 for (s,H)=(3,8):", I_of(3, 8, D) == Fr(D**4, 1536))

print("\n== 3. Parameter arithmetic")
def is_prime(n):
    if n < 2: return False
    i = 2
    while i * i <= n:
        if n % i == 0: return False
        i += 1
    return True
bad = 0; cnt = 0
for e_num in range(1, 101):
    eps = Fr(e_num, 100)
    for nk in range(1, 13):
        ell = 2
        while not (is_prime(ell) and ell - 1 >= 8 * nk / eps):
            ell += 1
        H = ceil(16 * nk * (ell - 1) / eps)
        C0 = ceil(4 / eps)
        # |G| = [K:k] is even, (ell-1)/nk <= |G| <= 2(ell-1)
        for G in range(2, 2 * (ell - 1) + 1, 2):
            if G * nk < ell - 1:
                continue
            s = G // 2 + 1
            cnt += 1
            c1 = Fr(s - 1, H) <= eps / 16
            c2 = Fr(s, G) <= Fr(1, 2) + eps / 8
            c3 = Fr(1, 2 * C0) <= eps / 8
            coef = (1 + Fr(s - 1, H)) * (Fr(s, G) + Fr(1, 2 * C0))
            if not (c1 and c2 and c3 and coef < Fr(1, 2) + eps):
                bad += 1
print(" Thm quadratic: checked", cnt, "(eps, n_k, |G|) triples; failures:", bad)
bad = 0; cnt = 0
for e_num in range(1, 101):
    eps = Fr(e_num, 100)
    for nK in range(1, 61):
        H = ceil(4 * nK / eps); C0 = ceil(nK / eps) + 1
        for nk in [d for d in range(1, nK + 1) if nK % d == 0]:
            G = nK // nk
            for s in {1, max(1, G // 2), G}:
                cnt += 1
                coef = (1 + Fr(nK, H)) * (Fr(s, G) + Fr(nK, 2 * C0))
                if not coef < Fr(s, G) + eps:
                    bad += 1
print(" Thm main: checked", cnt, "cases; failures:", bad)
eps = Fr(1, 8)
print(" Cor siegel: C0(1/8) =", ceil(4 / eps), "; split >= 1 - (1/2+1/8) - 1/8 =", 1 - (Fr(1, 2) + eps) - Fr(1, 8),
      "; minus c1 log d_F <= 1/8 leaves", 1 - (Fr(1, 2) + eps) - Fr(1, 8) - Fr(1, 8), "> 0")

print("\n== 4. Toy run of Lemma pigeon (must FAIL for every N, else Masser is contradicted)")
# K = Q(sqrt5, sqrt-3), k = Q, G = {1, s, t, st}, F = Q(sqrt5), S = {1, s, st}; delta_p = 1 iff
# p splits completely (p = +-1 mod 5 and p = 1 mod 3) or p is inert in F (p = +-2 mod 5).
# omega = 1, phi, w, phi*w (phi = golden ratio, w = (1+sqrt-3)/2): house <= phi.
def primes_upto(n):
    sv = bytearray([1]) * (n + 1); sv[0:2] = b"\x00\x00"
    for i in range(2, isqrt(n) + 1):
        if sv[i]: sv[i*i::i] = bytearray(len(sv[i*i::i]))
    return [i for i in range(n + 1) if sv[i]]
lam = (1 + 5 ** 0.5) / 2; nK = 4; s = 3
for H in [8, 200]:
    worst = 0.0
    for N in [10, 30, 100, 300, 1000, 3000]:
        D = ceil(((N - 1) / s) ** (nK / s)) - 1
        P = [p for p in primes_upto(D) if p >= H and p not in (3, 5)
             and ((p % 5 in (2, 3)) or (p % 5 in (1, 4) and p % 3 == 1))]
        A, W1, Ws, E = sums(s, H, D, P)
        lhs = sum(E[p] * log(p) for p in P)            # (1/n_k) sum delta_p E_p log p, n_k = 1
        B = nK * lam * N
        rhs = (W1 + Ws) * log(B) + A * log(4 * nK * A * lam)
        worst = max(worst, lhs / rhs)
        print(f" H={H:>3} N={N:>5} D={D:>7} #primes={len(P):>5}  lhs/rhs={lhs/rhs:.4f}")
    print(f" H={H}: max lhs/rhs = {worst:.4f} (< 1: no contradiction with Masser)")

