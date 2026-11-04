#!/usr/bin/env python3
"""Check of the preprint's Lemma 'Character comparison' (Section 4).

Claim: A = Z[eta] localized at a prime pi over p (p odd, p not | ord(eta)); F, K in A[C_{p^e}],
F*K = p^e * b with b a unit of A.  Then F(rho)/F(1) is a unit with residue 1 mod (rho-1),
rho a primitive p-th root of unity.

(a) exhaustive: every F in Z[C_9] with F F^* = 9.
(b) random: F = U*W*X^a in Z[C_{p^e}] and Z[i][C_{p^e}] (U a random unit, W a random product of divisors
    of p^e); K = p^e F^{-1} solved exactly and checked A-integral; the ratio test must pass for F and for K.
(c) sharpness: conclusion FAILS for p - 2J in Z[C_p] (scalar exponent 2 > group exponent 1,
    the paper's Remark) and for 3 - 2J_H in Z[C_3 x C_3] (non-cyclic 3-group).
Exact integer / Fraction arithmetic only.  Usage: python3 check_comparison.py [samples]
"""
import itertools, random, sys
from fractions import Fraction

def gmul(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gadd(a, b): return (a[0]+b[0], a[1]+b[1])
def gsub(a, b): return (a[0]-b[0], a[1]-b[1])
Z0 = (0, 0)

def vint(n, p):
    if n == 0: return 10**9
    v = 0
    while n % p == 0: n //= p; v += 1
    return v

def vpi(a, p):
    """valuation of Gaussian integer a at the prime pi over p used below:
    p = 3 mod 4: pi = p (inert); p = 1 mod 4: pi = the factor with i = r mod pi, r^2 = -1 mod p, r the smallest."""
    if a == Z0: return 10**9
    if p % 4 == 3:
        return min(vint(a[0], p), vint(a[1], p))
    r = min(x for x in range(p) if (x*x+1) % p == 0)
    # pi = p/gcd stuff: pi | a  <=>  a0 + r*a1 = 0 mod p for pi = (r - i)?  use explicit pi below
    pi = PI[p]; pib = (pi[0], -pi[1]); npi = pi[0]**2 + pi[1]**2
    v = 0
    while True:
        t = gmul(a, pib)  # a / pi = a * conj(pi) / p
        if t[0] % npi or t[1] % npi: return v
        a = (t[0]//npi, t[1]//npi); v += 1

PI = {5: (2, 1), 13: (3, 2)}
REP = None

def eval_mod_p(F, N, p):
    """F in Z[i][C_N] (list of Gaussian ints), p | N: image under X -> zeta_p, as list of length p (mod x^p - 1)."""
    z = [Z0]*p
    for j, f in enumerate(F): z[j % p] = gadd(z[j % p], f)
    return z

def pmul(a, b, p):
    c = [Z0]*p
    for i, x in enumerate(a):
        if x == Z0: continue
        for j, y in enumerate(b):
            if y == Z0: continue
            c[(i+j) % p] = gadd(c[(i+j) % p], gmul(x, y))
    return c

def relnorm(z, p):
    """N_{Q(i,zeta_p)/Q(i)} of z(zeta_p), z given mod x^p-1; returns Gaussian integer."""
    r = [Z0]*p; r[0] = (1, 0)
    for a in range(1, p):
        za = [Z0]*p
        for k, x in enumerate(z): za[(a*k) % p] = gadd(za[(a*k) % p], x)
        r = pmul(r, za, p)
    # reduce mod Phi_p: constant must remain
    top = r[p-1]
    red = [gsub(r[k], top) for k in range(p-1)]
    assert all(t == Z0 for t in red[1:]), red
    return red[0]

def v_elem(z, p):
    """normalized valuation (v(p)=1) at the unique prime over pi in Z[i][zeta_p]."""
    zz = [gsub(x, z[p-1]) for x in z]  # reduce mod Phi_p for zero test
    if all(t == Z0 for t in zz): return 10**9
    return Fraction(vpi(relnorm(z, p), p), p-1)

def ratio_ok(F, N, p):
    """F(rho)/F(1) unit with residue 1 mod (rho-1)  <=>  v(F(rho)-F(1)) > v(F(1))."""
    z = eval_mod_p(F, N, p)
    one = Z0
    for f in F: one = gadd(one, f)
    d = list(z); d[0] = gsub(d[0], one)
    return v_elem(d, p) > vpi(one, p), z, one

def conv(F, K, N):
    c = [Z0]*N
    for i, x in enumerate(F):
        if x == Z0: continue
        for j, y in enumerate(K):
            if y == Z0: continue
            c[(i+j) % N] = gadd(c[(i+j) % N], gmul(x, y))
    return c

def star(F, N):
    return [(F[(-j) % N][0], -F[(-j) % N][1]) for j in range(N)]

# ---------- (a) exhaustive F in Z[C_9], F F^* = 9 ----------
def vectors_sumsq(n, s):
    if n == 0:
        if s == 0: yield ()
        return
    m = int(s**0.5)
    for x in range(-m, m+1):
        if x*x <= s:
            for rest in vectors_sumsq(n-1, s-x*x): yield (x,) + rest

def part_a():
    N, p = 9, 3
    sols = 0; bad = 0; kinds = set()
    for v in vectors_sumsq(N, 9):
        F = [(x, 0) for x in v]
        P = conv(F, star(F, N), N)
        if P[0] == (9, 0) and all(t == Z0 for t in P[1:]):
            sols += 1
            ok, z, one = ratio_ok(F, N, p)
            kinds.add(tuple(sorted(abs(x) for x in v if x)))
            global REP
            if len([x for x in v if x]) == 6: REP = v
            if not ok: bad += 1
    print(f"(a) Z[C_9], FF*=9: {sols} solutions, supports {sorted(kinds)}, ratio-test failures: {bad}")
    print(f"    sample non-rational solution: {REP}")
    return bad == 0

# ---------- (b) random F dividing p^e ----------
def solve_inverse(F, N, pe):
    """K = pe * F^{-1} in Q(i)[C_N] via exact elimination on the 2N x 2N real system; None if singular."""
    n2 = 2*N
    M = [[Fraction(0)]*(n2+1) for _ in range(n2)]
    # unknown K = (kr, ki); (F*K)_g = sum_h F_h K_{g-h}
    for g in range(N):
        for h in range(N):
            fr, fi = F[h]; col = (g-h) % N
            M[g][col] += fr; M[g][N+col] -= fi          # real part
            M[N+g][col] += fi; M[N+g][N+col] += fr      # imag part
    M[0][n2] = Fraction(pe)
    # Gaussian elimination
    r = 0
    for c in range(n2):
        piv = next((i for i in range(r, n2) if M[i][c] != 0), None)
        if piv is None: return None
        M[r], M[piv] = M[piv], M[r]
        inv = 1/M[r][c]
        M[r] = [x*inv for x in M[r]]
        for i in range(n2):
            if i != r and M[i][c] != 0:
                f = M[i][c]; Mi = M[i]; Mr = M[r]
                M[i] = [a - f*b for a, b in zip(Mi, Mr)]
        r += 1
    sol = [M[i][n2] for i in range(n2)]
    return [(sol[g], sol[N+g]) for g in range(N)]

def integral_at_pi(K, p):
    for a, b in K:
        D = a.denominator * b.denominator // __import__('math').gcd(a.denominator, b.denominator)
        num = (int(a*D), int(b*D))
        if vpi(num, p) < vint(D, p) * (1 if p % 4 == 3 else 1): return False
    return True

def to_int(K):
    D = 1
    import math
    for a, b in K:
        D = D*a.denominator//math.gcd(D, a.denominator); D = D*b.denominator//math.gcd(D, b.denominator)
    return [(int(a*D), int(b*D)) for a, b in K], D

def part_b(samples, seed=1):
    """F = U * W * X^a with U a random unit of A[C_N] (augmentation prime to pi), W a random product of
    known divisors of p^e (p^m, p^j - 2 J_H with 2j <= e, pi^m for split p).  K = p^e F^{-1} is solved
    exactly; we keep F only if K is pi-integral and then test the ratio for F and for K."""
    rng = random.Random(seed)
    allok = True
    for (p, e, gauss, ns) in [(3, 2, False, samples), (3, 3, False, samples), (5, 2, False, samples),
                              (7, 2, False, samples // 2), (3, 2, True, samples // 2), (5, 2, True, samples // 3)]:
        N = p**e; pe = p**e
        tested = fails = 0; vals = set(); nonrat = 0
        for _ in range(ns):
            # unit U: random small coefficients, augmentation prime to pi
            while True:
                U = [(rng.randint(-1, 1), rng.randint(-1, 1) if gauss else 0) for _ in range(N)]
                aug = (sum(u[0] for u in U), sum(u[1] for u in U))
                if vpi(aug, p) == 0: break
            # divisor W of p^e
            W = [Z0]*N; W[rng.randrange(N)] = (1, 0)
            budget = e
            while budget > 0:
                choice = rng.random()
                if choice < 0.4 and 2 <= budget:              # p^j - 2 J_H with j = 1 (needs exponent 2)
                    H = [k for k in range(N) if k % (N // p) == 0]   # subgroup of order p
                    D = [Z0]*N; D[0] = (p, 0)
                    for k in H: D[k] = gsub(D[k], (2, 0))
                    W = conv(W, D, N); budget -= 2
                elif choice < 0.6 and gauss and p % 4 == 1:   # pi itself (pi-bar is a unit in A)
                    W = conv(W, [PI[p]] + [Z0]*(N-1), N); budget -= 1
                elif choice < 0.8:
                    W = [(p*a, p*b) for a, b in W]; budget -= 1
                else:
                    break
            F = conv(U, W, N)
            K = solve_inverse(F, N, pe)
            if K is None or not integral_at_pi(K, p): continue
            Ki, D = to_int(K)
            tested += 1
            one = (sum(f[0] for f in F), sum(f[1] for f in F))
            vals.add(vpi(one, p))
            z = eval_mod_p(F, N, p)
            if len(set(z[1:])) > 1: nonrat += 1
            okF = ratio_ok(F, N, p)[0]; okK = ratio_ok(Ki, N, p)[0]
            if not (okF and okK): fails += 1; allok = False
        ring = 'Z[i]' if gauss else 'Z'
        print(f"(b) {ring}[C_{N}] (p={p}, e={e}): {tested} pairs F*K = p^e*b tested, v(F(1)) values {sorted(vals)}, "
              f"{nonrat} with irrational F(rho); ratio-test failures: {fails}")
    return allok

# ---------- (c) sharpness ----------
def part_c():
    out = True
    for p in (3, 5, 7):
        F = [(p - 2, 0)] + [(-2, 0)]*(p-1)   # p - 2J in Z[C_p]
        P = conv(F, star(F, p), p)
        assert P[0] == (p*p, 0) and all(t == Z0 for t in P[1:])
        ok, z, one = ratio_ok(F, p, p)
        val = [x for x in z]
        print(f"(c) Z[C_{p}]: F = p - 2J, F F* = p^2 (exponent 2 > 1): F(1) = {one[0]}, "
              f"F(rho) = {val[0][0]-val[1][0]} (const), ratio test passes: {ok}  [expected False]")
        out &= (not ok)
    for (p, e) in ((3, 2), (5, 2)):
        N = p**e
        F = [(N - 2, 0)] + [(-2, 0)]*(N-1)   # p^e - 2J in Z[C_{p^e}], F F* = p^{2e}
        P = conv(F, star(F, N), N)
        assert P[0] == (N*N, 0) and all(t == Z0 for t in P[1:])
        ok = ratio_ok(F, N, p)[0]
        print(f"(c) Z[C_{N}]: F = p^e - 2J, F F* = p^(2e) (scalar exponent 2e > e): ratio test passes: {ok}  [expected False]")
        out &= (not ok)
    # C_3 x C_3, H = <(1,0)>, character phi(a,b) = rho^a (nontrivial on H)
    J = {(a, b): 0 for a in range(3) for b in range(3)}
    F = {g: (-2 if g[1] == 0 else 0) for g in J}; F[(0, 0)] += 3   # 3 - 2 J_H
    # F F^* = 9 check
    prod = {g: 0 for g in J}
    for g1, x in F.items():
        for g2, y in F.items():
            h = ((g1[0]-g2[0]) % 3, (g1[1]-g2[1]) % 3); prod[h] += x*y
    assert prod[(0, 0)] == 9 and all(prod[g] == 0 for g in J if g != (0, 0))
    F1 = sum(F.values())
    Frho = [0, 0, 0]
    for (a, b), x in F.items(): Frho[a % 3] += x
    assert Frho[1] == Frho[2]          # so F(phi) = Frho[0] - Frho[1] is rational
    Fphi = Frho[0] - Frho[1]
    print(f"(c) Z[C_3 x C_3]: F = 3 - 2J_H, F F* = 9: F(1) = {F1}, F(phi) = {Fphi} (phi nontrivial on H): "
          f"ratio = {Fraction(Fphi, F1)}  [non-cyclic 3-group: the lemma's conclusion fails]")
    return out and Fphi == -F1

if __name__ == '__main__':
    samples = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    a = part_a(); b = part_b(samples); c = part_c()
    print("SUMMARY comparison lemma:", "consistent" if (a and b and c) else "INCONSISTENT")

