#!/usr/bin/env python3
"""exp48 conductor.py: completeness of 2605.26020's tables E = 2, 4, 8 beyond our scan bound X1.
Run: python3 conductor.py SCANFILE X1      (pure Python, no packages)
Input: SCANFILE = output of scan.c for 3 <= |D| <= X1 (every discriminant D < 0 with exp Cl(D) | 8, with exact E).
Logic (unconditional given Theorem A):
  * D = D0 f^2, Cl(D) ->> Cl(D0) (Cox Cor. 7.17), so E(D) | E forces E(D0) | E: D0 is in the scan list (Theorem A:
    every fundamental D0 with E(D0) | 8 has |D0| <= 430950520 < X1).
  * 2605.26020 Lemma property-of-f (Fan-Pollack 2503.07801 Prop 2.3, algebraic): L(f) | 12 u E with
    L(f) = lcm_{p^k || f} p^(k-1) (p - (D0/p)), u = |O_K^x| (6, 4, 2).
  * Cl(D0 f^2) ->> Cl(D0 f'^2) for f' | f: if some divisor f' has |D0| f'^2 <= X1 and D0 f'^2 is not in the
    scan list with exponent dividing E, then E(D0 f^2) does not divide E.
Every candidate f with |D0| f^2 > X1 must be excluded this way; otherwise it is reported UNRESOLVED.
"""
import sys
from math import gcd
scanfile, X1 = sys.argv[1], int(sys.argv[2])
EXP = {}
for line in open(scanfile):
    D, fund, h, E, n2, n4 = map(int, line.split())
    EXP[D] = (E, fund)
def kron(D, p):
    if p == 2:
        if D % 2 == 0: return 0
        return 1 if D % 8 in (1, 7) else -1
    r = D % p
    if r == 0: return 0
    return 1 if pow(r, (p - 1) // 2, p) == 1 else -1
def primes(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]
PR = primes(2 * 12 * 8 * 6 + 10)
def lcm(a, b): return a // gcd(a, b) * b
total_cands = 0; unresolved = []; big = 0
def small_ok_global(Et, D0, fp):
    D = D0 * fp * fp
    return D in EXP and Et % EXP[D][0] == 0
for Etarget in (2, 4, 8):
    fund = sorted(D for D, (E, f) in EXP.items() if f and Etarget % E == 0)
    found = set(D for D, (E, f) in EXP.items() if E == Etarget)
    nc = 0; nbig = 0
    for D0 in fund:
        u = 6 if D0 == -3 else (4 if D0 == -4 else 2)
        N = 12 * u * Etarget
        # allowed prime powers: p^r with L(p^r) = p^(r-1)(p - chi(p)) dividing N
        allowed = []
        for p in PR:
            ch = kron(D0, p); L1 = p - ch
            if L1 <= 0 or N % L1: continue
            r, Lr, opts = 1, L1, []
            while N % Lr == 0:
                opts.append((r, Lr)); r += 1; Lr *= p
            allowed.append((p, opts))
        memo = {}
        def small_ok(fp):
            if fp not in memo: memo[fp] = (D0 * fp * fp in EXP and Etarget % EXP[D0 * fp * fp][0] == 0)
            return memo[fp]
        def bad(fp):  # fp is a divisor of every candidate in the branch: a certified exclusion if True
            return (-D0) * fp * fp <= X1 and not small_ok(fp)
        # DFS over candidate f = prod p^r, tracking lcm; a branch is cut as soon as one divisor of f is bad
        stack = [(0, 1, 1, ())]
        while stack:
            i, f, L, fac = stack.pop()
            if i == len(allowed):
                nc += 1
                if f == 1: continue
                D = D0 * f * f
                if -D <= X1: continue           # covered by the scan list directly
                nbig += 1
                # remaining divisors d of f with |D0| d^2 <= X1 (size-bounded generation)
                lim = int((X1 // (-D0)) ** 0.5) + 1
                divs = [1]
                for (q, r) in fac:
                    divs = [d * q ** j for d in divs for j in range(r + 1) if d * q ** j <= lim]
                if not any(d > 1 and bad(d) for d in divs): unresolved.append((Etarget, D0, f))
                continue
            q, opts = allowed[i]
            stack.append((i + 1, f, L, fac))
            for (r, Lr) in opts:
                L2 = lcm(L, Lr)
                if N % L2: continue
                f2 = f * q ** r
                if bad(q ** r) or bad(f2): continue   # every extension has the bad divisor q^r or f2
                stack.append((i + 1, f2, L2, fac + ((q, r),)))
    total_cands += nc; big += nbig
    print(f"E={Etarget}: {len(fund)} fundamental D0 with E(D0) | E; {nc} conductor candidates (f=1 included); {nbig} with |D0 f^2| > X1; unresolved so far {len([x for x in unresolved if x[0]==Etarget])}; discriminants with exponent exactly E and |D| <= X1: {len(found)}", flush=True)
# count of all candidates (no cuts) for D0 = -3, E = 8, to compare with ARSX's README ("0.8M candidates")
def count_all(D0, Et):
    u = 6 if D0 == -3 else (4 if D0 == -4 else 2); N = 12 * u * Et; allowed = []
    for p in PR:
        ch = kron(D0, p); L1 = p - ch
        if L1 <= 0 or N % L1: continue
        r, Lr, opts = 1, L1, []
        while N % Lr == 0: opts.append(Lr); r += 1; Lr *= p
        allowed.append(opts)
    cnt = 0; stack = [(0, 1)]
    while stack:
        i, L = stack.pop()
        if i == len(allowed): cnt += 1; continue
        stack.append((i + 1, L))
        for Lr in allowed[i]:
            L2 = lcm(L, Lr)
            if N % L2 == 0: stack.append((i + 1, L2))
    return cnt
print("all conductor candidates without cuts, D0 = -3, E = 8:", count_all(-3, 8), flush=True)
# stage 2: for D = D0 d^2 beyond X1, certified exclusion by exact form arithmetic (Python big integers):
#   a split prime q (q !| D, (D/q) = 1) with 4 q^8 < |D| (EKN Lemma 1, valid for orders), or a prime form f_q with
#   f_q^E' non-principal for E' = Etarget (explicit element whose order does not divide Etarget).
def egcd(a, b):
    x0, y0, x1, y1 = 1, 0, 0, 1
    while b:
        q = a // b; a, b = b, a - q * b; x0, x1 = x1, x0 - q * x1; y0, y1 = y1, y0 - q * y1
    if a < 0: a, x0, y0 = -a, -x0, -y0
    return a, x0, y0
def red(a, b, D):
    def norm(a, b):
        r = b % (2 * a)
        if r > a: r -= 2 * a
        return a, r, (r * r - D) // (4 * a)
    a, b, c = norm(a, b)
    while a > c or (a == c and b < 0):
        a, b, c = norm(c, -b)
    return a, b, c
def comp(f1, f2, D):  # Cohen Alg. 5.4.7
    if f1[0] > f2[0]: f1, f2 = f2, f1
    a1, b1, _ = f1; a2, b2, c2 = f2
    s = (b1 + b2) // 2; n = b2 - s
    if a2 % a1 == 0: y1, d = 0, a1
    else: d, y1, _v = egcd(a2, a1)
    if s % d == 0: y2, x2, d1 = -1, 0, d
    else: d1, x2, y2 = egcd(s, d); y2 = -y2
    v1, v2 = a1 // d1, a2 // d1
    r = (y1 * y2 * n - x2 * c2) % v1
    a3 = v1 * v2; b3 = b2 + 2 * v2 * r
    return red(a3, b3, D)
def sqrtmod(r, q):
    for t in range(q):
        if t * t % q == r % q: return t
    return None
def excl(D, Et, tries=8):
    n = -D; done = 0
    for q in PR + primes(5000)[len(PR):]:
        if done >= tries: break
        ch = kron(D, q)
        if ch != 1: continue
        if 4 * q ** 8 < n: return "L1(q=%d)" % q
        if q == 2: b = 1
        else:
            t = sqrtmod(D % q, q); b = t if (t - D) % 2 == 0 else t + q
        f = red(q, b, D); g = f
        for _ in range({2: 1, 4: 2, 8: 3}[Et]): g = comp(g, g, D)
        if g[0] != 1: return "pow(q=%d)" % q
        done += 1
    return None
# soundness self-test of excl(): it must never exclude a discriminant whose exponent divides Et (scan list)
bad_excl = [(D, Et) for D, (E, f) in EXP.items() for Et in (2, 4, 8) if Et % E == 0 and excl(D, Et)]
print("self-test: excl() on %d known (D, Et) with E(D) | Et: false exclusions = %d" % (sum(1 for D, (E, f) in EXP.items() for Et in (2, 4, 8) if Et % E == 0), len(bad_excl)), flush=True)
assert not bad_excl, bad_excl[:5]
# and it excludes the scan's exponent-8 discriminants for Et = 4 (order-8 element found or Lemma 1), informational
e8 = [D for D, (E, f) in EXP.items() if E == 8]
print("self-test: excl(D, 4) finds an exclusion for %d of the %d exponent-8 discriminants" % (sum(1 for D in e8 if excl(D, 4)), len(e8)), flush=True)
still = []; how = {}
for (Et, D0, f) in unresolved:
    fac = []; m = f
    for q in PR:
        if m % q == 0:
            r = 0
            while m % q == 0: m //= q; r += 1
            fac.append((q, r))
    assert m == 1
    divs = [1]
    for (q, r) in fac: divs = [d * q ** j for d in divs for j in range(r + 1)]
    res = None
    for d in sorted(divs)[1:]:
        if (-D0) * d * d <= X1:
            if not small_ok_global(Et, D0, d): res = "scan(d=%d)" % d; break
            continue
        e = excl(D0 * d * d, Et)
        if e: res = e + " at d=%d" % d; break
    if res is None: still.append((Et, D0, f))
    else: how[res.split("(")[0]] = how.get(res.split("(")[0], 0) + 1
print("stage 2 exclusions by kind:", how, flush=True)
print("UNRESOLVED after stage 2:", still[:20], len(still))
print("CONDUCTOR STEP COMPLETE: no discriminant beyond X1" if not still else "CONDUCTOR STEP INCOMPLETE")
