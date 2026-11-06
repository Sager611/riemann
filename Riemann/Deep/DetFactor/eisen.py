"""Theorem E check (DERIVATION section 4): small r-th power non-residues a + b*w of F_{p^2} = F_p[w]/(w^2+w+1),
p = 2 mod 3 (and a + b*i of F_p[i]/(i^2+1), p = 3 mod 4), against the bounds
  N(alpha) <= (0.49 (log 3 + r log|D_K| + 2 (r-1) log p))^32,  |D_K| = 3 resp. 4,
valid when every unit of K is an r-th power in F_{p^2} (otherwise a unit is the non-residue).
For each prime p <= PMAX in the class and r in {2, 3, 5, 7} with r | p^2 - 1: least norm of a non-residue with
a, b >= 0 small (search by increasing norm), maximum over p, and slack log10(bound / least norm).
Run: python3 eisen.py [PMAX]
"""
import math, sys
PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 30000
def primes(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]
def mul(x, y, p, c):          # (a + b u)(c + d u), u^2 = -c1 u - c0 with c = (c0, c1)
    a, b = x; cc, d = y; c0, c1 = c
    return ((a * cc - b * d * c0) % p, (a * d + b * cc - b * d * c1) % p)
def pw(x, e, p, c):
    r = (1, 0)
    while e:
        if e & 1: r = mul(r, x, p, c)
        x = mul(x, x, p, c); e >>= 1
    return r
def run(name, cond, c, normf, DK, units):
    print('== %s' % name)
    for r in [2, 3, 5, 7]:
        worst, worstp, minslack, cnt, unitcases = 0, 0, 1e9, 0, 0
        for p in primes(PMAX):
            if p < 5 or not cond(p) or (p * p - 1) % r: continue
            Q = p * p; e = (Q - 1) // r; cnt += 1
            if any(pw(u, e, p, c) != (1, 0) for u in units(p)):
                unitcases += 1; continue
            cands = sorted(((a, b) for a in range(0, 60) for b in range(0, 60) if (a, b) != (0, 0)), key=lambda t: normf(*t))
            best = None
            for a, b in cands:
                if a % p == 0 and b % p == 0: continue
                if pw((a % p, b % p), e, p, c) != (1, 0):
                    best = normf(a, b); break
            assert best is not None
            bound = 32 * math.log10(0.49 * (math.log(3) + r * math.log(DK) + 2 * (r - 1) * math.log(p)))
            minslack = min(minslack, bound - math.log10(best))
            if best > worst: worst, worstp = best, p
        print('r=%d: %5d primes (%d with a unit non-residue); max least norm %d (p=%d); min log10(bound/norm) = %.1f'
              % (r, cnt, unitcases, worst, worstp, minslack))
        assert minslack > 0
# K = Q(sqrt-3): w^2 = -w - 1 -> c = (1, 1); norm a^2 - ab + b^2; units +-1, +-w, +-w^2
run('F_p[w], p = 2 mod 3', lambda p: p % 3 == 2, (1, 1), lambda a, b: a * a - a * b + b * b, 3,
    lambda p: [(p - 1, 0), (0, 1), (0, p - 1), (p - 1, p - 1), (1, 1)])
# K = Q(i): i^2 = -1 -> c = (1, 0); norm a^2 + b^2; units +-1, +-i
run('F_p[i], p = 3 mod 4', lambda p: p % 4 == 3, (1, 0), lambda a, b: a * a + b * b, 4,
    lambda p: [(p - 1, 0), (0, 1), (0, p - 1)])
print('ALL THEOREM E CHECKS PASSED')
