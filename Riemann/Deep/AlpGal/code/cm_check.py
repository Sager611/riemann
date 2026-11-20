# cm_check.py: instantiate Row0b's unramified local-global condition (AwayB.tor) for E = Q(i), a = 1,
# pi = psi, the Hecke character of the CM elliptic curve y^2 = x^3 - x, and test it against Frobenius traces.
# Archive normalizations (read, not assumed): deg 1 = arithmetic Frobenius (IsFrob, Statement.lean:28; mem_Weil_of_isFrob,
# ExistsNotCongPow.lean:159); local Artin map a with ord(a s) = -deg s (IsLocalArtin.ord_eq, HasExponentZeroIff.lean:149);
# rc(character chi) = chi o a (LawsT.one, IsRecTOfIsRec.lean:58). Hence for unramified w: iota r(Frob_arith,w) = psi_w(pi_w)^(-1),
# i.e. r(Frob_geom,w) = psi(w) = primary generator of w. Test: H^1(E) restricted to G_Q(i) = r + r^c, so
#   split p = N(w):  a_p = alpha + conj(alpha),  alpha = a+bi primary (alpha = 1 mod (2+2i)), alpha*conj(alpha) = p;
#   inert p (w=(p), N w = p^2): trace Frob_{p^2} = a_p^2 - 2p must equal psi((p)) + psi^c((p)) = 2*(-p) (-p is primary).
# Exact integers; every quantity computed twice by independent methods.
import sys
PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 500
def primes(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]
def ap_count(p):            # method A: #E(F_p) by counting square roots via a table
    nsq = [0] * p
    for y in range(p): nsq[y * y % p] += 1
    n = 1 + sum(nsq[(x * x * x - x) % p] for x in range(p))
    return p + 1 - n
def ap_euler(p):            # method B: -sum of Legendre symbols by Euler's criterion
    s = 0
    for x in range(p):
        t = (x * x * x - x) % p
        if t: s += 1 if pow(t, (p - 1) // 2, p) == 1 else -1
    return -s
def gmod(z, m):             # is Gaussian integer z divisible by m (both as (re, im))?
    a, b = z; c, d = m; n = c * c + d * d
    re, im = a * c + b * d, b * c - a * d     # z * conj(m)
    return re % n == 0 and im % n == 0
def primary(a, b):          # the unique associate (among +-a+-bi, +-b+-ai) that is = 1 mod (2+2i); None if none
    cands = [(a, b), (-a, -b), (-b, a), (b, -a), (a, -b), (-a, b), (b, a), (-b, -a)]
    good = [z for z in cands if gmod((z[0] - 1, z[1]), (2, 2))]
    return good
def alpha_brute(p):         # method 1: brute force a^2 + b^2 = p
    for a in range(1, int(p ** 0.5) + 1):
        b2 = p - a * a; b = int(round(b2 ** 0.5))
        if b * b == b2: return a, b
def alpha_cornacchia(p):    # method 2: Hermite-Serret: x^2 = -1 mod p, Euclid on (p, x) until remainder < sqrt(p)
    for g in range(2, p):
        if pow(g, (p - 1) // 2, p) == p - 1: break
    x = pow(g, (p - 1) // 4, p)
    r0, r1 = p, x
    while r1 * r1 > p: r0, r1 = r1, r0 % r1
    a = r1; b = int(round((p - a * a) ** 0.5))
    assert a * a + b * b == p
    return a, b
def count_Fp2(p):           # #E(F_{p^2}), p = 3 mod 4, F_{p^2} = F_p[i]
    sq = {}
    for u in range(p):
        for v in range(p):
            k = ((u * u - v * v) % p, (2 * u * v) % p); sq[k] = sq.get(k, 0) + 1
    n = 1
    for u in range(p):
        for v in range(p):
            # x = u + v i ; x^3 - x
            x2 = ((u * u - v * v) % p, (2 * u * v) % p)
            x3 = ((x2[0] * u - x2[1] * v) % p, (x2[0] * v + x2[1] * u) % p)
            t = ((x3[0] - u) % p, (x3[1] - v) % p)
            n += sq.get(t, 0)
    return n
def count_Fp2_split(p, i0):  # #E(F_{p^2}) for p = 1 mod 4 via F_p[s]/(s^2 - g), g a non-residue
    for g in range(2, p):
        if pow(g, (p - 1) // 2, p) == p - 1: break
    mul = lambda z, w: ((z[0] * w[0] + g * z[1] * w[1]) % p, (z[0] * w[1] + z[1] * w[0]) % p)
    sq = {}
    for u in range(p):
        for v in range(p):
            k = mul((u, v), (u, v)); sq[k] = sq.get(k, 0) + 1
    n = 1
    for u in range(p):
        for v in range(p):
            z = (u, v); z3 = mul(mul(z, z), z)
            n += sq.get(((z3[0] - u) % p, (z3[1] - v) % p), 0)
    return n
bad = []; nsplit = ninert = nsq = 0
print('p  class  a_p(A)  a_p(B)  alpha(brute)  alpha(Cornacchia)  alpha+conj  check')
for p in primes(PMAX):
    if p == 2: continue          # bad reduction (conductor 32)
    A, B = ap_count(p), ap_euler(p)
    if A != B: bad.append((p, 'A!=B', A, B))
    if p % 4 == 3:
        ninert += 1
        if A != 0: bad.append((p, 'inert a_p != 0', A))
        prim = primary(-p, 0)
        if (-p, 0) not in prim: bad.append((p, '-p not primary', prim))
        tr2 = A * A - 2 * p      # trace of Frob_p^2 = Frob_(p) on H^1
        if tr2 != 2 * (-p): bad.append((p, 'inert trace', tr2))
        line = '%d inert %d %d  psi((p))=-p  tr Frob_(p) = %d = 2(-p)' % (p, A, B, tr2)
    else:
        nsplit += 1
        a1, b1 = alpha_brute(p); a2, b2 = alpha_cornacchia(p)
        P1, P2 = primary(a1, b1), primary(a2, b2)
        # the two primary associates are alpha and conj(alpha) (primes over p of the two factors) -> same real part
        re1 = {z[0] for z in P1}; re2 = {z[0] for z in P2}
        ok = len(P1) == 2 and len(P2) == 2 and re1 == re2 and len(re1) == 1 and A == 2 * next(iter(re1))
        if not ok: bad.append((p, 'split', A, P1, P2))
        line = '%d split %d %d  %s  %s  %d  %s' % (p, A, B, P1[0] if P1 else None, P2[0] if P2 else None, 2 * next(iter(re1)) if re1 else 0, 'ok' if ok else 'MISMATCH')
    print(line)
# extra: Frobenius-squared check over F_{p^2} for small p
print('--- F_{p^2} point counts (p <= 60): trace of Frob^2 on H^1 vs alpha^2 + conj(alpha)^2 (split) or 2(-p) (inert)')
for p in primes(60):
    if p == 2: continue
    if p % 4 == 3:
        N2 = count_Fp2(p); t = p * p + 1 - N2; want = -2 * p
    else:
        N2 = count_Fp2_split(p, None); t = p * p + 1 - N2
        a, b = primary(*alpha_brute(p))[0]; want = 2 * (a * a - b * b)   # alpha^2 + conj^2 = 2 Re(alpha^2)
    nsq += 1
    if t != want: bad.append((p, 'Fp2', t, want))
    print(p, 'F_p2 trace', t, 'expected', want, 'ok' if t == want else 'MISMATCH')
print('SUMMARY split=%d inert=%d Fp2=%d mismatches=%d' % (nsplit, ninert, nsq, len(bad)))
for x in bad: print('BAD', x)
