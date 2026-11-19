#!/usr/bin/env python3
"""falsify.py (unitaudit) -- attempts to break UnitCover Lemma C / sec 3.4 on small exact cases.
A: sec 3.4 sums  sum_{3<=q<=Q} phi(q)^m prod|G_q[d_i]|  vs the stated bounds, for every (m, d-tuple) of the k <= 4 closed
   sets, Q <= 4000 (|G_q[d]| by the CRT formula, itself checked against brute force for q <= 400, d <= 16).
B: tightest law for Lemma D: nu = (1-delta) uniform + delta point mass at 1 (Re nuhat(psi) = delta exactly).  Exact
   E[C_r^k] = sum_a (pi(a)+delta(1-pi(a)))^r vs the closed-set bound (3.1); exact Pr(S_N != G) by a DP over subsets
   vs the Lemma C bound min over alpha of E[C_r^k]/t^k + alpha H beta^s, for small groups (Z/qZ)^x."""
import itertools, math, json, os, sys, importlib.util
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("ic", os.path.join(HERE, "indep_classes.py"))
# reuse only the lattice code (not the module's main): exec the definitions up to classes_own
src = open(os.path.join(HERE, "indep_classes.py")).read().split("HERE = os.path")[0]
ns = {}; exec(src, ns)
def fac(n):
    f = {}; p = 2
    while p * p <= n:
        while n % p == 0: f[p] = f.get(p, 0) + 1; n //= p
        p += 1
    if n > 1: f[n] = f.get(n, 0) + 1
    return f
def Gd(q, d):
    r = 1
    for p, a in fac(q).items():
        if p == 2: r *= 1 if a == 1 else (math.gcd(2, d) if a == 2 else math.gcd(2, d) * math.gcd(2 ** (a - 2), d))
        else: r *= math.gcd(d, (p - 1) * p ** (a - 1))
    return r
for q in range(3, 401):
    for d in range(1, 17):
        assert Gd(q, d) == sum(1 for x in range(1, q) if math.gcd(x, q) == 1 and pow(x, d, q) == 1), (q, d)
print("A: CRT formula for |G_q[d]| = brute force, q <= 400, d <= 16: ok", flush=True)
phi = [0] * 4001
for q in range(1, 4001): phi[q] = sum(1 for x in range(1, q + 1) if math.gcd(x, q) == 1) if q < 50 else round(q * math.prod(1 - 1 / p for p in fac(q)))
types = set()
for k in (1, 2, 3, 4):
    for B in ns["lattices"](k):
        rho = len(B); ds = tuple(d for d in ns["smith"](list(B.values())) if d > 1) if rho else ()
        types.add((k - rho, ds))
worst = (0, None)
for (m, ds) in sorted(types):
    n = len(ds); D = math.prod(ds) if ds else 1
    S = 0
    for Q in range(3, 4001):
        S += phi[Q] ** m * math.prod(Gd(Q, d) for d in ds)
        if Q in (3, 4, 5, 6, 8, 10, 12, 30, 60, 100, 210, 500, 1000, 2310 // 2, 2000, 3000, 4000):
            bnd = (Q if m == 0 else Fr(Q ** (m + 1), m + 1)) if n == 0 else 2 ** n * Q ** (m + 1) * (1 + math.log(Q)) ** (D - 1)
            rt = float(S / bnd)
            if rt > worst[0]: worst = (rt, (m, ds, Q))
            assert rt <= 1, (m, ds, Q, rt)
print("A: %d (k - rho, d-tuple) types; sec 3.4 bound holds at all checkpoints Q <= 4000; max ratio %.4f at %s" % (len(types), worst[0], worst[1]), flush=True)

# B
def chars(q):
    """characters of (Z/qZ)^x as dicts via a generator decomposition (brute force: all homs to roots of unity)"""
    U = [x for x in range(1, q) if math.gcd(x, q) == 1]; H = len(U)
    # find cyclic decomposition by brute force search of a basis
    def order(x):
        o, y = 1, x % q
        while y != 1: y = y * x % q; o += 1
        return o
    gens = []; span = {1}
    while len(span) < H:
        best = max((x for x in U if x not in span), key=lambda x: len({s * pow(x, e, q) % q for s in span for e in range(order(x))}))
        o = order(best); new = {s * pow(best, e, q) % q for s in span for e in range(o)}
        # effective order relative to span
        eo = next(e for e in range(1, o + 1) if pow(best, e, q) in span)
        gens.append((best, eo)); span = new
    # represent every unit as exponent vector (search)
    rep = {}
    for ex in itertools.product(*[range(eo) for (_, eo) in gens]):
        x = 1
        for (g, _), e in zip(gens, ex): x = x * pow(g, e, q) % q
        rep.setdefault(x, ex)
    assert len(rep) == H
    return U, gens, rep
def exactB(q, k, delta, r, s):
    U, gens, rep = chars(q); H = len(U)
    ords = [eo for (_, eo) in gens]
    # characters = exponent vectors c; chi(x) = exp(2 pi i sum c_j e_j / o_j); relation test exact via fractions
    chis = list(itertools.product(*[range(o) for o in ords]))
    nonp = [c for c in chis if any(c)]
    cube = list(itertools.product((-1, 0, 1), repeat=k))
    total = Fr(0)
    for a in itertools.product(nonp, repeat=k):
        R = [e for e in cube if all(sum(ei * c[j] for ei, c in zip(e, a)) % ords[j] == 0 for j in range(len(ords)))]
        p = sum(Fr(1, 2 ** (k + sum(1 for x in e if x))) for e in R)
        total += (p + delta * (1 - p)) ** r
    # closed-set bound (3.1) for this group
    bound = Fr(0)
    for B in ns["lattices"](k):
        rho = len(B); ds = [d for d in ns["smith"](list(B.values())) if d > 1] if rho else []
        Rs = [e for e in cube if ns["member"](B, e, k)]
        p = sum(Fr(1, 2 ** (k + sum(1 for x in e if x))) for e in Rs)
        bound += H ** (k - rho) * math.prod(Gd(q, d) for d in ds) * (p + delta * (1 - p)) ** r
    # exact Pr(S_N != G) by DP over subsets (bitmasks over U), law (1-delta) uniform + delta at 1
    idx = {x: i for i, x in enumerate(U)}; full = (1 << H) - 1
    mult = [[idx[x * y % q] for y in U] for x in U]
    dist = {1 << idx[1]: Fr(1)}
    for _ in range(r + s):
        nd = {}
        for S, pr in dist.items():
            elems = [i for i in range(H) if S >> i & 1]
            for g in range(H):
                w = (1 - delta) / H + (delta if U[g] == 1 else 0)
                S2 = S
                for i in elems: S2 |= 1 << mult[g][i]
                nd[S2] = nd.get(S2, 0) + pr * w
        dist = nd
    fail = 1 - dist.get(full, Fr(0))
    best = None
    for al in [Fr(j, 100) for j in range(1, 100)]:
        t = al / (1 - al); beta = delta + (1 - delta) * al
        val = total / t ** k + al * H * beta ** s
        best = val if best is None or val < best else best
    return H, float(total), float(bound), float(fail), float(best)
for (q, k, delta, r, s) in [(5, 2, Fr(1, 10), 10, 4), (7, 2, Fr(1, 20), 12, 4), (8, 4, Fr(1, 10), 14, 4), (9, 3, Fr(1, 10), 14, 4),
                            (12, 4, Fr(1, 20), 14, 5), (15, 2, Fr(1, 10), 16, 5), (16, 3, Fr(1, 10), 14, 5), (5, 4, Fr(0), 12, 3)]:
    H, tot, bnd, fail, best = exactB(q, k, delta, r, s)
    print("B: q=%d (H=%d) k=%d delta=%s r=%d s=%d: E[C_r^k] = %.4g <= (3.1) %.4g: %s; exact Pr(S_N != G) = %.3g <= Lemma C bound %.3g: %s"
          % (q, H, k, delta, r, s, tot, bnd, tot <= bnd * (1 + 1e-12), fail, best, fail <= best), flush=True)
