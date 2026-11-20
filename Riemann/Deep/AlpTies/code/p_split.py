#!/usr/bin/env python3
# Probe C: which rational primes T can carry a prime set ST obeying the Split axioms of Built.Split
# (KerIsTerm.lean:23-37: T in P, T not in P^2, and nonconj: Q != conj P for all P, Q in ST; Q = P gives P != conj P).
# CM fields E = Q(zeta_n), n = 4, 3, 12 (complex conjugation = zeta -> zeta^-1). Two independent methods:
# (1) Galois: P = conj P iff -1 lies in the decomposition group <T mod n>; (2) factor Phi_n mod T over GF(T) and test whether
# each irreducible factor equals its reciprocal (conjugation maps the prime (T, f(zeta)) to (T, f*(zeta))).
from sympy import symbols, cyclotomic_poly, Poly, factor_list, primerange, GF
x = symbols('x')
def galois_selfconj(n, T):
    g, k = T % n, 1
    grp = set()
    v = 1
    while True:
        v = (v * g) % n
        grp.add(v)
        if v == 1: break
    return (n - 1) in grp
def factor_selfconj(n, T):
    f = Poly(cyclotomic_poly(n, x), x, modulus=T)
    facs = f.factor_list()[1]
    res = []
    for q, e in facs:
        d = q.degree()
        coeffs = q.all_coeffs()
        rev = Poly(list(reversed(coeffs)), x, modulus=T)
        lc = rev.LC()
        rev = rev.monic() if lc % T else rev
        res.append(rev == q.monic())
    return all(res), any(res), len(facs)
out = []
mism = 0
for n in (4, 3, 12):
    for T in primerange(2, 200):
        if n % T == 0: continue  # ramified: T in P^2 possible; excluded by T_nsq for P over T | n
        a = galois_selfconj(n, T)
        b_all, b_any, nf = factor_selfconj(n, T)
        if a != b_all or b_all != b_any: mism += 1
        out.append((n, T, nf, a))
for n in (4, 3, 12):
    split = [T for (m, T, nf, a) in out if m == n and not a]
    inert = [T for (m, T, nf, a) in out if m == n and a]
    print(f"E=Q(zeta_{n}): ST can be nonempty (P != conj P) for T = {split[:12]}...; impossible (P = conj P) for T = {inert[:12]}...")
    print(f"   residues mod {n}: allowed {sorted(set(T % n for T in split))}, excluded {sorted(set(T % n for T in inert))}")
print("method mismatches:", mism, "of", len(out))
