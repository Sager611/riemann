"""Independent check of UTS Proposition 7.1 in small degrees: for every n <= NMAX list the
self-conjugate lam |- n with g(lam,lam,nu) > 0 for all nu |- n (exact integer characters, own
Murnaghan-Nakayama code mn.py). Expected: none exactly for n in {2,4,9}.
usage: python3 uts_finite.py NMAX"""
import sys
from math import factorial
import mn

def universal_selfconj(n):
    P = list(mn.partitions(n))
    sc = [l for l in P if l == mn.conj(l)]
    table = {}
    for mu in P:
        mn._mn.cache_clear()
        table[mu] = [mn.chi(nu, mu) for nu in P]
    good = []
    for lam in sc:
        li = P.index(lam)
        acc = [0]*len(P)
        for mu in P:
            row = table[mu]
            w = mn.class_size(mu) * row[li]**2
            if w == 0: continue
            for k in range(len(P)):
                acc[k] += w * row[k]
        nf = factorial(n)
        g = []
        for a in acc:
            q, r = divmod(a, nf); assert r == 0 and q >= 0
            g.append(q)
        if min(g) > 0: good.append(lam)
    return sc, good

if __name__ == "__main__":
    NMAX = int(sys.argv[1])
    bad = []
    for n in range(1, NMAX+1):
        sc, good = universal_selfconj(n)
        if not good: bad.append(n)
        print(f"n={n}: self-conjugate {len(sc)}, universal {len(good)}: {good[:3]}{' ...' if len(good)>3 else ''}", flush=True)
    print("degrees with no universal self-conjugate square:", bad)
