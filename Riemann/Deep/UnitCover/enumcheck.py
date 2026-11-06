#!/usr/bin/env python3
"""enumcheck.py -- independent checks of the closed-relation-set enumeration of cert.py (DERIVATION sec 3.3).
(1) k = 2, 3: brute force over ALL subsets S of {-1,0,1}^k closed under negation and containing 0: S is closed iff
    span_Z(S) cap cube = S and span_Z(S) has no unit vector; the closed ones must be exactly cert.closed_sets(k).
    Lattice membership here by exact rational elimination (independent of cert.hnf).
(2) k = 4: relation sets of 30000 random and 3000 structured tuples in 12 abelian groups with elements of order up to 840
    (orders 5, 7, 9, 10 included) must all be enumerated sets, with the same pi.
Run: python3 enumcheck.py  (imports cert.py; pure integers/fractions; about 3 min)."""
import itertools, random, sys, os, time
from fractions import Fraction as Fr
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cert

def rank_and_member(gens, k):
    """independent lattice test: returns function v -> (v in Z-span of gens)"""
    # rational row echelon of gens; v in Q-span? then integrality via Smith-free method: enumerate? use HNF from scratch
    # (Euclid on columns with Fraction-free integer ops, written independently of cert.hnf)
    rows = [list(g) for g in gens if any(g)]
    basis = []
    for c in range(k):
        rows = [r for r in rows if any(r)]
        piv = [r for r in rows if r[c] != 0]
        while len(piv) > 1:
            piv.sort(key=lambda r: abs(r[c]))
            p0 = piv[0]
            for r in piv[1:]:
                m = r[c] // p0[c]
                for j in range(k): r[j] -= m * p0[j]
            rows = [r for r in rows if any(r)]
            piv = [r for r in rows if r[c] != 0]
        if piv:
            basis.append(piv[0]); rows = [r for r in rows if r is not piv[0]]
    def member(v):
        v = list(v)
        for bvec in basis:
            c = next(j for j in range(k) if bvec[j] != 0)
            if v[c] % bvec[c]: return False
            m = v[c] // bvec[c]
            v = [x - m * y for x, y in zip(v, bvec)]
        return not any(v)
    return member

def brute_closed(k):
    cube = list(itertools.product((-1, 0, 1), repeat=k))
    reps = [e for e in cube if any(e) and next(x for x in e if x) > 0]
    units = [tuple(1 if i == j else 0 for i in range(k)) for j in range(k)]
    found = set()
    for mask in range(1 << len(reps)):
        S = [tuple(0 for _ in range(k))]
        for i, e in enumerate(reps):
            if mask >> i & 1:
                S.append(e); S.append(tuple(-x for x in e))
        mem = rank_and_member(S, k)
        if any(mem(u) for u in units): continue
        if all((e in S) == mem(e) for e in cube):
            found.add(frozenset(S))
    return found

def group_relset(ns, a, k, cube):
    return frozenset(e for e in cube if all(sum(ei * ai[j] for ei, ai in zip(e, a)) % ns[j] == 0 for j in range(len(ns))))

if __name__ == "__main__":
    t0 = time.time(); bad = 0
    for k in (2, 3):
        enum = {R for (pi, rho, ed, R) in cert.closed_sets(k)}
        brute = brute_closed(k)
        ok = enum == brute
        bad += not ok
        print("enumcheck k=%d: brute force over all symmetric subsets finds %d closed sets, enumeration %d: %s"
              % (k, len(brute), len(enum), "identical" if ok else "MISMATCH"), flush=True)
    k = 4
    cl = cert.closed_sets(k)
    index = {R: pi for (pi, rho, ed, R) in cl}
    cube = list(itertools.product((-1, 0, 1), repeat=k))
    groups = [[840], [2, 420], [6, 60], [2, 2, 210], [10, 10], [9, 9], [7, 7], [2, 4, 8], [3, 3, 3], [5, 5, 2], [2, 2, 2, 2], [12, 12]]
    rnd = random.Random(20261106); seen = set(); n = 0
    for ns in groups:
        elems_nz = lambda: tuple(rnd.randrange(n_) for n_ in ns)
        for _ in range(2500):
            a = []
            while len(a) < k:
                x = elems_nz()
                if any(x): a.append(x)
            R = group_relset(ns, a, k, cube); n += 1
            if R not in index: print("MISSING (random)", ns, a); bad += 1
            seen.add(R)
        # structured: a_i = c_i * x for small integer vectors c, and a_i in span of two random elements
        for _ in range(250):
            x = elems_nz(); y = elems_nz()
            cs = [rnd.randint(-4, 4) for _ in range(k)]; ds = [rnd.randint(-3, 3) for _ in range(k)]
            a = [tuple((ci * xi + di * yi) % n_ for xi, yi, n_ in zip(x, y, ns)) for ci, di in zip(cs, ds)]
            if any(not any(ai) for ai in a): continue
            R = group_relset(ns, a, k, cube); n += 1
            if R not in index: print("MISSING (structured)", ns, a); bad += 1
            seen.add(R)
    print("enumcheck k=4: %d tuples in %d groups; %d distinct relation sets met, all among the %d enumerated: %s"
          % (n, len(groups), len(seen), len(cl), "ok" if bad == 0 else "FAIL"))
    print("enumcheck: %s in %.0f s" % ("ALL OK" if bad == 0 else "%d FAILURES" % bad, time.time() - t0))
