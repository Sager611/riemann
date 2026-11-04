"""Cyclic-module supports by brute-force isotypic projection (cyc.py), exact integers.
 S1  Saxl Thm 3.1 (cyclic staircase support): W_m = C[S_N] (v_R (x) v_C) contains every S^lam, lam |- N_m, m<=4
 S2  Saxl Prop 5.3 (two-band): U_{m,2} contains exactly the eta |- 2m-1 with <=4 rows, m=3,4,5
 S3  UTS Lemma 4.5 (odd-path support): the odd path tensor of length k contains every partition of k with <=4 rows, k<=9
 S4  UTS Lemma 4.6 (one extra copy): Z_E (x) sgn contains every alpha |- 2b+delta with <=4 rows and
     all rows even (delta=0) / exactly one odd row (delta=1), E=(b+delta,b), 2b+delta<=9
usage: python3 saxl_cyclic.py [S1 S2 S3 S4]"""
import sys, itertools, time
import mn, cyc
from saxl_band import triangle, rows_of, cols_of, alt_blocks, band_cells

def layer_terms(blocks, positions):
    idx = {c: k for k, c in enumerate(positions)}
    out = []
    for asg, c in alt_blocks(blocks):
        w = [0]*len(positions)
        for cell, a in asg.items(): w[idx[cell]] = a-1
        out.append((tuple(w), c))
    return out

def report(name, norms, expected):
    got = {g for g, v in norms.items() if v != 0}
    neg = [g for g, v in norms.items() if v < 0]
    assert not neg, ("negative norm", name, neg)
    missing = sorted(expected - got)
    extra = sorted(got - expected)
    print(f"{name}: support size {len(got)}/{len(norms)}; claimed {len(expected)}; missing={missing}; extra(not claimed)={extra}", flush=True)
    return missing

def S1(mmax=4):
    for m in range(1, mmax+1):
        B = triangle(m)
        L1 = layer_terms(rows_of(B), B); L2 = layer_terms(cols_of(B), B)
        N = len(B); t = time.time()
        norms = cyc.support([L1, L2], N)
        report(f"S1 Saxl W_{m} (N={N}, {time.time()-t:.0f}s)", norms, set(mn.partitions(N)))

def S2():
    for m in (3, 4, 5):
        F = band_cells(m)
        L1 = layer_terms(rows_of(F), F); L2 = layer_terms(cols_of(F), F)
        K = len(F)
        norms = cyc.support([L1, L2], K)
        report(f"S2 Saxl U_({m},2) (K={K})", norms, {e for e in mn.partitions(K) if len(e) <= 4})

def S3():
    for k in (1, 3, 5, 7, 9):
        pos = list(range(1, k+1))
        # layer 1: singleton at 1, pairs (2,3),(4,5),...; layer 2: pairs (1,2),(3,4),..., singleton at k
        b1 = [[1]] + [[i, i+1] for i in range(2, k, 2)]
        b2 = [[i, i+1] for i in range(1, k-1, 2)] + [[k]]
        L1 = layer_terms(b1, pos); L2 = layer_terms(b2, pos)
        norms = cyc.support([L1, L2], k)
        report(f"S3 UTS odd path k={k}", norms, {e for e in mn.partitions(k) if len(e) <= 4})

def S4():
    for b in range(0, 5):
        for delta in (0, 1):
            n = 2*b + delta
            if n == 0 or n > 9: continue
            E = [(1, j) for j in range(1, b+delta+1)] + [(2, j) for j in range(1, b+1)]
            L1 = layer_terms(rows_of(E), E); L2 = layer_terms(cols_of(E), E)
            norms = cyc.support([L1, L2], n)
            # Z_E (x) sgn contains alpha  <=>  Z_E contains alpha^t
            tw = {mn.conj(g): v for g, v in norms.items()}
            def ok(a):
                odd = sum(1 for x in a if x % 2)
                return len(a) <= 4 and (odd == 0 if delta == 0 else odd == 1)
            report(f"S4 UTS one-extra E=({b+delta},{b})", tw, {a for a in mn.partitions(n) if ok(a)})

if __name__ == "__main__":
    which = sys.argv[1:] or ["S1", "S2", "S3", "S4"]
    for w in which:
        globals()[w]()
