#!/usr/bin/env python3
"""Exact check of a cyclic non-separable (NS) family of positive homothetic triangles.

Standard library only, integer arithmetic only, about 0.2 s and 25 MB for config.txt.
Usage: python3 verify.py [config.txt]

Model (arXiv 2609.06558, Sections 2 and 5).  Work in the simplex
Sigma = {z in R^3 : z1+z2+z3 = 1, z >= 0}.  For x >= 0 with |x|_1 < 1 put
T(x) = {z in Sigma : z >= x}, a positive homothet of Sigma with scale t(x) = 1 - |x|_1.
The file lists p^(k) = (X_k, Y_k, Z_k)/DEN for k = 0..m; the family is
P_{k,i} = T(rho^i p^(k)), k = 0..m, i in Z/3, with rho(x1,x2,x3) = (x3,x1,x2);
C_i = P_{0,i} are the three supports.

Checks:
 1. all coordinates >= 0 and every scale > 0; the 3(m+1) members are distinct;
 2. each coordinate has minimum 0, so the smallest homothet of Sigma containing the
    union is Sigma itself (ratio 1 - sum_r min_i x_r = 1, paper Lemma 2.1); hence the
    covering factor is lambda = 1 / (sum of scales);
 3. the 3(2m+1) segment relations of paper Lemma 5.1:
      {P_{k-1,i}, C_{i-1}} -> P_{k,i},  {P_{k,i}, P_{m,i-1}} -> P_{k-1,i}  (1 <= k <= m),
      {P_{m,i}, C_{i-1}} -> P_{m,i+1},
    where {A,B} -> K means: K meets conv(A u B).  As theta*T(x) + (1-theta)*T(y) =
    T(theta x + (1-theta) y) lies in conv(T(x) u T(y)), it suffices to find theta in [0,1]
    with sum_r max(theta x_r + (1-theta) y_r, z_r) <= 1 (paper Lemma 2.3: then a point
    of Sigma is coordinatewise >= both vectors).  We test theta = 0, 1 and the breakpoints
    (z_r - y_r)/(x_r - y_r), exactly;
 4. non-separability by the colouring criterion (paper Lemma 2.2): a line that avoids the
    union and has members on both sides colours the members by side, and {A,B} -> K forces
    colour(K) = colour(A) whenever colour(A) = colour(B).  We run this forcing from every
    colouring of the three supports (up to swapping the colours): "all equal" must force
    every member to that colour, and each 2-1 split must force a contradiction.  So every
    compatible colouring is constant and the family is NS.
"""
import os
import sys
from fractions import Fraction

PAPER = Fraction(250000000, 235141779)   # arXiv 2609.06558, Theorem 1.1 (303 members)

def load(path):
    den, rows = None, {}
    with open(path) as f:
        for line in f:
            s = line.split()
            if not s or s[0].startswith('#'):
                continue
            if s[0] == 'DEN':
                den = int(s[1])
                continue
            k, a, b, c = map(int, s)
            rows[k] = (a, b, c)
    m = max(rows)
    assert den and den > 0 and sorted(rows) == list(range(m + 1)), 'bad file'
    return den, [rows[k] for k in range(m + 1)]

def rot(p, i):
    i %= 3
    return p if i == 0 else ((p[2], p[0], p[1]) if i == 1 else (p[1], p[2], p[0]))

def relation_slack(x, y, z, den):
    """Largest den - M(theta) over the tested theta, as a Fraction; >= 0 means the relation holds."""
    cands = {(0, 1), (1, 1)}
    for r in range(3):
        if x[r] != y[r]:
            p, q = z[r] - y[r], x[r] - y[r]
            if q < 0:
                p, q = -p, -q
            if 0 <= p <= q:
                cands.add((p, q))
    best = None
    for p, q in cands:   # theta = p/q
        qm = sum(max(p * x[r] + (q - p) * y[r], q * z[r]) for r in range(3))
        s = Fraction(q * den - qm, q)
        if best is None or s > best:
            best = s
    return best

def main(path):
    den, P = load(path)
    m = len(P) - 1
    idx, vec = {}, []
    for i in range(3):
        for k in range(m + 1):
            idx[(k, i)] = len(vec)
            vec.append(rot(P[k], i))
    n = len(vec)
    # 1. valid members, pairwise distinct
    assert all(c >= 0 for v in vec for c in v), 'negative coordinate'
    assert all(sum(v) < den for v in vec), 'non-positive scale'
    assert len(set(vec)) == n, 'repeated member'
    # 2. the smallest container is Sigma
    assert all(min(v[r] for v in vec) == 0 for r in range(3)), 'container is not Sigma'
    S = Fraction(sum(den - sum(v) for v in vec), den)
    lam = 1 / S
    # 3. segment relations
    rels = []
    for i in range(3):
        for k in range(1, m + 1):
            rels.append(((k - 1, i), (0, (i - 1) % 3), (k, i)))
            rels.append(((k, i), (m, (i - 1) % 3), (k - 1, i)))
        rels.append(((m, i), (0, (i - 1) % 3), (m, (i + 1) % 3)))
    rels = [tuple(idx[t] for t in r) for r in rels]
    worst = None
    for a, b, c in rels:
        s = relation_slack(vec[a], vec[b], vec[c], den)
        assert s >= 0, 'relation fails'
        worst = s if worst is None or s < worst else worst
    # 4. NS by forced colouring
    by_src = {}
    for a, b, c in rels:
        by_src.setdefault(a, []).append((b, c))
        by_src.setdefault(b, []).append((a, c))
    sup = [idx[(0, i)] for i in range(3)]

    def propagate(col):
        col = dict(col)
        stack = list(col)
        while stack:
            a = stack.pop()
            for b, c in by_src.get(a, ()):
                if b in col and col[b] == col[a]:
                    if c not in col:
                        col[c] = col[a]
                        stack.append(c)
                    elif col[c] != col[a]:
                        return None   # contradiction
        return col

    col = propagate({s: 0 for s in sup})
    assert col is not None and len(col) == n and set(col.values()) == {0}, 'all-equal case not forced'
    for odd in range(3):
        assert propagate({s: int(j == odd) for j, s in enumerate(sup)}) is None, '2-1 case not refuted'
    print('members n = %d (m = %d bridges per chain), relations checked = %d' % (n, m, len(rels)))
    print('min relation slack = %.3e' % float(worst / den))
    print('total scale S = %s' % S)
    print('covering factor lambda = 1/S = %s = %.12f' % (lam, float(lam)))
    print('previous record 250000000/235141779 = %.12f; lambda > previous: %s' % (float(PAPER), lam > PAPER))
    print('NS: proved (all-equal forced constant; all three 2-1 splits contradict)')

if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config.txt'))

