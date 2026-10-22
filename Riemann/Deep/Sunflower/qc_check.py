
import itertools, sys, time
from fractions import Fraction as Fr
from math import prod

def is_mono(a): return all(x == a[0] for x in a)
def is_rainbow(a): return len(set(a)) == len(a)

def build_qc(k, D):
    """Return dicts: qc values, list of terms Term_h (h=1..k-1), following l.421-532 literally."""
    A = range(1, D + 1)
    def e(a0, t): return 1 if t == a0 else 0
    def P(h, pre, t):
        # P_h(a_0..a_{h-1}; t): prod_{i=1}^{h-1} (t-a_i)/(a_0-a_i) if prefix rainbow else 0
        if not is_rainbow(pre): return Fr(0)
        return prod((Fr(t - pre[i], pre[0] - pre[i]) for i in range(1, h)), start=Fr(1))
    def R(h, a):  # R_h(a_0..a_h)
        r = Fr(1)
        for j in range(1, h + 1):
            r *= (e(a[0], a[j]) - P(j, a[:j], a[j]))
        return r
    def Delta(h, a):  # prod_{i=h}^{k-1} e(a_0, a_i)
        return prod((e(a[0], a[i]) for i in range(h, k)), start=1)
    q, terms = {}, {h: {} for h in range(1, k)}
    for a in itertools.product(A, repeat=k):
        q[a] = Delta(1, a) - R(k - 1, a)
        for h in range(1, k):
            terms[h][a] = R(h - 1, a) * P(h, a[:h], a[h]) * Delta(h + 1, a)
    return q, terms

def inv_vandermonde(D):
    # V[x][d] = x^d, x in 1..D, d in 0..D-1 ; return V^{-1} as Fractions (coef = Vinv @ values)
    n = D
    M = [[Fr(x) ** d for d in range(n)] + [Fr(int(i == r)) for i in range(n)] for r, x in enumerate(range(1, D + 1))]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        M[c] = [v / pv for v in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [vr - f * vc for vr, vc in zip(M[r], M[c])]
    Vinv = [row[n:] for row in M]  # this is inverse of V (rows indexed by x)... careful
    # We solved V * Vinv = I with V rows indexed by x, cols by d: Vinv rows indexed by d, cols by x
    return Vinv

def reduced_poly(f, k, D):
    """f: dict over [D]^k -> Fr. Returns dict exponent-tuple -> nonzero coefficient."""
    Vinv = inv_vandermonde(D)
    # tensor as nested array indexed by (x_0-1,...,x_{k-1}-1)
    import copy
    T = {tuple(x - 1 for x in a): v for a, v in f.items()}
    for axis in range(k):
        T2 = {}
        for idx in itertools.product(range(D), repeat=k):
            s = Fr(0)
            for xi in range(D):
                j = idx[:axis] + (xi,) + idx[axis + 1:]
                v = T[j]
                if v: s += Vinv[idx[axis]][xi] * v
            T2[idx] = s
        T = T2
    return {d: c for d, c in T.items() if c != 0}

def check(k, D, verbose=True):
    t0 = time.time()
    q, terms = build_qc(k, D)
    ok_vals = all((q[a] == 1) if is_mono(a) else ((q[a] != 0) if is_rainbow(a) else (q[a] == 0)) for a in q)
    ok_sum = all(q[a] == sum(terms[h][a] for h in terms) for a in q)
    poly = reduced_poly(q, k, D)
    strong = all(any(d[i] <= i - 1 for i in range(1, k)) for d in poly)
    weak = all(any(d[i] <= k - 2 for i in range(1, k)) for d in poly)
    term_ok = True
    for h in terms:
        ph = reduced_poly(terms[h], k, D)
        if not all(d[h] <= h - 1 for d in ph): term_ok = False
    worst = max(min(d[i] for i in range(1, k)) for d in poly)
    # negative control: naive tester (value 1 on rainbow)
    naive = {a: (1 if (is_mono(a) or is_rainbow(a)) else 0) for a in q}
    pn = reduced_poly(naive, k, D)
    worst_naive = max(min(d[i] for i in range(1, k)) for d in pn)
    print(f"k={k} D={D}: values(mono=1,rainbow!=0,mixed=0)={ok_vals} identity q=sum(terms)={ok_sum} "
          f"terms deg_{{a_h}}<=h-1={term_ok} monomials={len(poly)} strong(d_i<=i-1 some i)={strong} weak(d_i<=k-2 some i)={weak} "
          f"max_mon min_(i>=1) d_i={worst} | naive tester: monomials={len(pn)} max min_(i>=1) d_i={worst_naive} (D-1={D-1}) [{time.time()-t0:.1f}s]")
    return ok_vals and ok_sum and strong and weak and term_ok

if __name__ == "__main__":
    # hand-derived k=3 formula check: q = [a2=a0] - [a1!=a0](a2-a1)/(a0-a1)
    q3, _ = build_qc(3, 6)
    hand = all(q3[a] == (int(a[2] == a[0]) - (Fr(a[2] - a[1], a[0] - a[1]) if a[1] != a[0] else 0)) for a in q3)
    print("k=3 D=6 matches hand formula q=[a2=a0]-[a1!=a0](a2-a1)/(a0-a1):", hand)
    allok = hand
    for (k, D) in [(3, 3), (3, 5), (3, 7), (4, 4), (4, 5), (4, 6), (4, 7), (4, 9), (5, 5), (5, 6), (5, 7)]:
        allok &= check(k, D)
    print("ALL OK:", allok)

