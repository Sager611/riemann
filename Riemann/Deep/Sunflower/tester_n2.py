
import itertools, sys, time, random
from fractions import Fraction as Fr
from math import lcm, comb
sys.path.insert(0, '.')
from qc_check import build_qc, inv_vandermonde, is_mono, is_rainbow

def run(k, D, n=2):
    t0 = time.time()
    q, _ = build_qc(k, D)
    L = 1
    for v in q.values(): L = lcm(L, v.denominator)
    qi = {a: int(v * L) for a, v in q.items()}
    Vinv = inv_vandermonde(D)
    s = 1
    for row in Vinv:
        for v in row: s = lcm(s, v.denominator)
    W = [[int(v * s) for v in row] for row in Vinv]
    nv = k * n  # variables ordered (i, r) -> index i*n + r
    N = D ** nv
    # flat array of T values * L^n, index = sum digit_j * D^(nv-1-j), digit = x-1
    T = [0] * N
    tester_ok = True
    for idx in range(N):
        digs = []
        m = idx
        for _ in range(nv):
            digs.append(m % D); m //= D
        digs.reverse()
        x = [[digs[i * n + r] + 1 for r in range(n)] for i in range(k)]
        val = 1
        for r in range(n):
            val *= qi[tuple(x[i][r] for i in range(k))]
        T[idx] = val
        # tester property: nonzero iff every coordinate mono or rainbow
        cond = all(is_mono([x[i][r] for i in range(k)]) or is_rainbow([x[i][r] for i in range(k)]) for r in range(n))
        if (val != 0) != cond: tester_ok = False
    Tvals = T[:]
    # inverse Vandermonde along each axis
    for axis in range(nv):
        stride = D ** (nv - 1 - axis)
        T2 = [0] * N
        for idx in range(N):
            d = (idx // stride) % D
            base = idx - d * stride
            acc = 0
            Wd = W[d]
            for xi in range(D):
                v = T[base + xi * stride]
                if v: acc += Wd[xi] * v
            T2[idx] = acc
        T = T2
    coef = {}
    for idx in range(N):
        if T[idx]:
            digs = []; m = idx
            for _ in range(nv):
                digs.append(m % D); m //= D
            digs.reverse()
            coef[tuple(digs)] = T[idx]
    # reconstruction check at 300 random points: sum c_d prod x^d == T * s^nv
    scale = s ** nv
    rec_ok = True
    for _ in range(300):
        idx = random.randrange(N)
        m = idx; digs = []
        for _ in range(nv):
            digs.append(m % D); m //= D
        digs.reverse()
        xs = [dg + 1 for dg in digs]
        tot = 0
        for d, c in coef.items():
            p = c
            for j in range(nv): p *= xs[j] ** d[j]
            tot += p
        if tot != Tvals[idx] * scale: rec_ok = False; break
    # monomial structure
    per_coord = all(all(any(d[i * n + r] <= k - 2 for i in range(1, k)) for r in range(n)) for d in coef)
    need = -(-n // (k - 1))  # ceil(n/(k-1))
    pigeon = all(any(sum(d[i * n + r] <= k - 2 for r in range(n)) >= need for i in range(1, k)) for d in coef)
    smax = int((1 - 1 / (k - 1)) * n + 1e-12)
    # x_i-factor lies in P (at most smax coordinates with exponent >= k-1) for the pigeonhole-chosen i
    inP = all(any(sum(d[i * n + r] >= k - 1 for r in range(n)) <= smax for i in range(1, k)) for d in coef)
    dimP = sum(comb(n, t) * (D - k + 1) ** t * (k - 1) ** (n - t) for t in range(smax + 1))
    print(f"k={k} D={D} n={n}: tuples={N} tester_property={tester_ok} reconstruct={rec_ok} monomials={len(coef)} "
          f"per-coord(d_(i,r)<=k-2 some i>=1)={per_coord} pigeonhole(>=ceil(n/(k-1)) low coords in one x_i)={pigeon} "
          f"x_i-factor in P={inP} dimP={dimP} D^n={D**n} [{time.time()-t0:.1f}s]")
    return tester_ok and rec_ok and per_coord and pigeon and inP

if __name__ == "__main__":
    ok = run(3, 4, 2) and run(4, 4, 2) and run(4, 5, 2)
    print("ALL OK:", ok)

