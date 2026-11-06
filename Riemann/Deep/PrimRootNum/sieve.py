#!/usr/bin/env python3
"""
sieve.py -- exact main terms and costs of lower-bound sieve weights on the first w primes (DERIVATION.md sec 3).

For a level y (integer) the weight set D(y) is generated depth-first: d = q_1 q_2 ... q_r with primes q_1 > q_2 > ... > q_r;
a prime at an ODD position may always be appended, a prime q at an EVEN position m only if  q_1 ... q_{m-1} q^3 < y
(Rosser's linear-sieve rule, beta = 2; integer comparison, exact).  lambda_d = mu(d) = (-1)^r on D(y).
D(y) is prefix-closed and every boundary element has even length, so  sum_{d | n, d in D} mu(d) <= [n = 1]  for every
squarefree n (DERIVATION.md Lemma S1); restricting to divisors of P_w = p_1 ... p_w keeps this property.
Outputs, for every w <= Wmax (q_1 = P^+(d) <= p_w):
    T_w = sum_{d in D, d | P_w} mu(d)/d   as an exact lower bound  (integer numerator, unit 2^-K, K = 64),
    C_w = #{d in D : d | P_w}             (exact integer).
Run:  python3 sieve.py Wmax y1 y2 ... > table.json   (pure Python, one thread)
"""
import sys, math, time, json, bisect

K = 64
ONE = 1 << K

def primes_first(w):
    n = max(30, int(w * (math.log(w + 2) + math.log(math.log(w + 3)) + 3)))
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(range(i * i, n + 1, i)))
    P = [i for i in range(n + 1) if s[i]]
    assert len(P) >= w
    return P[:w]

def enumerate_weights(y, P):
    """per index i of q_1 = P[i]: exact lower-bound numerators (unit 2^-K) of sum mu(d)/d and counts.
    The d = 1 term (mu = +1, value exactly ONE) is NOT included."""
    n = len(P)
    preC = [0] * (n + 1); preF = [0] * (n + 1)
    for i, q in enumerate(P):
        preF[i + 1] = preF[i] + ONE // q            # floor(2^K/q)  (lower bounds)
        preC[i + 1] = preC[i] + (-((-ONE) // q))    # ceil(2^K/q)   (upper bounds)
    Tq = [0] * n; Cq = [0] * n
    sys.setrecursionlimit(100000)
    P3 = [q * q * q for q in P]

    def even_node(prod, li, sgn, i1):
        # node d = prod of even length r >= 2 (counted by the caller); children q = P[j], j < li, at an odd position
        # (no condition), sign -sgn.  Child is a leaf if no grandchild q' < q with prod*q*q'^3 < y exists;
        # sufficient: prod*q*8 >= y.  Leaves: j in [jl, li), P[j] >= ceil(y/(8 prod)).
        thr = -((-y) // (8 * prod))
        jl = bisect.bisect_left(P, thr, 0, li)
        if jl < li:
            if sgn > 0:      # children mu = -1: subtract an upper bound of sum 1/(prod q)
                Tq[i1] -= -((-(preC[li] - preC[jl])) // prod)
            else:            # children mu = +1: add a lower bound
                Tq[i1] += (preF[li] - preF[jl]) // prod
            Cq[i1] += li - jl
        for j in range(jl):
            odd_node(prod * P[j], j, -sgn, i1)

    def odd_node(prod, li, sgn, i1):
        # node d = prod of odd length, mu(d) = sgn; count it
        if sgn > 0:
            Tq[i1] += ONE // prod
        else:
            Tq[i1] -= -((-ONE) // prod)
        Cq[i1] += 1
        # children q = P[j], j < li, at an even position: allowed iff prod*q^3 < y (P3 increasing: binary search)
        lo, hi = 0, li
        while lo < hi:
            mid = (lo + hi) >> 1
            if prod * P3[mid] < y:
                lo = mid + 1
            else:
                hi = mid
        for j in range(lo):
            d = prod * P[j]
            if sgn < 0:      # even node d has mu = -sgn = +1
                Tq[i1] += ONE // d
            else:
                Tq[i1] -= -((-ONE) // d)
            Cq[i1] += 1
            even_node(d, j, -sgn, i1)

    for i in range(n):
        odd_node(P[i], i, -1, i)
    return Tq, Cq

def cumulative(Tq, Cq):
    T = ONE; C = 1; outT = []; outC = []
    for i in range(len(Tq)):
        T += Tq[i]; C += Cq[i]
        outT.append(T); outC.append(C)
    return outT, outC

if __name__ == "__main__":
    Wmax = int(sys.argv[1])
    ys = [int(float(a)) for a in sys.argv[2:]]
    P = primes_first(Wmax)
    res = {"K": K, "Wmax": Wmax, "pmax": P[-1], "levels": {}}
    for y in ys:
        t0 = time.time()
        Tq, Cq = enumerate_weights(y, P)
        T, C = cumulative(Tq, Cq)
        res["levels"][str(y)] = {"T": [str(t) for t in T], "C": C}
        print("y=%d done in %.1f s; C_Wmax=%d T_Wmax=%.6g" % (y, time.time() - t0, C[-1], T[-1] / ONE), file=sys.stderr, flush=True)
    json.dump(res, sys.stdout)
