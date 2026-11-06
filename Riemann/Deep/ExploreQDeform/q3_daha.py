#!/usr/bin/env python3
"""q3_daha.py -- explore16, transfer T3 (N12, N13): the mechanism behind Cherednik's 'Riemann hypothesis for DAHA
superpolynomials' (1709.07589), tested on the torus knots T(n, n+1), n = 2..7.

At a = 0 the DAHA superpolynomial of T(n, n+1) is the q,t-Catalan polynomial C_n (Garsia-Haiman; area/dinv over Dyck paths).
Normalisation (checked against Cherednik's trefoil zeta (1 + q t^2 + a q t)/(1 - t) and his limit S^0 = (1 - t^{2d+2})/(1 - t^2)):
  H0(q,t) = sum over monomials q^a t^b of C_n of q^{d-a} t^b  (d = delta = n(n-1)/2),
  L(q,t)  = H0(q t, t)  (the 'L-function', functional equation t -> 1/(q t)),
  RH      : zeros of L in t on |t| = q^{-1/2}, i.e. zeros of P_q(u) = L(q, u/sqrt q) = sum q^{(d-a-b)/2} u^{d-a+b} on |u| = 1.
At q -> 0 only the a + b = d terms survive: P_0(u) = 1 + u^2 + ... + u^{2d} (simple zeros at roots of unity).
Test: q_c = sup{q : all 2d zeros of P_q on the unit circle}; whether the RH set is an interval; what happens at q_c (a collision
of two unit-circle zeros); RH at q = 2, 3, 4, 5 (possible field sizes; a = 0 specialisation, not the Galkin-Stohr zeta).
Run: uv run --no-project --with numpy python3 -u q3_daha.py   (seconds)
"""
import math
import numpy as np

def dyck_area_sequences(n):
    out = []
    def rec(seq):
        if len(seq) == n:
            out.append(tuple(seq)); return
        for a in range(0, seq[-1] + 2):
            rec(seq + [a])
    rec([0])
    return out

def qt_catalan(n):
    C = {}
    for s in dyck_area_sequences(n):
        area = sum(s)
        dinv = sum(1 for i in range(n) for j in range(i + 1, n) if s[i] == s[j] or s[i] == s[j] + 1)
        C[(area, dinv)] = C.get((area, dinv), 0) + 1
    return C

def unit_circle_count(coef_by_pow, d, q, M=200000):
    """number of sign changes of the real cosine polynomial e^{-i d th} P_q(e^{i th}) on (0, 2 pi): number of simple zeros
    of P_q on the unit circle"""
    th = (np.arange(M) + 0.5) * (2 * math.pi / M)
    T = np.zeros(M)
    for (p, e), c in coef_by_pow.items():
        T += c * q ** (e / 2) * np.cos((p - d) * th)
    s = np.sign(T)
    return int(np.sum(s[1:] * s[:-1] < 0) + (s[-1] * s[0] < 0))

def main():
    print(' n  delta  #paths  limit P_0 ok   q_c (RH for q < q_c)   1/q_c    RH set interval?  RH at q=2,3,4,5')
    for n in range(2, 8):
        C = qt_catalan(n); d = n * (n - 1) // 2
        coef = {}
        for (a, b), c in C.items():
            key = (d - a + b, d - a - b)            # (power of u, exponent of q times 2)
            coef[key] = coef.get(key, 0) + c
        # self-reciprocity check (functional equation)
        fe = all(coef.get((2 * d - p, e), 0) == c for (p, e), c in coef.items())
        lim = sorted(p for (p, e) in coef if e == 0) == list(range(0, 2 * d + 1, 2))
        grid = [k / 200 for k in range(1, 1001)]                  # q in (0, 5]
        ok = [unit_circle_count(coef, d, q, 40000) == 2 * d for q in grid]
        first_bad = next((q for q, o in zip(grid, ok) if not o), None)
        interval = first_bad is None or not any(o for q, o in zip(grid, ok) if q > first_bad)
        if first_bad is None:
            qc = float('inf')
        else:
            lo, hi = first_bad - 1 / 200, first_bad
            for _ in range(40):
                mid = (lo + hi) / 2
                if unit_circle_count(coef, d, mid) == 2 * d: lo = mid
                else: hi = mid
            qc = lo
        phys = ''.join('Y' if unit_circle_count(coef, d, q) == 2 * d else 'n' for q in (2, 3, 4, 5))
        print(f'{n:2d} {d:5d} {sum(C.values()):6d}   {lim and fe!s:>5}        {qc:14.6f}     {1 / qc:8.5f}   {interval!s:>5}            {phys}')
        if first_bad is not None:
            # zeros just past q_c: which pair left the circle?
            P = np.zeros(2 * d + 1)
            for (p, e), c in coef.items(): P[p] += c * (qc * 1.001) ** (e / 2)
            z = np.roots(P[::-1])
            off = sorted(abs(abs(x) - 1) for x in z)[-4:]
            zc = np.roots(np.array([sum(c * qc ** (e / 2) for (p, e), c in coef.items() if p == k) for k in range(2 * d, -1, -1)]))
            mind = min(abs(zc[i] - zc[j]) for i in range(len(zc)) for j in range(i + 1, len(zc)))
            print(f'      at 1.001 q_c: largest | |u| - 1 | = {[round(x, 5) for x in off]};  min zero distance at q_c = {mind:.2e}')

if __name__ == '__main__':
    main()
