#!/usr/bin/env python3
"""q5audit: independent re-evaluation of the three Question-5 certificate nodes (no code from Degree2Mol).
Reads only the data file q5nodes.json (R_a, R_b, h). Builds Q_a(y) = 1 + sum_n h_n (T_{2n-1}(1-2y) - 1) with h_last reset so that
Q_a(1) = 0 (as the certificate does), and evaluates I_{R_a}(S_a) by mpmath tanh-sinh quadrature (not the exact-sum formula).
Then checks L(R_b) = (2/sqrt15)(I_a + e^{-R_a}) + (1/2 - 3/sqrt15) e^{-R_b} < rho = exp(-(2/3)sqrt(3/5)).
Run: uv run --no-project --with mpmath python3 indep_nodes.py <path to q5nodes.json>
"""
import json, sys
from mpmath import mp, mpf, sqrt, exp, quad, log
mp.dps = 40
nodes = json.load(open(sys.argv[1]))
s15 = sqrt(15); a0 = 2 / s15; b = mpf(1) / 2 - 3 / s15; rho = exp(-(mpf(2) / 3) * sqrt(mpf(3) / 5))
print("b =", mp.nstr(b, 12), " rho =", mp.nstr(rho, 15), " 2sqrt(3/5) =", mp.nstr(2 * sqrt(mpf(3) / 5), 15))
def TU(k, z):  # T_k(z) and T_k'(z) by recurrence (T_k' = k U_{k-1})
    T0, T1 = mpf(1), z; U0, U1 = mpf(1), 2 * z
    if k == 0: return T0, mpf(0)
    Ts = [T0, T1]; Us = [U0, U1]
    for j in range(2, k + 1):
        Ts.append(2 * z * Ts[-1] - Ts[-2]); Us.append(2 * z * Us[-1] - Us[-2])
    return Ts[k], k * Us[k - 1]
for nd in nodes:
    Ra = mpf(nd['R']); h = [mpf(x) for x in nd['h']]
    old = h[-1]; h[-1] = mpf(1) / 2 - sum(h[:-1])
    N = len(h)
    def QdQ(y):
        z = 1 - 2 * y; q = mpf(1); dq = mpf(0)
        # T_{2n-1}(z) and derivative, all n at once
        T = [mpf(1), z]; U = [mpf(1), 2 * z]
        for j in range(2, 2 * N): T.append(2 * z * T[-1] - T[-2]); U.append(2 * z * U[-1] - U[-2])
        for n in range(1, N + 1):
            k = 2 * n - 1; q += h[n - 1] * (T[k] - 1); dq += h[n - 1] * k * U[k - 1] * (-2)
        return q, dq
    q0, _ = QdQ(mpf(0)); q1, _ = QdQ(mpf(1)); qh, _ = QdQ(mpf(1) / 2)
    _, d1 = QdQ(mpf('0.3')); _, d2 = QdQ(mpf('0.7'))
    def integrand(t):
        y = mpf(1) / 2 + t / (2 * Ra); S, dQ = QdQ(y); Sp = dQ / (2 * Ra)
        return exp(t) * (S ** 2 + Sp ** 2) + exp(-t) * ((1 - S) ** 2 + Sp ** 2)
    pts = [Ra * k / 8 for k in range(9)]
    Ia = quad(integrand, pts)
    Rb = None if nd['Rb'] in (None, 'None') else mpf(nd['Rb'])
    L = a0 * (Ia + exp(-Ra)) + (b * exp(-Rb) if Rb is not None else 0)
    print(f"node R_a={nd['R']} R_b={nd['Rb']} N={N} h_last {mp.nstr(old,6)} -> {mp.nstr(h[-1],6)}; Q(0)={mp.nstr(q0,20)} Q(1)={mp.nstr(q1,5)} Q(1/2)={mp.nstr(qh,20)}"
          f" Q'(.3)-Q'(.7)={mp.nstr(d1-d2,3)}")
    print(f"   I_a (quad) = {mp.nstr(Ia, 18)};  L(R_b) = {mp.nstr(L, 15)}  rho - L = {mp.nstr(rho - L, 6)}  ok={L < rho}")
    # kappa at R_a from this S_a, for comparison with the certificate's info line
    c = mpf(1) / 2 - 1 / s15 + a0 * exp(Ra) * Ia; th = sqrt(mpf(3) / 5) / Ra
    print(f"   theta_a = {mp.nstr(th, 10)}  kappa(S_a) = {mp.nstr(1 - log(c) / Ra, 12)}  2theta/3 = {mp.nstr(2 * th / 3, 12)}")
