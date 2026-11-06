#!/usr/bin/env python3
"""q4_qlp.py -- explore16, transfer T2' (N6): Dimitrov-Shapiro's weak q-Laguerre-Polya class (arXiv 2606.17864) for Xi and DH.

B_q(sum a_k z^k/k!) = sum a_k q^{k(k-1)/2} (1-q)^k/(q;q)_k z^k maps LP into LP (their Thm), B_q f -> f as q -> 1, so
RH <=> B_q F_Xi in LP for every q < 1 (Hurwitz), with F(w) = sum gamma(k) w^k/k! (GORZ normalisation).
(a) Hutchinson: b_k^2/(b_{k-1} b_{k+1}) = T_k/q * (1-q^{k+1})/(1-q^k) >= T_k/q, T_k = gamma(k)^2/(gamma(k-1) gamma(k+1)) (Turan;
    >= 1 for Xi by Csordas-Norfolk-Varga). So B_q F in LP for q <= 1/4 whenever T_k >= 1 for all k: test T_k for Xi and DH, k <= 140.
(b) The q-Jensen polynomials J_{n,q}(F; z) = sum C(n,k) gamma(k) delta_k(q) z^k, delta_k = k! q^{k(k-1)/2}(1-q)^k/(q;q)_k:
    certified hyperbolicity for n <= 40 at q = 0.25 ... 0.999, for Xi and DH.
(c) The N6 record and the Laguerre/theta constant carry the same number: Delta_0 = max |arg(-w_rho)| <= atan(2aT/(T^2 - a^2)),
    a = 0.37338196661320752, T = 3e12; Chasse 1/sin^2 Delta_0 versus the cylinder (PSW) threshold 1/(-2 log cos Delta_0).
Run: uv run --no-project --with python-flint --with mpmath python3 -u q4_qlp.py
"""
import time
from math import comb
from flint import arb, ctx
from q1_margins import gam_xi, gam_dh, hyp

def main():
    t0 = time.time()
    N = 140
    ctx.prec = 6000; ctx.cap = 2 * N + 2
    gx = gam_xi(N); gd = gam_dh(N)
    ctx.prec = 1600
    for name, g in (('Xi', gx), ('DH', gd)):
        T = [g[k] ** 2 / (g[k - 1] * g[k + 1]) for k in range(1, N)]
        mn = min(range(len(T)), key=lambda i: float(T[i].mid()))
        allok = all(t > 1 for t in T)
        print(f'{name}: Turan ratios T_k, k = 1..{N - 1}: all > 1: {allok}; min T_{mn + 1} = {T[mn].str(8)}; T_139 = {T[-1].str(8)}')
    print(f'({time.time() - t0:.0f} s)')
    for name, g in (('Xi', gx), ('DH', gd)):
        res = []
        for q in ('0.25', '0.5', '0.8', '0.95', '0.99', '0.999'):
            Q = arb(q); ok = 0; tot = 0
            poch = arb(1); delta = [arb(1)]
            for k in range(1, 41):
                poch *= (1 - Q ** k)
                delta.append(arb(1) * __import__('math').factorial(k) * Q ** (k * (k - 1) // 2) * (1 - Q) ** k / poch)
            for n in range(2, 41):
                tot += 1; ok += hyp([comb(n, k) * g[k] * delta[k] for k in range(n + 1)])
            res.append(f'q={q}: {ok}/{tot}')
        print(f'{name}: q-Jensen J_(n,q), n = 2..40 certified hyperbolic: ' + ', '.join(res))
    import mpmath as mp
    mp.mp.dps = 50
    a = mp.mpf('0.37338196661320752'); T = mp.mpf('3e12')
    D0 = mp.atan(2 * a * T / (T ** 2 - a ** 2))
    print(f'Delta_0 = {mp.nstr(D0, 12)}; Chasse 1/sin^2 = {mp.nstr(1 / mp.sin(D0) ** 2, 12)}; cylinder 1/(-2 log cos) = {mp.nstr(1 / (-2 * mp.log(mp.cos(D0))), 12)}; -log cos Delta_0 = {mp.nstr(-mp.log(mp.cos(D0)), 6)}')
    print(f'done ({time.time() - t0:.0f} s)')

if __name__ == '__main__':
    main()
