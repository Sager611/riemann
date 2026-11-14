#!/usr/bin/env python3
"""Exact exponent bookkeeping for the Bombieri-Iwaniec/Bourgain bound (Bourgain 1408.5794, (3.12)-(3.18)), Q = R binding.
Exponents of T: M = T^a, N = T^n, R = T^((3a-n-1)/2), H = N/R.
 trivial  E_t = 6a - 3n            (|S| << M/N^{1/2})
 diagonal E_d = a + 2n + 3/2        (B1 >= |I| ~ M/N)
 off-diag E_o = 9a/2 - 3n/2 + 1/2 - (s/4)(3n - 3a + 1)   (B1_off << [H-W] * H^{-s})
 V-route (x2 tightened by V = N/Q ~ H): diagonal inflated by V: E_dV = E_d + (3n - 3a + 1)/4.
|S| << M^{1/2} T^theta with 6 theta = min_n max(E) - 3a."""
from fractions import Fraction as F
def best(a, s, vroute=False, grid=None):
    # exact optimum: the max of linear functions of n is minimized at an intersection point or an endpoint
    lo, hi = a - F(1, 3), (3*a - 1)/2            # R <= N <= R^2
    Et = lambda n: 6*a - 3*n
    Ed = lambda n: a + 2*n + F(3, 2) + ((3*n - 3*a + 1)/4 if vroute else 0)
    Eo = lambda n: F(9, 2)*a - F(3, 2)*n + F(1, 2) - s*(3*n - 3*a + 1)/4
    fs = [Et, Ed, Eo]
    cands = [lo, hi]
    for i in range(3):
        for j in range(i+1, 3):
            # f_i(n) = f_j(n): both linear in n
            fi0, fi1 = fs[i](F(0)), fs[i](F(1)) - fs[i](F(0)); fj0, fj1 = fs[j](F(0)), fs[j](F(1)) - fs[j](F(0))
            if fi1 != fj1:
                n = (fj0 - fi0)/(fi1 - fj1)
                if lo <= n <= hi: cands.append(n)
    n = min(cands, key=lambda n: max(f(n) for f in fs))
    E = max(f(n) for f in fs)
    return n, (E - 3*a)/6
print("theta(sigma) at a = 1/2 (V-free route); closed form (26+5s)/(12(14+3s)):")
for s in [F(0), F(1, 2), F(2, 3), F(3, 4), F(1), F(3, 2), F(2)]:
    n, th = best(F(1, 2), s)
    print(f"  s={str(s):>4}  n*={str(n):>8}  theta={str(th):>8} = {float(th):.6f}  closed={str((26+5*s)/(12*(14+3*s)))}")
print("alpha-independence check (s = 1): ", [str(best(a, F(1))[1]) for a in [F(7, 17), F(9, 20), F(1, 2)]])
print("V-route at a = 1/2 (Huxley resonance curves, V = N/Q): theta_V(s); 13/84 =", float(F(13, 84)))
for s in [F(0), F(2, 3), F(3, 4), F(1), F(2)]:
    n, th = best(F(1, 2), s, vroute=True)
    print(f"  s={str(s):>4}  n*={str(n):>8}  theta_V={str(th):>8} = {float(th):.6f}")
