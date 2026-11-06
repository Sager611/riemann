#!/usr/bin/env python3
"""g_screen.py -- explore15: the screening law behind the excess over Chasse's bound (node N6).

Heuristic.  To first order in 1/d, J^{d,0}(X/d) = sum_k t_k X^k prod_{i<k}(1 - i/d) ~ Theta(X) - (1/(2d)) X^2 Theta''(X).
With X = z^2 and z = i x near a zero at height gamma this is the de Bruijn-Newman flow d_t H = -H_xx run for the
time t = gamma^2/(8d).  Under that flow an isolated pair x +- i delta becomes real at t = delta^2/2, i.e. at
d = gamma^2/(4 delta^2) = Chasse's D_sec (leading order).  Real neighbours x_k add attraction:
y' = -1/y - 2 y S,  S = sum_k 1/(x - x_k)^2  (y -> 0 limit; the mirrors -x_k and -x - i y included), so the pair becomes real at
t_c = log(1 + 2 S delta^2)/(4 S), and the first failing degree is d* ~ D_sec * 2 S delta^2 / log(1 + 2 S delta^2).
Here: S from the actual zeros (Arb) for the toys of c_reach.py, compared with the measured d*; then the
same law at the N6 record's height T = 3.000175e12 with the mean neighbour sum S = rho^2 I_GUE, rho = log(T/2pi)/2pi.
Run: uv run --no-project --with python-flint --with mpmath python3 -u g_screen.py   (seconds)
"""
import math
import mpmath as mp
from flint import acb, ctx

def main():
    ctx.prec = 64
    zz = [float(z.imag.mid()) for z in acb.zeta_zeros(1, 400)]
    rho = lambda t: math.log(t / (2 * math.pi)) / (2 * math.pi)
    def S_of(k):
        g = zz[k - 1]
        s = sum(1 / (g - x) ** 2 + 1 / (g + x) ** 2 for j, x in enumerate(zz) if j != k - 1)
        s += 1 / (2 * g * g)                           # mirror zero -gamma - i delta: Im 2/(2 gamma + 2 i y) = -y/gamma^2
        top = zz[-1]                                   # tail beyond the 400th zero by the density
        s += mp.quad(lambda x: rho(x) / (x - g) ** 2, [top, mp.inf])
        return float(s)
    measured = {1: 249, 2: 558, 3: 795, 4: 1197, 5: 1412, 7: 2207, 10: None}
    import sys
    if len(sys.argv) > 1: measured[10] = int(sys.argv[1])
    print(' k  gamma_k   S       D_sec   predicted d*   measured d*   (delta = 0.45)')
    for k in [1, 2, 3, 4, 5, 7, 10]:
        g = zz[k - 1]; d = 0.45
        S = S_of(k)
        D = int(1 / math.sin(2 * math.atan(d / g)) ** 2)
        a = 2 * S * d * d
        pred = D * a / math.log1p(a)
        print(f' {k:2d} {g:8.3f} {S:7.4f} {D:6d}   {pred:9.0f} (+{100 * (pred / D - 1):.1f}%)   {measured[k]}'
              + (f' (+{100 * (measured[k] / D - 1):.1f}%)' if measured[k] else ''))
    # mean neighbour sum for GUE statistics: S = rho^2 * I, I = int_R (1 - sinc^2(pi u))/u^2 du
    I = 2 * mp.quad(lambda u: (1 - (mp.sin(mp.pi * u) / (mp.pi * u)) ** 2) / u ** 2, [0, 1, 10, mp.inf])
    print(f'GUE: I = int (1 - (sin pi u/pi u)^2)/u^2 du = {mp.nstr(I, 8)}  (= 2 pi^2/3 = {mp.nstr(2 * mp.pi ** 2 / 3, 8)}?)')
    T = 3000175332800.0; a8 = 3 / 8
    r = rho(T)
    for name, S in [('zeros on both sides', r * r * float(I)), ('zeros below only', r * r * float(I) / 2)]:
        a = 2 * S * a8 * a8
        f = a / math.log1p(a)
        D = (16 * T * T - 4) / 9
        print(f'N6 record height T = {T:.6e}, rho = {r:.3f}, {name}: S = {S:.1f}, 2 S delta^2 = {a:.2f},'
              f' heuristic reach = {f:.2f} x (16T^2-4)/9 = {f * D:.3e}')

if __name__ == '__main__':
    main()
