#!/usr/bin/env python3
"""e_li_reach.py -- explore15, transfer C (N3, N6): Barvinok interpolation = Taylor-coefficient criteria; reach.

Li: lambda_n = sum_rho [1 - (1 - 1/rho)^n] are the Taylor coefficients of d/dz log xi(1/(1-z)) (Barvinok's
interpolation object for the disc |z| < 1 = {Re s > 1/2}).  An off-line zero rho = 1/2 - delta + i T maps to
|1 - 1/rho|^2 = 1 + 2 delta / ((1/2 - delta)^2 + T^2) > 1, so it is seen only at n ~ (T^2/(2 delta)) log n.
Toy (as c_reach.py): zero pair 1/2 +- i gamma_k replaced by the quartet 1/2 +- delta +- i gamma_k.
lambda_n(RH part) ~ (n/2)(log n - 1 - log 2 pi + gamma_E) (main term; the remainder is O(sqrt n log n) under RH).
Also: the worst-case reach of a zero just above the verified height T = 3.000175e12 with Re rho <= 7/8, for
Li (compare the N3 record 77.03 T^2) and for Jensen (Chasse: (16 T^2 - 4)/9, the N6 record).
Run: uv run --no-project --with mpmath python3 -u e_li_reach.py   (seconds)
"""
import mpmath as mp
mp.mp.dps = 40

def main_term(n):
    return n / 2 * (mp.log(n) - 1 - mp.log(2 * mp.pi) + mp.euler)

def delta_n(n, g, d):
    """lambda_n(toy) - lambda_n(zeta): quartet 1/2 +- d +- i g minus pair 1/2 +- i g"""
    tot = mp.mpf(0)
    for r in [mp.mpc(0.5 + d, g), mp.mpc(0.5 - d, g)]:
        tot += 2 * mp.re(1 - (1 - 1 / r) ** n)
    tot -= 2 * mp.re(1 - (1 - 1 / mp.mpc(0.5, g)) ** n)
    return tot

def first_negative(g, d):
    """first n >= 50 (the main term is positive from n = 9 on) with main_term(n) + Delta_n < 0"""
    n, step = 50, 64
    while True:
        if main_term(n) + delta_n(n, g, d) < 0:
            if step == 1: return n
            n = max(50, n - step); step //= 4
            continue
        n += step

def main():
    zeros = {1: mp.mpf('14.134725141734693790'), 2: mp.mpf('21.022039638771554993'), 3: mp.mpf('25.010857580145688763')}
    print(' k delta | Li: first n with lambda_n(toy) < 0 | (T^2/delta) log(main/2) | Jensen D_sec = 1/sin^2(2 atan(delta/T)) (c_reach: d* = D_sec + 2..23)')
    for k, d in [(1, 0.45), (1, 0.35), (1, 0.25), (1, 0.1), (2, 0.45), (3, 0.45)]:
        g = zeros[k]; d = mp.mpf(d)
        n = first_negative(g, d)
        pred = ((g ** 2 + (0.5 - d) ** 2) / d) * mp.log(main_term(n) / 2)   # |1-1/rho|^n = e^{n delta/|rho|^2 (1+o(1))}
        D = int(1 / mp.sin(2 * mp.atan(d / g)) ** 2)
        print(f' {k} {float(d):.2f}  | {n:8d} | {float(pred):10.0f} | {D}')
    # worst case above the verified height
    T = mp.mpf(3000175332800); a = mp.mpf(3) / 8
    r = mp.mpc(0.5 - a, T)
    lg = mp.log(abs(1 - 1 / r))                      # growth rate per n
    # solve n/2 (log n - 1 - log 2pi + gamma) = 2 exp(n lg)
    f = lambda n: mp.log(main_term(n) / 2) - n * lg
    n = mp.findroot(f, mp.mpf(80) * T ** 2)
    print(f' worst case Re rho = 7/8 just above T = {int(T)}: Li reach n* = {mp.nstr(n, 6)} = {mp.nstr(n / T ** 2, 5)} T^2'
          f'  (N3 record: 77.03 T^2 from RH to T; with the 7/8 half-plane ours 105.77 T^2)')
    # many zeros just above T, all at Re rho = 7/8 with aligned phases: K = N(2T) - N(T) ~ (T/2pi) log(T/2pi)
    Kz = T / (2 * mp.pi) * mp.log(T / (2 * mp.pi))
    f2 = lambda n: mp.log(main_term(n) / 2) - mp.log(Kz) - n * lg
    n2 = mp.findroot(f2, mp.mpf(80) * T ** 2)
    print(f'   same with K = (T/2pi) log(T/2pi) = {mp.nstr(Kz, 3)} such zeros in phase: n* = {mp.nstr(n2 / T ** 2, 5)} T^2')
    D = (16 * T ** 2 - 4) / 9
    print(f'   Jensen reach (Chasse, sharp up to O(1) by c_reach.py): (16T^2 - 4)/9 = {mp.nstr(D, 6)} = {mp.nstr(D / T ** 2, 5)} T^2 (N6 record)')

if __name__ == '__main__':
    main()
