# explore18 test T1: the class-D (Kitaev / Pfaffian-parity) Z2 invariant on windows = Hardy's sign-change count.
#   nu(I) = sign A(a) * sign A(b) = (-1)^{#zeros of A in I}, A = F(1/2+it) real.
#   Windows chosen without knowing the zeros: Gram intervals (Gamma-factor phase theta alone) and equal windows of
#   width w * (local mean spacing). Certified zeros = number of odd windows. Bulk count N(T): argument principle
#   (mpmath nzeros for zeta; rectangle count for DH).
#   zeta: Z(t) by the Riemann-Siegel formula (main sum + C0 + C1 terms) in numpy, checked against mpmath siegelz.
import sys, time
import numpy as np
from topo_common import *

def theta_np(t):
    return t / 2 * np.log(t / (2 * np.pi)) - t / 2 - np.pi / 8 + 1 / (48 * t) + 7 / (5760 * t ** 3)

def Z_rs(t):
    t = np.asarray(t, dtype=float)
    a = np.sqrt(t / (2 * np.pi)); N = np.floor(a).astype(int); p = a - N
    th = theta_np(t)
    out = np.zeros_like(t)
    for n in range(1, int(N.max()) + 1):
        m = n <= N
        out[m] += 2 * np.cos(th[m] - t[m] * np.log(n)) / np.sqrt(n)
    c0 = np.cos(2 * np.pi * (p * p - p - 1 / 16)) / np.cos(2 * np.pi * p)
    # C1 = -(1/(96 pi^2)) Psi'''(p), Psi(p) = cos(2 pi (p^2 - p - 1/16))/cos(2 pi p); numerical third derivative
    h = 1e-3
    Psi = lambda q: np.cos(2 * np.pi * (q * q - q - 1 / 16)) / np.cos(2 * np.pi * q)
    d3 = (Psi(p + 2 * h) - 2 * Psi(p + h) + 2 * Psi(p - h) - Psi(p - 2 * h)) / (2 * h ** 3)
    c1 = -d3 / (96 * np.pi ** 2)
    sgn = np.where(N % 2 == 1, 1.0, -1.0)
    out += sgn * (t / (2 * np.pi)) ** (-0.25) * (c0 + c1 * (t / (2 * np.pi)) ** (-0.5))
    return out

def gram_np(T):
    gs = []; n = -1
    while True:
        g = 2 * np.pi * np.exp(1 + np.real(complex(mp.lambertw((8 * n + 1) / (8 * np.e)))))
        for _ in range(30):
            g -= (theta_np(g) - n * np.pi) / (0.5 * np.log(g / (2 * np.pi)))
        if g > T: break
        if g >= 10: gs.append(g)
        n += 1
    return np.array(gs)

def odd_windows(signs_at, edges):
    s = signs_at(edges)
    return int(np.sum(s[:-1] * s[1:] < 0)), len(edges) - 1

def main(which):
    widths = [0.25, 0.5, 1.0, 2.0]
    if which == 'xi':
        T = 2000.0; T0 = time.time()
        chk = [(t, float(mp.siegelz(t)), float(Z_rs(np.array([t]))[0])) for t in (20.5, 101.3, 555.5, 1999.1)]
        print("RS check (t, mpmath siegelz, numpy RS):", [(a, round(b, 6), round(c, 6)) for a, b, c in chk])
        tg = np.arange(10.0, T, 0.002); zg = Z_rs(tg)
        n0 = int(np.sum(np.sign(zg[:-1]) * np.sign(zg[1:]) < 0))
        NT = int(mp.nzeros(T))
        print(f"zeta on [10, {T:.0f}]: N_0 (sign changes, step 0.002) = {n0}; N(T) (mpmath nzeros, argument principle) = {NT}  ({time.time()-T0:.0f} s)")
        sg = lambda x: np.sign(Z_rs(np.asarray(x)))
        G = gram_np(T)
        odd, nw = odd_windows(sg, G)
        print(f"   Gram intervals: {nw} windows, odd parity {odd} -> certified fraction {odd / NT:.4f}")
        for w in widths:
            e = [10.0]
            while True:
                u = e[-1] + w * 2 * np.pi / np.log(e[-1] / (2 * np.pi))
                if u > T: break
                e.append(u)
            odd, nw = odd_windows(sg, np.array(e))
            print(f"   equal windows {w:4.2f} x mean spacing: {nw} windows, odd parity {odd} -> certified fraction {odd / NT:.4f}")
        ia = int(np.searchsorted(G, float(GA))) - 1; ib = int(np.searchsorted(G, float(GB))) - 1
        print(f"   moved pair: gamma_6 = {float(GA):.4f} in Gram interval [{G[ia]:.4f}, {G[ia+1]:.4f}], gamma_7 = {float(GB):.4f} in [{G[ib]:.4f}, {G[ib+1]:.4f}]:"
              f" {'same' if ia == ib else 'different'} intervals; for mv both lose their odd parity (2 certified zeros fewer), while N(T) is unchanged")
    elif which == 'dh':
        mp.mp.dps = 15
        T = 120; T0 = time.time()
        tg = [1 + mp.mpf(k) / 20 for k in range(int((T - 1) * 20) + 1)]
        sv = [mp.sign(mp.re(dh(mp.mpc(0.5, t)))) for t in tg]
        n0 = sum(1 for k in range(len(sv) - 1) if sv[k] * sv[k + 1] < 0)
        print(f"DH on [1, {T}]: N_0 (sign changes, step 0.05) = {n0}  ({time.time()-T0:.0f} s)", flush=True)
        nr = zeros_in_rect(dh, 0.505, 1.6, 1, T, 0.25)
        NT = n0 + 2 * int(mp.nint(nr))
        print(f"   zeros with 0.505 < Re s < 1.6: {mp.nstr(nr, 6)}  ->  N_total = N_0 + 2 * that = {NT}  ({time.time()-T0:.0f} s)", flush=True)
        th = lambda t: t / 2 * mp.log(5 / mp.pi) + mp.im(mp.loggamma(mp.mpf(3) / 4 + 1j * t / 2))
        grid = [1 + mp.mpf(k) / 20 for k in range(int((T - 1) * 20) + 1)]
        tv = [th(t) for t in grid]; grams = []
        for k in range(len(grid) - 1):
            lo, hi = mp.floor(tv[k] / mp.pi), mp.floor(tv[k + 1] / mp.pi)
            if hi > lo and tv[k + 1] > tv[k]:
                grams.append(mp.findroot(lambda t: th(t) - hi * mp.pi, (grid[k], grid[k + 1]), solver='bisect'))
        gs = [mp.sign(mp.re(dh(mp.mpc(0.5, g)))) for g in grams]
        odd = sum(1 for k in range(len(gs) - 1) if gs[k] * gs[k + 1] < 0)
        print(f"   Gram intervals: {len(gs) - 1} windows, odd parity {odd} -> certified fraction {odd / NT:.4f}")
        idx = lambda t: int(mp.nint((t - 1) * 20))
        for w in widths:
            e = [mp.mpf(1)]
            while True:
                u = e[-1] + w * 2 * mp.pi / mp.log(5 * max(e[-1], 3) / (2 * mp.pi))
                if u > T: break
                e.append(u)
            s = [sv[idx(x)] for x in e]
            odd = sum(1 for k in range(len(s) - 1) if s[k] * s[k + 1] < 0)
            print(f"   equal windows {w:4.2f} x mean spacing: {len(s) - 1} windows, odd parity {odd} -> certified fraction {odd / NT:.4f}")
        print(f"   ({time.time()-T0:.0f} s)")

if __name__ == '__main__':
    main(sys.argv[1])
