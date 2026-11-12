#!/usr/bin/env python3
"""explore36 / lecam_stats: median Le Cam distance (relative operator norm) over many pairs, per L/L*.
usage: lecam_stats.py <set> <npairs> <M> <fracs> [delta|rel] [cap_spacings]"""
import sys, os, math
import numpy as np
os.environ.setdefault('CAP_SP', sys.argv[6] if len(sys.argv) > 6 else '0.45')
import lecam1 as lc

def main():
    which, npairs, M = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    fracs = [float(f) for f in sys.argv[4].split(',')]
    dmode = sys.argv[5] if len(sys.argv) > 5 else 'rel'
    if which == 'z1e8':
        gam = lc.load_gram(f'{lc.ROOT}/ExploreDPP/logs/zeros_1e8.txt')
        centre = len(gam) // 2
    elif which == 'z1e11':
        gam = lc.load_gram(f'{lc.ROOT}/ExploreDPP/logs/zeros_1e11.txt')
        centre = len(gam) // 2
    else:
        gam = lc.load_zeta_low()
        centre = int(np.searchsorted(gam, float(which[4:])))
    # pairs: spread over the available range keeping M/2 + 2 neighbours on each side
    half = M // 2 + 3
    lo, hi = max(half, centre - 400), min(len(gam) - half - 2, centre + 400)
    i0s = np.unique(np.linspace(lo, hi, npairs).astype(int))
    t = gam[centre]
    Lstar = 0.5 * math.log(t / (2 * math.pi))
    D = 2 * Lstar / (2 * math.pi)
    print(f"# set={which} t~{t:.4e} L*={Lstar:.4f} X*={math.exp(Lstar):.3f} M={M} pairs={len(i0s)} cap={os.environ['CAP_SP']} spacings dmode={dmode}", flush=True)
    for fr in fracs:
        L = fr * Lstar
        rh, ro, mv, gp = [], [], [], []
        for i0 in i0s:
            delta = (0.5 / Lstar) if dmode == 'rel' else float(dmode)
            res, (th, eps, th1, th2) = lc.compensate(gam, int(i0), delta, L, M)
            a = 0.5 * L
            thmax = float(np.max(np.abs(th))) + 1.0
            nx = int(min(700, max(160, 2.0 * thmax * a / math.pi + 100)))
            op0 = lc.opnorm(a, th, np.zeros_like(eps), delta, th1, th2, nx)[0]
            op1 = lc.opnorm(a, th, eps, delta, th1, th2, nx)[0]
            rh.append(res['rho_hs']); ro.append(op1 / op0); mv.append(res['maxeps'] * D); gp.append(res['gap'] * D)
        ro = np.array(ro); rh = np.array(rh)
        q = np.percentile(ro, [10, 50, 90])
        print(f"L/L*={fr:5.2f} c=e^L={math.exp(L):11.3f} rho_op p10/p50/p90 = {q[0]:.2e} {q[1]:.2e} {q[2]:.2e} "
              f"rho_HS p50 {np.median(rh):.2e} maxmove p50 {np.median(mv):.2f} sp; frac(rho_op<1e-3)={np.mean(ro<1e-3):.2f} "
              f"rss={lc.rss_mb():.0f}MB", flush=True)

if __name__ == '__main__':
    main()
