#!/usr/bin/env python3
"""explore36 / lecam_dh: Le Cam distance for the 9 real off-line pairs of Davenport-Heilbronn (conductor 5) below 420,
against on-line alternatives (pair replaced by two on-line zeros one mean spacing apart + moved neighbours);
and for zeta pairs moved off the line by the same delta at the same heights.  usage: lecam_dh.py M fracs"""
import sys, os, json, math
import numpy as np
os.environ.setdefault('CAP_SP', '1.0')
import lecam1 as lc

def run(gam, i0, delta, Lstar, fracs, M, tag, cond):
    os.environ['COND'] = str(cond)
    out = []
    for fr in fracs:
        L = fr * Lstar
        res, (th, eps, th1, th2) = lc.compensate(gam, i0, delta, L, M)
        a = 0.5 * L
        thmax = float(np.max(np.abs(th))) + 1.0
        nx = int(min(700, max(160, 2.0 * thmax * a / math.pi + 100)))
        op0 = lc.opnorm(a, th, np.zeros_like(eps), delta, th1, th2, nx)[0]
        op1 = lc.opnorm(a, th, eps, delta, th1, th2, nx)[0]
        out.append(op1 / op0)
    return out

def main():
    M = int(sys.argv[1]); fracs = [float(f) for f in sys.argv[2].split(',')]
    d = json.load(open(f'{lc.ROOT}/ExploreHelio/results/dh_zeros.json'))
    on = np.array(sorted(d['online']))
    zz = lc.load_zeta_low()
    print(f"# DH online zeros {len(on)} to {d['TMAX']}, offline pairs {len(d['offline'])}; M={M}; rho_op at L/L* = {fracs}")
    for beta, g0 in d['offline']:
        delta = beta - 0.5
        Ddh = math.log(5 * g0 / (2 * math.pi)) / (2 * math.pi)
        Ls_dh = 0.5 * math.log(5 * g0 / (2 * math.pi))
        # on-line alternative: two zeros at g0 -+ 1/(2 Ddh)
        gam = np.sort(np.concatenate([on, [g0 - 0.5 / Ddh, g0 + 0.5 / Ddh]]))
        i0 = int(np.argmin(np.abs(gam - (g0 - 0.5 / Ddh))))
        r_dh = run(gam, i0, delta, Ls_dh, fracs, M, 'DH', 5)
        # zeta: move the consecutive pair nearest g0 off the line by the same delta
        j0 = int(np.searchsorted(zz, g0)) - 1
        Ls_z = 0.5 * math.log(g0 / (2 * math.pi))
        r_z = run(zz, j0, delta, Ls_z, fracs, M, 'zeta', 1)
        print(f"t={g0:8.3f} delta={delta:.4f}  DH: L*={Ls_dh:.3f} X*={math.exp(Ls_dh):6.2f} rho_op " + ' '.join(f'{x:.1e}' for x in r_dh)
              + f"  | zeta moved pair (t={zz[j0]:.2f}): L*={Ls_z:.3f} X*={math.exp(Ls_z):5.2f} rho_op " + ' '.join(f'{x:.1e}' for x in r_z), flush=True)
    print(f"rss={lc.rss_mb():.0f}MB")

if __name__ == '__main__':
    main()
