# D: explore36 threshold test near r = log X / L* in [0.9, 1.0], L* = (1/2) log(t/2pi), using explore36's own compensate/opnorm
import sys, os, math, time, json
sys.path.insert(0, '/Users/adriansager/playground/riemann/Riemann/Deep/ExploreMinimax')
import numpy as np
import lecam1 as LC
gam = LC.load_zeta_low()
t0 = float(os.environ.get('T0', '1.0e4'))
pairs = [int(np.searchsorted(gam, t0)) + k for k in (0, 7)]
rs = [float(x) for x in os.environ.get('RS', '0.90,0.95,1.00').split(',')]
Ms = [int(x) for x in os.environ.get('MS', '100,400,800').split(',')]
print(f"zeros loaded: {len(gam)}, max {gam[-1]:.1f}; pairs at indices {pairs}, CAP_SP={os.environ.get('CAP_SP','0.45')}", flush=True)
for i0 in pairs:
    gbar = 0.5*(gam[i0] + gam[i0+1]); Ls = 0.5*math.log(gbar/(2*math.pi)); delta = 1/(2*Ls)
    for r in rs:
        L = r*Ls
        for M in Ms:
            ts = time.time()
            res, (th, eps, th1, th2) = LC.compensate(gam, i0, delta, L, M)
            a = L/2
            thmax = float(np.max(np.abs(th))) + 1.0
            nx = int(min(900, max(200, 2.5*thmax*L/math.pi + 100)))
            op1 = LC.opnorm(a, th, eps, delta, th1, th2, nx)[0]
            op0 = LC.opnorm(a, th, np.zeros_like(eps), delta, th1, th2, nx)[0]
            print(f"t={gbar:.2f} gap={res['gap']:.3f} r={r:.2f} M={res['M']} delta={delta:.4f} rho_hs={res['rho_hs']:.3e} rho_op={op1/op0:.3e} maxeps/sp={res['maxeps']*math.log(gbar/(2*math.pi))/(2*math.pi):.3f} nx={nx} sec={time.time()-ts:.1f} rssMB={LC.rss_mb():.0f}", flush=True)
