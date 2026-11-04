# t3.py: Hankel-transform check of the Fourier pairing (independent of the spectral-measure formula for K)
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
import core
core.setup(40)
b, h, B = core.consts()
PL = core.laurent_P()
tab = core.read_table()
(L1, C1), (L2, C2) = core.lists_from_table(tab)
t0 = time.time()
M1 = core.Measure(L1, C1, PL, gl_n=64); M2 = core.Measure(L2, C2, PL, gl_n=64)
print('measures', round(time.time()-t0,1), flush=True)
U = 75
def panels(npts):
    xs, ws = core.gl_nodes(npts)
    pts = []
    edges = [mp.mpf(k)/4 for k in range(0, 21)] + [5 + mp.mpf(k)/2 for k in range(1, 2*(U-5)+1)]
    for lo, hi in zip(edges[:-1], edges[1:]):
        mid = (lo+hi)/2; half = (hi-lo)/2
        for x, w in zip(xs, ws):
            pts.append((mid + half*x, half*w))
    return pts
radii = [mp.mpf(x) for x in ['0','0.31','0.77','1','1.5','2.2','3','3.6','4','5.5','7','8.3','9','10.4','12','13.5','15','16.6','18','20','21.7','24']]
for npts in (24, 32):
    pts = panels(npts)
    t1 = time.time()
    # weights e^{-pi h u} p_i(u) and K_i(u) (for the full F) at the quadrature points
    W1 = []; W2 = []; G1 = []; G2 = []
    for u, w in pts:
        e = mp.exp(-mp.pi*h*u) * w
        p1 = core.p_direct(L1, C1, u); p2 = core.p_direct(L2, C2, u)
        W1.append(e*p1); W2.append(e*p2)
        if npts == 32:
            K1 = core.K_derivs(M1, u, 0)[0]; K2 = core.K_derivs(M2, u, 0)[0]
            G1.append(e*(p1 + K2)); G2.append(e*(p2 + K1))   # F and Fhat radial profiles (times weights)
    print('npts', npts, 'points', len(pts), 'prep', round(time.time()-t1,1), flush=True)
    maxrel = [0, 0, 0, 0]
    for s in radii:
        J = [mp.besselj(0, 2*mp.pi*mp.sqrt(u*s)/b) for u, w in pts]
        hk1 = mp.pi/b * mp.fsum(a*j for a, j in zip(W1, J))
        hk2 = mp.pi/b * mp.fsum(a*j for a, j in zip(W2, J))
        e = mp.exp(-mp.pi*h*s)
        K1 = core.K_derivs(M1, s, 0)[0]; K2 = core.K_derivs(M2, s, 0)[0]
        d1 = (hk1 - e*K1); d2 = (hk2 - e*K2)
        line = ' s=%-5s  Hankel[q_p1]-e*K1=%9s (e*K1=%s)  Hankel[q_p2]-e*K2=%9s' % (mp.nstr(s,4), mp.nstr(d1,3), mp.nstr(e*K1,8), mp.nstr(d2,3))
        maxrel[0] = max(maxrel[0], abs(d1)); maxrel[1] = max(maxrel[1], abs(d2))
        if npts == 32:
            p1s = core.p_direct(L1, C1, s) if s not in (0,) and int(s) != s else core.p_derivs(M1, s, 0)[0]
            p2s = core.p_direct(L2, C2, s) if s not in (0,) and int(s) != s else core.p_derivs(M2, s, 0)[0]
            HF = mp.pi/b * mp.fsum(a*j for a, j in zip(G1, J))    # Hankel transform of F
            HFh = mp.pi/b * mp.fsum(a*j for a, j in zip(G2, J))   # Hankel transform of Fhat
            Fh = e*(p2s + K1); F = e*(p1s + K2)
            line += '  Hankel[F]-Fhat=%9s  Hankel[Fhat]-F=%9s (Fhat=%s)' % (mp.nstr(HF-Fh,3), mp.nstr(HFh-F,3), mp.nstr(Fh,8))
            maxrel[2] = max(maxrel[2], abs(HF-Fh)); maxrel[3] = max(maxrel[3], abs(HFh-F))
        print(line, flush=True)
    print('npts', npts, 'max abs diffs', [mp.nstr(x,3) for x in maxrel], 'time', round(time.time()-t1,1), flush=True)

