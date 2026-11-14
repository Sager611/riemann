# Sign of the printed Dimitrov-Xu off-axis transform  R(x,y) = Re(Xi'(z)^2 - Xi(z)Xi''(z)), z = x+iy,
# = Re(xi(s)^2 * (xi'/xi)'(s)), s = 1/2 - y + i x.  Printed Thm 1.1 (arXiv 1606.05011 l.198) needs R(.,y) zero-free on R
# for every y in (0,1/2) under RH.  Float-level scan (mpmath), not a certificate.
import sys, mpmath as mp
mp.mp.dps = 25
def R(x, y):
    s = mp.mpf(1)/2 - y + 1j*x
    if mp.re(s) < 0.5:          # use xi(s) = xi(1-s): (xi'/xi)'(s) = (xi'/xi)'(1-s), xi(s)^2 = xi(1-s)^2
        s = 1 - s
    z0 = mp.zeta(s); z1 = mp.zeta(s, 1, 1); z2 = mp.zeta(s, 1, 2)
    L2 = -1/s**2 - 1/(s-1)**2 + mp.psi(1, s/2)/4 + (z2/z0 - (z1/z0)**2)
    lx = mp.log(s) + mp.log(s-1) - s/2*mp.log(mp.pi) + mp.loggamma(s/2) + mp.log(z0)
    ph = mp.exp(2j*mp.im(lx))
    L1 = 1/s + 1/(s-1) - mp.log(mp.pi)/2 + mp.psi(0, s/2)/2 + z1/z0
    J = abs(L1)**2 + mp.re(L2 + L1**2)   # Jensen's J/|xi|^2, must be >= 0 under RH
    global JMIN; JMIN = min(JMIN, float(J))
    return mp.re(ph*L2) / abs(L2)      # sign-faithful, normalized
x0 = float(sys.argv[1]); x1 = float(sys.argv[2]); n = int(sys.argv[3]); ys = [float(v) for v in sys.argv[4].split(',')]
for y in ys:
    JMIN = 1e99
    prev = None; ch = []; mn = 9; 
    for k in range(n + 1):
        x = x0 + (x1 - x0) * k / n
        v = float(R(x, y)); mn = min(mn, v)
        if prev is not None and (v > 0) != (prev > 0): ch.append(round(x, 4))
        prev = v
    print('y=%.3f x in [%g,%g] n=%d  min=%.4f  sign changes=%d first=%s  Jensen J/|xi|^2 min=%.4g' % (y, x0, x1, n, mn, len(ch), ch[:4], JMIN), flush=True)
