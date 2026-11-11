
# sanity2: Theta^{(n)}(-s)/gamma(n) = E_n[0F1(;n+1/2;-t^2 s/4)] on a fine grid of [0, R_n] (R_n = 0.9 nu^2/(u+1)^2), and the
# first sign change beyond R_n on a 3% grid. Small n only: the theorem needs u >= 12; this only tests the scale. Not part of the proof.
import mpmath as mp
mp.mp.dps = 18
pi = mp.pi
def phi(u):
    e2 = mp.e**(2*u)
    return sum((2*pi**2*k**4*mp.e**(mp.mpf(9)/2*u)-3*pi*k**2*mp.e**(mp.mpf(5)/2*u))*mp.e**(-pi*k**2*e2) for k in range(1,4))
def ustar(n): return mp.findroot(lambda u: u*(pi*mp.e**(2*u)-mp.mpf(9)/4)-n, mp.log(n)/2-1)
theta = mp.mpf('0.9')
import sys
NOSCAN = '--noscan' in sys.argv
NS = [int(x) for x in sys.argv[1:] if x != '--noscan'] or [50, 200]
for n in NS:
    u = ustar(n); nu = n-mp.mpf(1)/2
    w_ = 12*mp.sqrt(u/(4*n))+mp.mpf('0.05')
    lo, hi = max(mp.mpf('0.02'), u-max(w_, mp.mpf('0.15'))-0.1), u+max(w_, mp.mpf('0.15'))
    lpk = 2*n*mp.log(u)+mp.log(phi(u))
    dens = lambda t: mp.e**(2*n*mp.log(t)+mp.log(phi(t))-lpk)
    nodes = mp.linspace(lo, hi, 31)
    M = mp.quad(dens, nodes)
    F = lambda s: mp.quad(lambda t: dens(t)*mp.hyp0f1(n+mp.mpf(1)/2, -t**2*s/4), nodes)/M
    Rn = theta*nu**2/(u+1)**2
    grid = [Rn*k/60 for k in range(1, 61)]
    vals = [F(s) for s in grid]
    mono = all(vals[i+1] < vals[i] for i in range(len(vals)-1))
    print('n=%d u*=%s R_n=%s: F>0 on 60-point grid of (0,R_n]: %s, decreasing: %s, F(R_n)=%s' % (n, mp.nstr(u,6), mp.nstr(Rn,6), all(v > 0 for v in vals), mono, mp.nstr(vals[-1],5)), flush=True)
    s, prev, z = Rn, vals[-1], None
    if NOSCAN: continue
    while s < 15*Rn:
        s2 = s*mp.mpf('1.03'); v = F(s2)
        if v*prev < 0: z = (s, s2); break
        s, prev = s2, v
    print('   first sign change on the 3%% grid beyond R_n: %s (ratio %s)' % (None if z is None else [mp.nstr(z[0],7), mp.nstr(z[1],7)], None if z is None else mp.nstr(z[0]/Rn,4)), flush=True)
