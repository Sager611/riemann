import os
BASE = os.path.abspath(os.environ.get('RAMAUDIT', '.scratch/ramaudit')) + '/'
import mpmath as mp, numpy as np, time
mp.mp.dps = 20
t0 = time.time()
N = 200000
sieve = np.ones(N+1, dtype=bool); sieve[:2] = False
for i in range(2, int(N**0.5)+1):
    if sieve[i]: sieve[i*i::i] = False
P = np.nonzero(sieve)[0].astype(np.float64); LP = np.log(P)
vs = 0.02; gE = float(mp.euler)
a = 0.2
def integrand(t):
    s = complex(a, t)
    x1 = np.exp((-1+s)*LP); x2 = np.exp((-1-s)*LP); z = 1/P
    Ap = (1 - x1 - x2 + z)*(1 - z)/((1 - x1)*(1 - x2))
    A = np.prod(Ap)
    dlogA = np.sum(-LP*x2*(1/(1 - x2) - 1/(1 - x1 - x2 + z)))
    sm = mp.mpc(a, t)
    z1m = complex(mp.zeta(1 - sm)); z1p = complex(mp.zeta(1 + sm)); dz1p = complex(mp.zeta(1 + sm, derivative=1))
    E2 = np.exp(vs*vs*s*s)
    g = A*E2/(z1m*z1p*s**4)
    D = dlogA + vs*vs*s - dz1p/z1p - 2/s
    return g, g*D
xg, wg = np.polynomial.legendre.leggauss(16)
panels = [(i*0.5, (i+1)*0.5) for i in range(0, 40)] + [(20+5*i, 25+5*i) for i in range(26)]
I0 = 0.0; I1 = 0.0
for lo, hi in panels:
    for xi, wi in zip(xg, wg):
        t = (hi-lo)/2*xi + (hi+lo)/2
        g, gD = integrand(t)
        I0 += wi*(hi-lo)/2*g.real
        I1 += wi*(hi-lo)/2*gD.real
# (1/2pi) int_{-inf}^{inf} = (1/pi) int_0^inf of real part
k1 = I0/np.pi; Hp = I1/np.pi
k0 = Hp + gE*k1
print('a=%.2f  k1=H_c(0)=%.6f  H_c\'(0)=%.6f  k0=%.6f   (%.0fs)' % (a, k1, Hp, k0, time.time()-t0))
print('tail check: |g| at t=150:', abs(integrand(150.0)[0]), ' at t=60:', abs(integrand(60.0)[0]))

