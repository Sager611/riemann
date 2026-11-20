# z9_pair.py — explore66: exact collision time (Taylor-coefficient heat flow, 40 digits) for a pair given by indices in zeros_1_20000.npy,
# to check the frozen-model P5 exit time of z8_scale.py at a second Lehmer pair.  usage: z9_pair.py k   (pair = zeros k, k+1; 1-based)
import mpmath as mp, numpy as np, sys, time
mp.mp.dps = 40
T0 = time.time()
LOG = '/Users/adriansager/playground/riemann/Riemann/Deep/ExploreConley/logs/'
k = int(sys.argv[1]); g = np.load(LOG + 'zeros_1_20000.npy')
g1 = mp.findroot(lambda t: mp.siegelz(t), mp.mpf(g[k-1])); g2 = mp.findroot(lambda t: mp.siegelz(t), mp.mpf(g[k]))
c = g1 + g2
def H0(z):
    s = mp.mpf(1)/2 + 1j*z/2
    return s*(s-1)*mp.power(mp.pi, -s/2)*mp.gamma(s/2)*mp.zeta(s)/16
R = mp.mpf('2.0'); M = 96; K = M//2
vals = [H0(c + R*mp.expjpi(2*mp.mpf(j)/M)) for j in range(M)]
a = [sum(vals[j]*mp.expjpi(-2*mp.mpf(j*kk)/M) for j in range(M))/(M*R**kk) for kk in range(K)]
sc = a[2]; a = [mp.re(x/sc) for x in a]
print('pair', k, k+1, 'gamma', mp.nstr(g1, 15), mp.nstr(g2, 15), ' tail |a_k|R^k:', mp.nstr(abs(a[K-1])*R**(K-1), 3), ' time', round(time.time()-T0, 1), flush=True)
fact = [mp.factorial(n) for n in range(K+1)]
def coeffs(t):
    t = mp.mpf(t)
    return [sum(a[m+2*j]*(-t)**j*fact[m+2*j]/(fact[j]*fact[m]) for j in range((K-1-m)//2+1)) for m in range(K)]
def ev(b, w): return mp.polyval(b[::-1], w)
def der(b): return [b[m+1]*(m+1) for m in range(len(b)-1)]
w0 = g2 - g1; u0 = w0**2
S_all = (a[1]/a[0])**2 - 2*a[2]/a[0]; Sfar = S_all - 2/u0
def crit(t):
    b = coeffs(t); ws = mp.findroot(lambda w: ev(der(b), w), mp.mpf(0)); return ev(b, ws)
tc = mp.findroot(crit, -u0/2, tol=mp.mpf(10)**-30)
print('u0 =', mp.nstr(u0, 10), ' S_far =', mp.nstr(Sfar, 8), ' frozen-model t_c =', mp.nstr(-mp.log(1/(1-2*Sfar*u0))/(4*Sfar), 10), ' EXACT t_c =', mp.nstr(tc, 12), ' time', round(time.time()-T0, 1))
