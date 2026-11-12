# Lemma A part 1 checks: closed form I(n) = sum_{k odd<n}(n-k)/k - (n/2)log 4pi vs t(n); int r = -1/8; lambda_n^(0) by quadrature.
import mpmath as mp
mp.mp.dps = 25
kap = 1 + mp.log(2*mp.pi) - mp.euler
t = lambda n: mp.mpf(n)/2*(mp.log(n) - kap)
def I(n):
    return mp.fsum(mp.mpf(n - k)/k for k in range(1, n, 2)) - mp.mpf(n)/2*mp.log(4*mp.pi)
for n in (10, 11, 100, 101, 1000, 1001):
    print('n=%5d  n*(I(n)-t(n)) = %s' % (n, mp.nstr(n*(I(n) - t(n)), 8)))
for T in (10**3, 10**6):
    print('T=%g: theta_RS(T)/pi - (T/2pi)log(T/2pi e) + 1/8 = %s' % (T, mp.nstr(mp.siegeltheta(T)/mp.pi - T/(2*mp.pi)*mp.log(T/(2*mp.pi*mp.e)) + mp.mpf(1)/8, 5)))
def lam0(n):
    f = lambda v: 2*(1 - mp.cos(n*v))*mp.siegeltheta(mp.cot(v/2)/2, 1)/mp.pi/(4*mp.sin(v/2)**2)
    pts = [mp.pi*j/(4*n) for j in range(4*n + 1)]
    return mp.quad(f, pts)
for n in (5, 10, 20, 50):
    d = lam0(n) - t(n)
    print('n=%3d  lambda0 - t = %s ; +1/4+1/(24n) = %s' % (n, mp.nstr(d, 12), mp.nstr(d + mp.mpf(1)/4 + mp.mpf(1)/(24*n), 4)))
