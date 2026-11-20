# z1_lehmer.py — explore66: heat flow of H_t near the first Lehmer pair (zeros 6709, 6710 of zeta), via the EXACT
# linear action of e^{-t d^2} on Taylor coefficients: H_0(c+w) = sum a_k w^k  =>  H_t(c+w) = sum_m b_m(t) w^m,
# b_m(t) = sum_j a_{m+2j} (-t)^j (m+2j)!/(j! m!).  a_k from M samples of H_0 on |w| = R (trapezoid = DFT).
# H_0(z) = xi(1/2 + i z/2)/8, xi(s) = s(s-1)/2 pi^{-s/2} Gamma(s/2) zeta(s).  Zeros of H_0: z = 2 gamma.
import mpmath as mp, json, time, sys
mp.mp.dps = 40
T0 = time.time()
LOG = '/Users/adriansager/playground/riemann/Riemann/Deep/ExploreConley/logs/'
zs = json.load(open(LOG + 'zeros_6700_6720.json')); g = {n: mp.mpf(s) for n, s in zs}
c = g[6709] + g[6710]          # midpoint in z = 2 gamma units
def H0(z):
    s = mp.mpf(1)/2 + 1j*z/2
    return s*(s-1)*mp.power(mp.pi, -s/2)*mp.gamma(s/2)*mp.zeta(s)/16
R = mp.mpf('2.5'); M = int(sys.argv[1]) if len(sys.argv) > 1 else 96
vals = [H0(c + R*mp.expjpi(2*mp.mpf(j)/M)) for j in range(M)]
print('samples', M, 'time', round(time.time()-T0, 1), flush=True)
K = M//2
a = []
for k in range(K):
    a.append(sum(vals[j]*mp.expjpi(-2*mp.mpf(j*k)/M) for j in range(M))/(M*R**k))
scale = a[2]; a = [x/scale for x in a]          # normalise: a_2 = 1
print('max |Im a_k|/|a_k| (H real on R => a_k real):', mp.nstr(max(abs(mp.im(x))/(abs(x)+mp.mpf(10)**-60) for x in a[:30]), 3))
a = [mp.re(x) for x in a]
print('tail |a_k| R^k for k = K-6..K-1:', [mp.nstr(abs(a[k])*R**k, 3) for k in range(K-6, K)])
# check series vs direct at an interior point
wt = mp.mpc('0.3', '0.2'); ser = sum(a[k]*wt**k for k in range(K)); dirv = H0(c + wt)/scale
print('series vs direct at w=0.3+0.2i: rel err', mp.nstr(abs(ser-dirv)/abs(dirv), 3), flush=True)
fact = [mp.factorial(n) for n in range(K+1)]
def coeffs(t):
    t = mp.mpf(t)
    return [sum(a[m+2*j]*(-t)**j*fact[m+2*j]/(fact[j]*fact[m]) for j in range((K-1-m)//2+1)) for m in range(K)]
def ev(b, w): return mp.polyval(b[::-1], w)
def der(b): return [b[m+1]*(m+1) for m in range(len(b)-1)]
json.dump({'c': str(c), 'R': str(R), 'M': M, 'scale': str(scale), 'a': [str(x) for x in a]}, open(LOG + 'taylor_lehmer.json', 'w'))
# zeros of H_0 near 0 vs zetazero
b0 = coeffs(0)
w1 = mp.findroot(lambda w: ev(b0, w), -(g[6710]-g[6709])); w2 = mp.findroot(lambda w: ev(b0, w), (g[6710]-g[6709]))
print('t=0 pair from series: w1, w2 =', mp.nstr(w1, 15), mp.nstr(w2, 15), ' zetazero:', mp.nstr(2*g[6709]-c, 15), mp.nstr(2*g[6710]-c, 15))
u0 = ((w2-w1)/2)**2
# tidal constant: S_all = -(log H)''(0) = (a1/a0)^2 - 2 a2/a0 ; pair part 1/w1^2 + 1/w2^2
S_all = (a[1]/a[0])**2 - 2*a[2]/a[0]; S_pair = 1/w1**2 + 1/w2**2; Sfar = S_all - S_pair
# direct sum over the listed neighbours (6700..6720 minus the pair), both signs +-z_j
Sdir = sum(1/(c - 2*g[n])**2 + 1/(c + 2*g[n])**2 for n in g if n not in (6709, 6710))
print('u0 =', mp.nstr(u0, 12), ' S_far =', mp.nstr(Sfar, 10), ' (listed 19 neighbours alone:', mp.nstr(Sdir, 8), ')')
print('Wazewski/comparison prediction: du/dt = 2 - 4 S_far u  =>  t_c = -(1/(4 S_far)) log(1/(1 - 2 S_far u0)) =',
      mp.nstr(-mp.log(1/(1-2*Sfar*u0))/(4*Sfar), 12), ';  bare -u0/2 =', mp.nstr(-u0/2, 12))
# exact collision: saddle w*(t) (zero of H_t'), critical value v(t) = H_t(w*); v(t_c) = 0
def crit(t):
    b = coeffs(t); bd = der(b)
    ws = mp.findroot(lambda w: ev(bd, w), mp.mpf(0))
    return ev(b, ws), ws, b
tc = mp.findroot(lambda t: crit(t)[0], -u0/2, tol=mp.mpf(10)**-30)
vc, wsc, bc = crit(tc)
print('EXACT collision time t_c =', mp.nstr(tc, 15), ' saddle at w* =', mp.nstr(wsc, 12), ' |v| =', mp.nstr(abs(vc), 3), flush=True)
def pair(t):
    b = coeffs(t); d = mp.sqrt(abs(u0 + 2*(t - 0)))
    if t > tc: guess = (wsc - d, wsc + d)
    else: guess = (wsc - 1j*mp.sqrt(abs(2*(t-tc))), wsc + 1j*mp.sqrt(abs(2*(t-tc))))
    z1 = mp.findroot(lambda w: ev(b, w), guess[0]); z2 = mp.findroot(lambda w: ev(b, w), guess[1])
    return z1, z2, ((z2-z1)/2)**2
for dl in ['1e-6', '1e-5', '1e-4']:
    dl = mp.mpf(dl)
    up = pair(tc + dl)[2]; um = pair(tc - dl)[2]
    print('du/dt at t_c (central diff, delta=%s): %s   [u(t_c+d) = %s, u(t_c-d) = %s]' % (mp.nstr(dl, 2), mp.nstr((up-um)/(2*dl), 12), mp.nstr(up, 6), mp.nstr(um, 6)))
res = {'tc': str(tc), 'u0': str(u0), 'Sfar': str(Sfar), 'wstar_tc': str(wsc)}
json.dump(res, open(LOG + 'lehmer_tc.json', 'w'))
print('total time', round(time.time()-T0, 1))
