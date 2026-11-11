# T1 check: graphical Krein signature = phase velocity of the Riemann-Siegel half E(t)=e^{i theta} zeta_N(1/2+it)
import mpmath as mp
mp.mp.dps = 25
N = 33
def phidot(t):
    s = mp.mpf(1)/2 + 1j*t
    zN = mp.fsum(mp.power(n, -s) for n in range(1, N+1))
    dzN = mp.fsum(-mp.log(n)*mp.power(n, -s) for n in range(1, N+1))
    return mp.diff(mp.siegeltheta, t) + mp.re(dzN/zN)
ts = [mp.mpf(7004.90) + k*mp.mpf('0.01') for k in range(41)]
vals = [(t, phidot(t)) for t in ts]
neg = [t for t, v in vals if v < 0]
print('phase velocity of E on [7004.90, 7005.30]: min %s at t=%s; negative on %s..%s (%d of 41 grid points)' % (
  mp.nstr(min(v for t, v in vals), 5), mp.nstr(min(vals, key=lambda r: r[1])[0], 8),
  mp.nstr(neg[0], 7) if neg else '-', mp.nstr(neg[-1], 7) if neg else '-', len(neg)))
print("theta'(t) =", mp.nstr(mp.diff(mp.siegeltheta, mp.mpf(7005.08)), 6))
for t in (mp.mpf('7008.11474206'), mp.mpf('7005.0628661749'), mp.mpf('7005.1005646726')):
    print('phidot(%s) = %s' % (mp.nstr(t, 12), mp.nstr(phidot(t), 5)))
