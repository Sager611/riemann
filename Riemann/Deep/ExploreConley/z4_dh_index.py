# z4_dh_index.py — explore66: Conley data of DH's off-line pair across its landing t* (DH's own heat flow), block |w| <= r0.
import mpmath as mp, json, sys
mp.mp.dps = 40
LOG = '/Users/adriansager/playground/riemann/Riemann/Deep/ExploreConley/logs/'
d = json.load(open(LOG + 'taylor_dh.json')); a = [mp.mpf(x) for x in d['a']]; K = len(a)
ts = mp.mpf(json.load(open(LOG + 'dh_tstar.json'))['tstar'])
fact = [mp.factorial(n) for n in range(K+1)]
def coeffs(t):
    t = mp.mpf(t)
    return [sum(a[m+2*j]*(-t)**j*fact[m+2*j]/(fact[j]*fact[m]) for j in range((K-1-m)//2+1)) for m in range(K)]
def ev(b, w):
    s = mp.mpc(0)
    for x in reversed(b): s = s*w + x
    return s
def der(b): return [b[m+1]*(m+1) for m in range(len(b)-1)]
r0 = mp.mpf(sys.argv[1]) if len(sys.argv) > 1 else mp.mpf('1.2'); N = int(sys.argv[3]) if len(sys.argv) > 3 else 1440; tmax = mp.mpf('0.2'); nst = int(sys.argv[2]) if len(sys.argv) > 2 else 40
ws = [r0*mp.expjpi(2*mp.mpf(k)/N) for k in range(N)]
bt = [sum(abs(a[m+2*j])*tmax**j*fact[m+2*j]/(fact[j]*fact[m]) for j in range((K-1-m)//2+1)) for m in range(K)]
B0 = sum(bt[m]*r0**m for m in range(K)); B1 = sum(m*bt[m]*r0**(m-1) for m in range(1, K))
B2 = sum(m*(m-1)*bt[m]*r0**(m-2) for m in range(2, K)); B3 = sum(m*(m-1)*(m-2)*bt[m]*r0**(m-3) for m in range(3, K))
print('r0 =', r0, ' majorants (|t| <= 0.2): B0..B3 =', [mp.nstr(x, 5) for x in (B0, B1, B2, B3)])
arc = 2*mp.pi*r0/N
zoff = mp.mpc(json.load(open(LOG + 'dh_tstar.json'))['zoff']); cc = mp.mpf(json.load(open(LOG + 'dh_tstar.json'))['c'])
y0 = mp.im(zoff); w0 = zoff - cc
S_all = (a[1]/a[0])**2 - 2*a[2]/a[0]; S_pair = 1/w0**2 + 1/mp.conj(w0)**2; Sfar = S_all - S_pair
print('DH pair w0 =', mp.nstr(w0, 10), ' S_far =', mp.nstr(Sfar, 10), ' frozen-tidal comparison: t* = log(1 + 2 S y0^2)/(4 S) =', mp.nstr(mp.log(1 + 2*Sfar*y0**2)/(4*Sfar), 10), ' exact t* =', mp.nstr(ts, 12), ' de Bruijn y0^2/2 =', mp.nstr(y0**2/2, 10))
def wind(vals):
    return int(mp.nint(sum(mp.arg(vals[(k+1) % N]/vals[k]) for k in range(N))/(2*mp.pi)))
grid = sorted(set([tmax*k/nst for k in range(nst+1)] + [ts]))
minH = minHp = minrad = mp.inf; minprod = mp.inf; prev = None
for t in grid:
    b = coeffs(t); bd = der(b)
    Hv = [ev(b, w) for w in ws]; Hpv = [ev(bd, w) for w in ws]
    mH = min(abs(v) for v in Hv); mHp = min(abs(v) for v in Hpv)
    radv = [mp.re(mp.conj(Hv[k])*Hpv[k]*ws[k]/r0) for k in range(N)]
    rad = min(radv[k]/(abs(Hv[k])*abs(Hpv[k])) for k in range(N)); minprod = min(minprod, min(radv))
    nH, nHp = wind(Hv), wind(Hpv)
    wsd = mp.findroot(lambda w: ev(bd, w), mp.mpf('-0.1'))
    if t > ts: dd = mp.sqrt(2*(t-ts)); g1, g2 = wsd-dd, wsd+dd
    elif t < ts: dd = mp.sqrt(2*(ts-t)); g1, g2 = wsd-1j*dd, wsd+1j*dd
    else: g1 = g2 = wsd
    if t != ts:
        z1 = mp.findroot(lambda w: ev(b, w), g1); z2 = mp.findroot(lambda w: ev(b, w), g2)
        kind = 'real pair' if abs(mp.im(z1)) < mp.mpf(10)**-25 else 'conjugate pair'
    else:
        z1 = z2 = wsd; kind = 'double'
    prev = (z1, z2, kind)
    s1, s2 = mp.sign(mp.re(ev(b, -r0))), mp.sign(mp.re(ev(b, r0)))
    minH = min(minH, mH); minHp = min(minHp, mHp); minrad = min(minrad, rad)
    print('t=%-12s zeros(H)=%d zeros(H\')=%d %-14s z=%s, %s saddle=%s min|H|=%s min|H\'|=%s min cos=%s signs %d,%d' % (
        mp.nstr(t, 9), nH, nHp, kind, mp.nstr(z1, 7), mp.nstr(z2, 7), mp.nstr(mp.re(wsd), 6), mp.nstr(mH, 4), mp.nstr(mHp, 4), mp.nstr(rad, 4), s1, s2), flush=True)
Delta = tmax/nst
slH = B2*Delta/2 + B1*arc/2; slHp = B3*Delta/2 + B2*arc/2
slR = (B2*B1 + B0*B3)*Delta/2 + (B1*B1 + B0*B2)*arc/2
print('CERT no zero of H_t on dD, all t in [0,0.2]:', minH > slH, mp.nstr(minH, 5), '>', mp.nstr(slH, 5))
print("CERT no zero of H_t' on dD, all t in [0,0.2]:", minHp > slHp, mp.nstr(minHp, 5), '>', mp.nstr(slHp, 5))
print('CERT exit set empty (Re conj(H) H_r > 0 on dD), all t in [0,0.2]:', minprod > slR, mp.nstr(minprod, 5), '>', mp.nstr(slR, 5))
