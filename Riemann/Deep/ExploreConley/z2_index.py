# z2_index.py — explore66: Conley data of the Lehmer pair across its collision t_c (zeta's own heat flow).
# Block D = {|w| <= r0} around c (midpoint of zeros 6709, 6710), t in [-2e-3, 0].
# (i) 2D gradient-descent flow of |H_t|^2 (= Newton flow up to time change): exit set E_t = {w in dD : d_r |H_t|^2 < 0};
#     E_t empty  =>  D is an attracting isolating block, Conley index = pointed S^0 (D/empty = D u {*}).
# (ii) zeros of H_t (minima) and of H_t' (saddles) inside D: counts by winding number, positions, real or conjugate pair.
# (iii) 1D flow xdot = H_t(x) on [-r0, r0]: index from the signs of H_t(+-r0).
# Certification (40-digit floating, Taylor tail < 1e-28 relative): no zero of H_t, H_t' on dD for ANY t in [-2e-3, 0]
#   via min_dD |H_{t_i}| > (Delta/2) * B2 + sampling slack, B2 >= max |d_t H| = max |H''| over the whole slab (majorant series).
import mpmath as mp, json, time
mp.mp.dps = 40
LOG = '/Users/adriansager/playground/riemann/Riemann/Deep/ExploreConley/logs/'
d = json.load(open(LOG + 'taylor_lehmer.json')); a = [mp.mpf(x) for x in d['a']]; K = len(a)
tc = mp.mpf(json.load(open(LOG + 'lehmer_tc.json'))['tc'])
fact = [mp.factorial(n) for n in range(K+1)]
def coeffs(t):
    t = mp.mpf(t)
    return [sum(a[m+2*j]*(-t)**j*fact[m+2*j]/(fact[j]*fact[m]) for j in range((K-1-m)//2+1)) for m in range(K)]
def ev(b, w):
    s = mp.mpc(0)
    for x in reversed(b): s = s*w + x
    return s
def der(b): return [b[m+1]*(m+1) for m in range(len(b)-1)]
r0 = mp.mpf('0.6'); N = 1440; smax = mp.mpf('0.0021')
ws = [r0*mp.expjpi(2*mp.mpf(k)/N) for k in range(N)]
# majorants valid for every t with |t| <= smax and |w| <= r0
bt = [sum(abs(a[m+2*j])*smax**j*fact[m+2*j]/(fact[j]*fact[m]) for j in range((K-1-m)//2+1)) for m in range(K)]
B0 = sum(bt[m]*r0**m for m in range(K))                     # >= |H_t|
B1 = sum(m*bt[m]*r0**(m-1) for m in range(1, K))            # >= |H_t'|
B2 = sum(m*(m-1)*bt[m]*r0**(m-2) for m in range(2, K))      # >= |H_t''| = |d_t H_t|
B3 = sum(m*(m-1)*(m-2)*bt[m]*r0**(m-3) for m in range(3, K))  # >= |H_t'''| = |d_t H_t'|
print('majorants on the slab |w|<=0.6, |t|<=2.1e-3:  B1 =', mp.nstr(B1, 6), ' B2 =', mp.nstr(B2, 6), ' B3 =', mp.nstr(B3, 6))
arc = 2*mp.pi*r0/N
def wind(vals):
    tot = mp.mpf(0)
    for k in range(N):
        tot += mp.arg(vals[(k+1) % N]/vals[k])
    return int(mp.nint(tot/(2*mp.pi)))
grid = [mp.mpf(-k)/10000 for k in range(0, 21)] + [tc]
grid = sorted(set(grid), reverse=True)
rows = []; minH = mp.inf; minHp = mp.inf; minrad = mp.inf; minprod = mp.inf
for t in grid:
    b = coeffs(t); bd = der(b)
    Hv = [ev(b, w) for w in ws]; Hpv = [ev(bd, w) for w in ws]
    mH = min(abs(v) for v in Hv); mHp = min(abs(v) for v in Hpv)
    # radial derivative of |H|^2: 2 Re(conj(H) H' e^{i theta}) ; normalise by |H| |H'|
    radv = [mp.re(mp.conj(Hv[k])*Hpv[k]*ws[k]/r0) for k in range(N)]; minprod = min(minprod, min(radv))
    rad = min(radv[k]/(abs(Hv[k])*abs(Hpv[k])) for k in range(N))
    nH, nHp = wind(Hv), wind(Hpv)
    # positions: saddle (zero of H') near 0, and the two zeros
    wsd = mp.findroot(lambda w: ev(bd, w), mp.mpf(0))
    if t >= tc:
        dd = mp.sqrt(max(2*(t-tc), mp.mpf(10)**-30)); g1, g2 = wsd-dd, wsd+dd
    else:
        dd = mp.sqrt(2*(tc-t)); g1, g2 = wsd-1j*dd, wsd+1j*dd
    try:
        z1 = mp.findroot(lambda w: ev(b, w), g1); z2 = mp.findroot(lambda w: ev(b, w), g2)
    except Exception as e:
        z1 = z2 = mp.mpc(wsd)
    kind = 'double' if t == tc else ('real pair' if abs(mp.im(z1)) < mp.mpf(10)**-25 else 'conjugate pair')
    s1, s2 = mp.sign(mp.re(ev(b, -r0))), mp.sign(mp.re(ev(b, r0)))
    minH = min(minH, mH); minHp = min(minHp, mHp); minrad = min(minrad, rad)
    rows.append((t, nH, nHp, kind, z1, z2, wsd, mH, mHp, rad, s1, s2))
    print('t=%-22s zeros(H)=%d zeros(H\')=%d  %-14s z=%s, %s  saddle=%s  min|H|=%s min|H\'|=%s  min cos(rad)=%s  sign H(-r0),H(r0) = %d,%d' % (
        mp.nstr(t, 10), nH, nHp, kind, mp.nstr(z1, 8), mp.nstr(z2, 8), mp.nstr(mp.re(wsd), 6), mp.nstr(mH, 5), mp.nstr(mHp, 5), mp.nstr(rad, 4), s1, s2), flush=True)
Delta = mp.mpf('1e-4')
slackH = B2*Delta/2 + B1*arc/2; slackHp = B3*Delta/2 + B2*arc/2
print('CERT no zero of H_t on dD for all t in [-2e-3,0]:', minH > slackH, ' (min %s > slack %s)' % (mp.nstr(minH, 5), mp.nstr(slackH, 5)))
print('CERT no zero of H_t\' on dD for all t in [-2e-3,0]:', minHp > slackHp, ' (min %s > slack %s)' % (mp.nstr(minHp, 5), mp.nstr(slackHp, 5)))
slR = (B2*B1 + B0*B3)*Delta/2 + (B1*B1 + B0*B2)*arc/2
print('CERT exit set empty (Re conj(H) H_r > 0 on dD), all t in [-2e-3,0]:', minprod > slR, mp.nstr(minprod, 5), '>', mp.nstr(slR, 5), ' (min normalised cos %s)' % mp.nstr(minrad, 5))
json.dump([[str(r[0]), r[1], r[2], r[3], str(r[4]), str(r[5]), str(r[6]), str(r[7]), str(r[8]), str(r[9]), int(r[10]), int(r[11])] for r in rows], open(LOG + 'z2_index.json', 'w'))
