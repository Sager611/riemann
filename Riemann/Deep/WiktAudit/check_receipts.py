# Independent recomputation (mpmath only, no flint, no author code) of the a-dependent analytic constants of the
# Wiktorowicz n100 certificate and of the final budget from the deposited receipts.
import json, sys
from mpmath import mp, mpf, mpc, log, sqrt, pi, sin, cos, sinh, cosh, exp, digamma, fabs, floor, ceil
mp.dps = 60
R = sys.argv[1]   # results dir with a23029_{even,odd}_pipeline.json, _ldl_1.json
a = mpf(23029)/10000
ppow = []
for m in range(2, 101):
    p = next(d for d in range(2, m+1) if m % d == 0); x = m
    while x % p == 0: x //= p
    if x == 1: ppow.append((m, p))
assert len(ppow) == 35 and exp(2*a) > 100 and exp(2*a) < 101
w = {m: log(p)/sqrt(m) for m, p in ppow}
phi = {m: pi*log(m)/a for m, _ in ppow}
RHO_W = mpf(11101037276778)/10**12
rhoH = a/2*RHO_W
def dfree(k): return a/2*(digamma(mpc(mpf(1)/4, pi*k/(2*a))).real - log(pi))
def phase_consts():
    V = sum(w[m]**2 for m in w)/2
    E1 = sum(w[m]/fabs(sin(phi[m]/2)) for m in w)
    E = mpf(0)
    for m in w:
        for n in w:
            if m != n: E += w[m]*w[n]/(2*fabs(sin((phi[m]-phi[n])/2)))
            E += w[m]*w[n]/(2*fabs(sin((phi[m]+phi[n])/2)))
    c = a/4; b = a/pi
    return c*c + b*b*V, 2*c*b*E1 + b*b*E, c + b*sum(w.values())
Vh, Eh, Hc = phase_consts()
print('V_h', mp.nstr(Vh, 15), 'E_h', mp.nstr(Eh, 15), 'H', mp.nstr(Hc, 15))
out = {}
for par, name in ((0, 'even'), (1, 'odd')):
    pj = json.load(open(f'{R}/a23029_{name}_pipeline.json'))
    lj = json.load(open(f'{R}/a23029_{name}_ldl_1.json'))
    far = pj['far']; asm = pj['assembly']
    def ball_mid(s): return mpf(s.strip('[]').split('+/-')[0])
    N, M, K = pj['N'], pj['M'], pj['K']
    om = mpf(par)/2
    kN = N + om; kM = M + om
    J2 = lambda kk: Vh/(kk-1) + Eh/kk**2
    J4 = lambda kk: 1/(3*(kk-1)**3)
    alpha = sinh(a/2)**2 if par == 0 else -cosh(a/2)**2
    rem = a**2/pi**2*(1 + exp(-a)/(1-exp(-4*a)) + fabs(alpha)/(mpf(1)/4 + (pi*kN/a)**2))
    pole = 0 if par == 0 else (a/pi)**4*fabs(alpha)/(3*(kN-1)**3)
    cov = ball_mid(far['covariance'])
    assert cov + pole < mpf('0.227')
    reserve = mpf('0.227')
    tN = dfree(kN) - rhoH - reserve; tM = dfree(kM) - rhoH - reserve
    aM, bM, vM = ball_mid(far['aM']), ball_mid(far['bM']), ball_mid(far['vM'])
    SU, SV = ball_mid(asm['SU']), ball_mid(asm['SV'])
    J2M, J4M = J2(kM), J4(kM)
    kN1 = N - 1 + om
    omitted = (sqrt(SU)*vM*sqrt(J4M) + sqrt(SV)*(aM*sqrt(J2M) + (aM+bM)*sqrt(J4M)))**2/(tM*(1-(kN1/kM)**2)**2)
    # (7.5) assembly error from componentwise errors and magnitudes
    errs = [ball_mid(s) for s in asm['errors']]; mags = [ball_mid(s) for s in asm['magnitudes']]
    Hs = (1 + log(N-K))/K
    ae = errs[0] + 2*Hs*sum(errs[iu]*mags[iv] + errs[iv]*mags[iu] + errs[iu]*errs[iv] for iu, iv in ((1,2),(3,4)))
    # (8.5) LDL budget, re-implemented from the paper (8.2)-(8.5) with the receipt maxima
    def up(s): x = mpf(s); return x + (abs(x)+1)*mpf('1e-22')
    pmin = mpf(lj['pmin']) - (abs(mpf(lj['pmin']))+1)*mpf('1e-22')
    p0 = (floor(pmin*1000)-1)/1000
    C = ceil(up(lj['cmax'])) + 1; D = ceil(up(lj['dmax'])) + 1
    G4 = [ceil(up(s)) + 1 for s in lj['generator_max']]
    G = 2*(G4[0]*G4[2] + G4[1]*G4[3])
    u = mpf(2)**(-(lj['mantissa_bits']-1))
    n = N - K
    resid = n*(64*u*G*(1+C/p0)*Hs + 4*u*(D + C*C/p0)) + mpf('1e-12')
    shift = mpf('0.30')
    margin = shift - resid - ae - omitted
    rec = {
      'r_bound': (rem, ball_mid(far['r_bound'])), 'pole': (pole, ball_mid(far['pole'])),
      'tN': (tN, ball_mid(far['tN'])), 'tM': (tM, ball_mid(far['tM'])),
      'J2M': (J2M, ball_mid(far['J2M'])), 'J4M': (J4M, ball_mid(far['J4M'])),
      'omitted': (omitted, ball_mid(asm['omitted'])), 'assembly_error': (ae, ball_mid(asm['assembly_error'])),
      'residual(%s)' % lj['precision']: (resid, None), 'margin': (margin, None)}
    print('==', name, 'N', N, 'M', M, 'K', K, 'precision', lj['precision'], 'pivots', lj['pivots'], 'completed', lj['completed'])
    for k2, (ours, theirs) in rec.items():
        if theirs is None: print('  %-16s ours %s' % (k2, mp.nstr(ours, 20)))
        else:
            rel = fabs(ours-theirs)/max(fabs(theirs), mpf('1e-300')) if theirs != 0 else fabs(ours)
            print('  %-16s ours %s  theirs %s  rel.diff %s' % (k2, mp.nstr(ours, 20), mp.nstr(theirs, 20), mp.nstr(rel, 3)))
    print('  p0 %s C %s D %s G %s u 2^-%d' % (p0, C, D, G4, lj['mantissa_bits']-1))
