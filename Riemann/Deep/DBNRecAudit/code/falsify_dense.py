#!/usr/bin/env python3
# dbnrecaudit falsify_dense.py: independent float/mpmath re-implementation of the actual signed head Q_i(t,x) of LEMMA.md and of the
# shell-0 bound function Q+ (own formulas; inputs read from the Zenodo plan/affine files and the logged certs only), then dense sampling
# (C helper 'dense', 10^6 points per window) at each row's binding cell: tests (L2)  Q_i - s log(x/X) <= (1-theta)A' + theta B',
# Q_i(t,x) <= Q+(x') + (1-theta) D, max_cell Q+ <= certify's Arb UB, and the termwise amplitude bounds. Falsification only: nothing certified.
# usage: falsify_dense.py ROW [ROW ...]   (cwd = repo root; writes Riemann/Deep/DBNRecAudit/logs/dense_r<ROW>.json)
import sys, json, subprocess, os
import mpmath as mp
from fractions import Fraction as Fr
mp.mp.dps = 34
DL = '.scratch/dbnrecaudit/dl/Proof_and_Replay__Lambda_Peer_Review__'
OUT = 'Riemann/Deep/DBNRecAudit/logs/'; TMP = '/tmp/dbnra_keep/'; DENSE = TMP + 'dense'
plan = {int(r['row']): r for r in json.load(open(DL + 'inputs__workspace__lambda_round14__material__m8_late__M8_LATE_PLAN.json'))['rows']}
aff = {int(r['row']): r for r in json.load(open(DL + 'verification__referee__ALL_67_ORIGINAL_AFFINE_NATIVE.json'))['rows']}
K = 4096; X = mp.mpf(6000000185827); PI = mp.pi; Dd = 4*3.141592653589793/0.6931471805599453   # bb anchor spacing (the same double)
import math; assert Dd == 4*math.pi/math.log(2)
ln = [mp.mpf(0), mp.mpf(0)] + [mp.log(n) for n in range(2, K + 1)]
def dcoef(t):   # own implementation of U.tex l.3311-3316: d_t = (b_t log) * b_t^{-1}
    b = [mp.mpf(0), mp.mpf(1)] + [mp.exp(t*ln[n]**2/4) for n in range(2, K + 1)]
    d = [mp.mpf(0)]*(K + 1)
    for n in range(2, K + 1):
        d[n] = ln[n]*b[n]
    for n in range(2, K + 1):          # d[n] final once all divisors r < n are final (increasing n)
        r = 2
        while r*r <= n:
            if n % r == 0:
                q = n//r
                d[n] -= d[r]*b[q]
                if q != r: d[n] -= d[q]*b[r]
            r += 1
    return d
def isprime(n): return n > 1 and all(n % q for q in range(2, int(n**0.5) + 1))
def run(f, mode, a, b, npts):
    out = subprocess.run([DENSE, f, mode, repr(a), repr(b), str(npts)], capture_output=True, text=True, check=True).stdout.split()
    return dict(max=float(out[0]), argmax=float(out[1]), min=float(out[2]), mean=float(out[3]))
for row in map(int, sys.argv[1:]):
    p = plan[row]; ar = aff[row]
    assert p.get('joint_jet_overlay') is None
    tm, tp = (mp.mpf(Fr(v).numerator)/Fr(v).denominator for v in p['time']); tmq, tpq = map(Fr, p['time'])
    nodes = [(Fr(n['height']), Fr(n['shift']), Fr(n['weight'])) for n in p['nodes']]
    M = sum(w for _, _, w in nodes); beta = Fr(ar['beta']); assert M == Fr(ar['mass'])
    s = (1 - beta)*M/4; sf = mp.mpf(s.numerator)/s.denominator
    cdir = 'Riemann/Deep/DBNRec/logs/' if row <= 3 else 'Riemann/Deep/DBNRec2/logs/'
    cert = json.load(open(cdir + 'cert_r%d.json' % row)); Ap, Bp = Fr(cert['A_new']), Fr(cert['B_new'])
    sh0 = cert['shells'][0]; w = sh0['worst']['cell'].split(); ub_arb = sh0['worst']['ub_arb']
    xs = X + PI*tm/4                                                    # shell 0 (u0 = 0) scan start, prep_end.py l.87
    if len(w) == 4: xc = xs + mp.mpf(int(w[0]))*mp.mpf(Dd) + mp.mpf(w[1]); half = float(w[2])/2; slack = 1e-9
    else: xc = xs + mp.mpf(w[0]); half = Dd/2048/32; slack = 1e-4          # row 1: old dump (location error < 1e-4)
    F = [mp.mpf(0)]*(K + 1)
    for n in range(2, K + 1): F[n] = sum(mp.mpf(w_.numerator)/w_.denominator*mp.exp(-(mp.mpf(h.numerator)/h.denominator)*ln[n]/2)*mp.cos((mp.mpf(a.numerator)/a.denominator)*ln[n]/2) for h, a, w_ in nodes)
    dm, dp = dcoef(tm), dcoef(tp)
    L0 = mp.log(X/(4*PI)); Llo, Lhi = L0, L0 + mp.mpf(1)/100         # shell 0 = [X, X e^{1/100}]
    phase = [0.0]*(K + 1)
    for n in range(2, K + 1): phase[n] = float(mp.fmod(xc*ln[n]/2, 2*PI))
    Hp, Hm, lo_ = {}, {}, {}
    for n in range(2, K + 1):
        aF = abs(F[n])
        if isprime(n):
            base = ln[n]*n**mp.mpf(-0.5)/2
            Hp[n] = base*mp.exp(tp*ln[n]*(ln[n] - Llo)/4)*aF; Hm[n] = base*mp.exp(tm*ln[n]*(ln[n] - Llo)/4)*aF; lo_[n] = base*mp.exp(tp*ln[n]*(ln[n] - Lhi)/4)*aF
        else:
            Hp[n] = dp[n]*mp.exp(-(mp.mpf(1)/2 + tp*Llo/4)*ln[n])*aF/2; Hm[n] = dp[n]*mp.exp(-(mp.mpf(1)/2 + tm*Llo/4)*ln[n])*aF/2
            lo_[n] = dm[n]*mp.exp(-(mp.mpf(1)/2 + tp*Lhi/4)*ln[n])*aF/2
    Dch = sum(Hm[n] - Hp[n] for n in range(2, K + 1))
    res = dict(row=row, t=[str(tmq), str(tpq)], s=float(sf), A_new=float(Ap), B_new=float(Bp), xc_minus_xs=float(xc - xs), u=float(mp.log((xc - PI*tp/4)/X)),
               D_shell0_own=float(Dch), D_shell0_cert=sh0['D'], ub_arb=ub_arb, runs=[])
    # termwise amplitude checks at the binding cell for several t (actual x = x' - pi t/4)
    worst_up, worst_lo = -1e9, -1e9
    tl = [tm + (tp - tm)*mp.mpf(j)/8 for j in range(9)]
    dts = {0: dm, 8: dp}
    for j, t in enumerate(tl):
        dt = dts.get(j) or dcoef(t); th = (t - tm)/(tp - tm)
        x = xc - PI*t/4; L = mp.log(x/(4*PI))
        for n in range(2, K + 1):
            m = dt[n]*mp.exp(-(mp.mpf(1)/2 + t*L/4)*ln[n])*abs(F[n])/2
            ch = (1 - th)*Hm[n] + th*Hp[n]
            if ch > 0: worst_up = max(worst_up, float((m - ch)/ch))
            if m > 0: worst_lo = max(worst_lo, float((lo_[n] - m)/m))
        if j in ((0, 4, 8) if row == 1 else (0, 8)):
            f = TMP + 'a_r%d_t%d.txt' % (row, j)
            with open(f, 'w') as fo:
                fo.write('%d\n' % (K - 1))
                for n in range(2, K + 1):
                    amp = dt[n]*mp.exp(-(mp.mpf(1)/2 + t*L/4)*ln[n])*F[n]/2
                    fo.write('%.17g %.17g 0 %.17g\n' % (float(ln[n]/2), float(amp), phase[n]))
            rhs = float((1 - th)*(mp.mpf(Ap.numerator)/Ap.denominator) + th*(mp.mpf(Bp.numerator)/Bp.denominator))
            slope = float(sf*mp.log((xc - 2 - PI*t/4)/X))   # smallest x over the wide window (conservative for falsification)
            for name, a_, b_, npts in ((('cell', -half - slack, half + slack, 1000000), ('wide', -2.0, 2.0, 200000)) if (row in (1, 8) and j in (0, 8)) else (('cell', -half - slack, half + slack, 1000000),)):
                r_ = run(f, 'a', a_, b_, npts)
                res['runs'].append(dict(kind='actual', theta=float(th), window=name, Qmax=r_['max'], at=r_['argmax'], G=r_['max'] - slope, rhs=rhs, margin=rhs - (r_['max'] - slope)))
            os.remove(f)
    res['amp_upper_worst_rel'] = worst_up; res['amp_lower_worst_rel'] = worst_lo   # must be <= 0 (up to 1e-30)
    f = TMP + 'b_r%d.txt' % row
    with open(f, 'w') as fo:
        fo.write('%d\n' % (K - 1))
        for n in range(2, K + 1):
            if F[n] > 0: lo, hi = lo_[n], Hp[n]
            else: lo, hi = -Hp[n], -lo_[n]
            fo.write('%.17g %.17g %.17g %.17g\n' % (float(ln[n]/2), float(lo), float(hi), phase[n]))
    r_ = run(f, 'b', -half - slack, half + slack, 1000000)
    res['Qplus_cell'] = dict(max=r_['max'], at=r_['argmax'], ub_arb=ub_arb, slack=ub_arb - r_['max'])
    r_ = run(f, 'b', -2.0, 2.0, 200000); res['Qplus_wide'] = dict(max=r_['max'], at=r_['argmax'], npts=200000)
    os.remove(f)
    json.dump(res, open(OUT + 'dense_r%d.json' % row, 'w'), indent=1)
    print(json.dumps(dict(row=row, minmargin=min(r['margin'] for r in res['runs']), Qplus_cell=res['Qplus_cell'], amp=[worst_up, worst_lo], D=[res['D_shell0_own'], res['D_shell0_cert']])), flush=True)

