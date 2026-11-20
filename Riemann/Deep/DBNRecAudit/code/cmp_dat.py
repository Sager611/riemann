#!/usr/bin/env python3
# dbnrecaudit cmp_dat.py: rounding directions of prep_end.py (TABLE B4): the doubles lo/hi of a .dat (row ROW, shell [u0, u1]) against own
# mpmath values (34 digits) of the true interval s_n[l_n, H+_n]; also D (chord difference), tri+/tri-, xs, K, the uint64 phase data p0/dK/d1/d2/d3
# and the exact phase drift at the last anchor (needs <= EPSU = 2^40 units for levels 1-3, <= 2^34 units = 2^-30 turn for PHERR).
# usage: cmp_dat.py DAT PREPJSON ROW   (cwd = repo root; prints JSON)
import sys, json, math
sys.argv, args = sys.argv[:1], sys.argv[1:]
sys.path.insert(0, 'Riemann/Deep/DBNRecAudit/code')
import falsify_dense as fd
mp = fd.mp
from fractions import Fraction as Fr
dat, pjf, row = args[0], args[1], int(args[2])
pj = json.loads(open(pjf).readline())
p = fd.plan[row]; K_ = fd.K; ln = fd.ln; X = fd.X; PI = fd.PI
q = lambda r: mp.mpf(r.numerator)/r.denominator
tm, tp = (q(Fr(v)) for v in p['time']); u0, u1 = q(Fr(pj['u0'])), q(Fr(pj['u1']))
nodes = [(Fr(n['height']), Fr(n['shift']), Fr(n['weight'])) for n in p['nodes']]
F = [mp.mpf(0)]*(K_ + 1)
for n in range(2, K_ + 1): F[n] = sum(q(w)*mp.exp(-q(h)*ln[n]/2)*mp.cos(q(a)*ln[n]/2) for h, a, w in nodes)
dm, dp = fd.dcoef(tm), fd.dcoef(tp)
L0 = mp.log(X/(4*PI)); Llo, Lhi = L0 + u0, L0 + u1
tlo, thi, Dsum, tri_p, tri_m = {}, {}, mp.mpf(0), mp.mpf(0), mp.mpf(0)
for n in range(2, K_ + 1):
    aF = abs(F[n])
    if fd.isprime(n):
        base = ln[n]*n**mp.mpf(-0.5)/2; Hp = base*mp.exp(tp*ln[n]*(ln[n] - Llo)/4)*aF; Hm = base*mp.exp(tm*ln[n]*(ln[n] - Llo)/4)*aF
        l = base*mp.exp(tp*ln[n]*(ln[n] - Lhi)/4)*aF
    else:
        Hp = dp[n]*mp.exp(-(mp.mpf(1)/2 + tp*Llo/4)*ln[n])*aF/2; Hm = dp[n]*mp.exp(-(mp.mpf(1)/2 + tm*Llo/4)*ln[n])*aF/2
        l = dm[n]*mp.exp(-(mp.mpf(1)/2 + tp*Lhi/4)*ln[n])*aF/2
    Dsum += Hm - Hp; tri_p += Hp; tri_m += Hm
    tlo[n], thi[n] = (l, Hp) if F[n] > 0 else (-Hp, -l)
f = open(dat); hdr = f.readline().split(); NG, Kan = int(hdr[0]), int(hdr[1]); D = float(hdr[2])
two64 = mp.mpf(2)**64
xs = X*mp.exp(u0) + PI*tm/4; xe = X*mp.exp(u1) + PI*tp/4
worst_hi, worst_lo, nn, ph_err, drift = 1e300, 1e300, 0, 0, 0
M1, M2, M3 = int(hdr[3]), int(hdr[4]), int(hdr[5]); h1 = mp.mpf(D)/M1
for _ in range(NG):
    w = f.readline().split(); b, nh = int(w[0]), int(w[1]); p0, dK, d1, d2, d3 = map(int, w[2:7])
    wt = mp.log(b)/(4*PI)
    ex0 = mp.frac(xs*wt)*two64; ph_err = max(ph_err, float(abs(p0 - ex0)))
    exK = mp.frac((xs + (Kan - 1)*mp.mpf(D))*wt)*two64              # exact phase at the last anchor, in units
    got = (p0 + (Kan - 1)*dK) % (1 << 64); dlt = float(min(abs(got - exK), two64 - abs(got - exK))); drift = max(drift, dlt)
    for _ in range(nh):
        k, lo, hi, ac = f.readline().split(); k = int(k); n = b**k if k > 1 or fd.isprime(b) else b
        if not fd.isprime(b): n = b
        lo, hi = float(lo), float(hi)
        worst_hi = min(worst_hi, float((mp.mpf(hi) - thi[n])/abs(thi[n])) if thi[n] != 0 else 0.0)
        worst_lo = min(worst_lo, float((tlo[n] - mp.mpf(lo))/abs(tlo[n])) if tlo[n] != 0 else float(tlo[n] - mp.mpf(lo)))
        nn += 1
print(json.dumps(dict(dat=dat, row=row, n_checked=nn, NG=NG, K_dat=Kan, K_own=int(mp.ceil((xe - xs)/mp.mpf(D))) + 1,
                      min_rel_hi_minus_true=worst_hi, min_rel_true_minus_lo=worst_lo,
                      D_up_logged=pj['D_up'], D_own=float(Dsum), tri_plus_logged=pj['tri_plus'], tri_plus_own=float(tri_p), tri_minus_logged=pj['tri_minus'], tri_minus_own=float(tri_m),
                      xs_logged=pj['xs'], xs_own=mp.nstr(xs, 30), p0_err_units=ph_err, drift_last_anchor_units=drift, drift_log2=math.log2(drift) if drift > 0 else None)))
