#!/usr/bin/env python3
# dbnrec prep_end.py: Arb (python-flint) coefficient enclosures of row i's end-of-cell head Q^+ on the shell x in [X e^u0, X e^u1],
# exact uint64 phase data for bb_rec.c (format of dbnhead prep.py), and the chord difference D (LEMMA.md part B).
# c_n(t,x) = (1/2) d_t(n) n^{-1/2-tL/4} F_i(n), L = log(x/4pi), F_i(n) = sum_j w_j n^{-h_j/2} cos(a_j log n/2)  (U.tex l.1044-1047, 1595-1604)
# mags (all n <= 4096, Lemma src:coeff: d_t >= 0 nondecreasing in t; primes: d_t(p) = log p b_t(p) exactly, l.1272):
#   hi+ = upper bound of m_n on {t = t+, L >= Llo}; hi- = same at t = t-; lo = lower bound of m_n on the whole cell and L <= Lhi.
#   primes: m_p = (1/2) log p p^{-1/2} exp(t log p (log p - L)/4) |F|: hi = value at (t, Llo), lo = value at (t+, Lhi).
#   others: hi = (1/2) d_{t+}(n) n^{-1/2 - t Llo/4} |F|_hi (chord in t is above; t = t- or t+), lo = (1/2) d_{t-}(n) n^{-1/2 - t+ Lhi/4} |F|_lo.
# Every bound is an Arb ball endpoint, rounded outward to double by nextafter.
# usage: prep_end.py ARC AFFINE row u0 u1 T out.dat   (u0, u1, T decimal strings; thr = T + s u0 written for bb)
import sys, json, math
import mpmath as mp
from flint import arb, fmpq, ctx
ctx.prec = 160; mp.mp.dps = 50
arcf, aff, ri, u0s, u1s, Ts, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7]
K = 4096
row = json.load(open(arcf))[ri]
assert row['joint_jet_overlay'] is None, 'jet rows not handled'
ar = [r for r in json.load(open(aff))['rows'] if int(r['row']) == ri][0]
M = fmpq(*map(int, ar['mass'].split('/'))); beta = fmpq(*map(int, ar['beta'].split('/')))
s = (1 - beta)*M/4
assert s == fmpq(*map(int, ar['retained_main_mass'].split('/')))/4
def Q(v): return fmpq(*map(int, v.split('/'))) if '/' in v else fmpq(int(v))
nodes = [(Q(n['height']), Q(n['shift']), Q(n['weight'])) for n in row['nodes']]
assert sum(w for _, _, w in nodes) == M, 'mass = signed node mass'
tm, tp = Q(row['time'][0]), Q(row['time'][1])
def dec(v): return fmpq(int(round(float(v)*10**9)), 10**9) if v != '0' else fmpq(0)
def fq(v): return int(v.p)/int(v.q)
def mq(v): return mp.mpf(int(v.p))/int(v.q)
u0, u1, T = dec(u0s), dec(u1s), dec(Ts)
assert u0 >= 0 and u1 > u0
X = arb(6000000185827); L0 = (X/(4*arb.pi())).log()
Llo, Lhi = L0 + arb(u0), L0 + arb(u1)
ln = [arb(0), arb(0)] + [arb(n).log() for n in range(2, K + 1)]
def bt(t, n): return (arb(t)*ln[n]*ln[n]/4).exp()
def dcoef(t):
    d = [arb(0), arb(0)] + [ln[n]*bt(t, n) for n in range(2, K + 1)]
    b = [arb(0), arb(0)] + [bt(t, n) for n in range(2, K + 1)]
    for r in range(2, K//2 + 1):
        for m in range(2, K//r + 1): d[r*m] -= d[r]*b[m]
    return d
dm, dp = dcoef(tm), dcoef(tp)
isprime = [False, False] + [all(n % q for q in range(2, int(n**0.5) + 1)) for n in range(2, K + 1)]
def up(a): return math.nextafter(float(a.upper()), math.inf)
def dn(a): return math.nextafter(float(a.lower()), -math.inf)
lo, hi, act, hm_minus, hm_plus = {}, {}, {}, arb(0), arb(0)
nsign0 = 0
for n in range(2, K + 1):
    F = sum((arb(w)*(-arb(h)*ln[n]/2).exp()*(arb(a)*ln[n]/2).cos() for h, a, w in nodes), arb(0))
    Fhi = abs(F).upper()
    if F > 0: sg, Flo = 1, F.lower()
    elif F < 0: sg, Flo = -1, (-F).lower()
    else: sg, Flo = 0, arb(0); nsign0 += 1
    if isprime[n]:
        base = ln[n]*(-ln[n]/2).exp()/2
        hp = base*(arb(tp)*ln[n]*(ln[n] - Llo)/4).exp()*Fhi
        hmn = base*(arb(tm)*ln[n]*(ln[n] - Llo)/4).exp()*Fhi
        lw = base*(arb(tp)*ln[n]*(ln[n] - Lhi)/4).exp()*Flo
    else:
        hp = dp[n]*(-(arb(1)/2 + arb(tp)*Llo/4)*ln[n]).exp()*Fhi/2
        hmn = dp[n]*(-(arb(1)/2 + arb(tm)*Llo/4)*ln[n]).exp()*Fhi/2
        lw = dm[n]*(-(arb(1)/2 + arb(tp)*Lhi/4)*ln[n]).exp()*Flo/2
    hm_plus += hp; hm_minus += hmn
    H, Lw = up(hp), max(0.0, dn(lw))
    if sg > 0: lo[n], hi[n] = Lw, H
    elif sg < 0: lo[n], hi[n] = -H, -Lw
    else: lo[n], hi[n] = -H, H
    act[n] = float(dp[n]*(-(arb(1)/2 + arb(tp)*Llo/4)*ln[n]).exp()*F/2)   # witness only (t+, x = X e^u0)
Dch = hm_minus - hm_plus                     # chord difference D = sum(hi- - hi+) >= 0
tri_p, tri_m = up(hm_plus), up(hm_minus)
thr = fq(T + s*u0)
def factor(m):
    f = {}; p = 2
    while p*p <= m:
        while m % p == 0: f[p] = f.get(p, 0) + 1; m //= p
        p += 1
    if m > 1: f[m] = f.get(m, 0) + 1
    return f
groups = {}
for n in range(2, K + 1):
    f = factor(n)
    if len(f) == 1: (p, e), = f.items(); groups.setdefault(p, []).append((e, lo[n], hi[n], act[n]))
    else: groups[n] = [(1, lo[n], hi[n], act[n])]
D = 4*math.pi/math.log(2)
M1, M2, M3 = 32, 8, 8
Xm = mp.mpf(6000000185827)
xs = Xm*mp.exp(mq(u0)) + mp.pi*mq(tm)/4; xe = Xm*mp.exp(mq(u1)) + mp.pi*mq(tp)/4
Kan = int(mp.ceil((xe - xs)/D)) + 1
Dm = mp.mpf(D); h1 = Dm/M1; h2 = h1/M2; h3 = h2/M3
two64 = mp.mpf(2)**64
def u64(v): return int(mp.floor(mp.frac(v)*two64)) % (1 << 64)
def u64r(v): return int(mp.nint(mp.frac(v)*two64)) % (1 << 64)
order = sorted(groups, key=lambda b: -sum(max(abs(h[1]), abs(h[2])) for h in groups[b]))
order.remove(2); order = [2] + order
with open(out, 'w') as fo:
    fo.write('%d %d %.17g %d %d %d %.17g %.17g %.17g %.17g\n' % (len(order), Kan, D, M1, M2, M3, thr, 0.0, tri_p, up(Dch)))
    for b in order:
        w = mp.log(b)/(4*mp.pi)
        fo.write('%d %d %d %d %d %d %d %.17g\n' % (b, len(groups[b]), u64(xs*w), u64r(Dm*w), u64r(h1*w), u64r(h2*w), u64r(h3*w), math.log(b)/2))
        for h in groups[b]: fo.write('%d %.17g %.17g %.17g\n' % h)
print(json.dumps(dict(row=ri, t=[str(tm), str(tp)], u0=str(u0), u1=str(u1), T=str(T), s=str(s), s_float=fq(s), thr=thr, K=Kan, groups=len(order),
                      xs=mp.nstr(xs, 30), tri_plus=tri_p, tri_minus=tri_m, D_up=up(Dch), D_lo=dn(Dch), A=fq(Q(row['first_hi'])), B=fq(Q(row['last_hi'])),
                      nsign0=nsign0, prec=ctx.prec)))
