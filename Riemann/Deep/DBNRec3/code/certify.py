#!/usr/bin/env python3
# dbnrec certify.py: assemble the certified head fields (A', B') of one row from its bb_rec shells (LEMMA.md Part B), in Arb.
# Per shell j (tags r<row>_<jj>): every unresolved final cell of bb_rec is re-evaluated in interval arithmetic (Arb): UB_C = sum_g g_g(theta_g(x'_c))
#   (upper endpoints of the balls lo*cos, hi*cos) + LIP (r + eps_loc); pruned/resolved cells are <= thr_j (bb_rec float margin 1e-9).
#   S_j = max(thr_j - s u0_j, max_C [UB_C - s log((x'_c - r - eps_loc - pi t+/4)/X)]); a_j = S_j + D_j, b_j = S_j.
# Far range x >= X e^{uF}: a_F = tri-(last shell) - s uF, b_F = tri+(last shell) - s uF (chord mags decrease in L_lo).
# A' = max a, B' = max b, rounded up to the 2^-96 grid.  usage: certify.py ROW tag1 tag2 ...   (prints JSON)
import sys, json, math
from flint import arb, fmpq, fmpz, ctx
ctx.prec = 128
L = 'Riemann/Deep/DBNRec3/logs/'; S_ = '.scratch/dbnrec3/'
row = int(sys.argv[1]); tags = sys.argv[2:]
src = json.load(open(S_ + 'ws/lambda_round15/flow_main_v2/MAIN_CLOCK_SOURCE.json'))['records'][row]
assert int(src['row']) == row and src['mode'] == 'retained'
def Q(v): return fmpq(*map(int, v.split('/'))) if '/' in v else fmpq(int(v))
def exact(x):   # a double as an exact rational
    n, d = float(x).as_integer_ratio(); return fmpq(n, d)
def upmax(a, b):
    ua, ub = a.upper(), b.upper(); return ua if ua >= ub else ub
A0, B0 = Q(src['first_hi']), Q(src['last_hi'])
X = arb(6000000185827); PI = arb.pi()
def load_dat(path):
    f = open(path); hdr = f.readline().split(); NG = int(hdr[0]); D = float(hdr[2]); gs = []
    for _ in range(NG):
        b, nh = map(int, f.readline().split()[:2]); hs = []
        for _ in range(nh):
            k, lo, hi, ac = f.readline().split(); hs.append((int(k), exact(lo), exact(hi)))
        gs.append((b, hs))
    return D, gs
res = dict(row=row, A_raymond=str(A0), B_raymond=str(B0), shells=[])
prev = None; avals = []; bvals = []
for tag in tags:
    pj = json.loads(open(L + tag + '.prep.json').readline())
    bbo = open(L + tag + '.bb.out').read()
    assert 'CERT sup' in bbo and ('DONE K %d ' % pj['K']) in bbo, tag + ': bb run incomplete or K mismatch'
    s = Q(pj['s']); u0 = Q(pj['u0']); u1 = Q(pj['u1']); tp = Q(pj['t'][1])
    assert (prev is None and u0 == 0) or prev == u0, tag + ': shells not contiguous'
    prev = u1
    thr = exact(pj['thr']); D, gs = load_dat(S_ + tag + '.dat'); Dx = exact(D)
    om = [arb(b).log()/2 for b, _ in gs]
    LIP = arb(0)
    for i, (b, hs) in enumerate(gs):
        LIP += arb(sum((k*max(abs(lo), abs(hi)) for k, lo, hi in hs), fmpq(0)))*om[i]
    xs = arb(pj['xs'])
    Sj = arb(thr - s*u0).upper(); ncell = 0; worst = None; dmax = 0.0
    for line in open(L + tag + '.unres.txt'):
        w = line.split()
        if len(w) == 2:   # old dump: centre offset from xs (double; error < 1e-4), length h3/16
            xc = xs + arb(w[0]); r = arb(Dx)/2048/32; eps = arb(fmpq(1, 10**4)); ubf = float(w[1])
        else:             # new dump: anchor k, in-anchor centre offset, length
            xc = xs + arb(Dx)*int(w[0]) + arb(w[1]); r = arb(w[2])/2 + arb(fmpq(1, 10**9)); eps = arb(fmpq(1, 10**12)); ubf = float(w[3])
        tot = arb(0)
        for i, (b, hs) in enumerate(gs):
            th = xc*om[i]
            for k, lo, hi in hs:
                c = (th*k).cos(); tot += upmax(arb(lo)*c, arb(hi)*c)
        ub = tot + LIP*(r + eps)
        lg = ((xc - r - eps - PI*arb(tp)/4)/X).log()
        v = (ub - arb(s)*lg).upper()
        dmax = max(dmax, float(ub.upper()) - ubf)
        if v > Sj: worst = dict(cell=line.strip(), u=float(lg.mid()), ub_arb=float(ub.upper()), ub_float=ubf)
        Sj = upmax(Sj, v)
        ncell += 1
    Dj = exact(pj['D_up'])
    a = (Sj + arb(Dj)).upper(); b = Sj.upper()
    avals.append(a); bvals.append(b)
    res['shells'].append(dict(tag=tag, u0=str(u0), u1=str(u1), thr=pj['thr'], cells=ncell, S=float(b), D=float(Dj), a=float(a), b=float(b),
                              max_ub_arb_minus_float=dmax, worst=worst))
    last = pj
uF = Q(last['u1']); s = Q(last['s'])
aF = (arb(exact(last['tri_minus'])) - arb(s*uF)).upper(); bF = (arb(exact(last['tri_plus'])) - arb(s*uF)).upper()
avals.append(aF); bvals.append(bF)
res['far'] = dict(uF=str(uF), a=float(aF), b=float(bF))
def ceil96(x):   # strictly above x, on the 2^-96 grid
    n = (x*arb(2)**96).upper().floor().unique_fmpz()
    assert n is not None
    return fmpq(int(n) + 1, 2**96)
A1, B1 = avals[0], bvals[0]   # exact points (Arb upper endpoints): exact maxima
for z in avals[1:]: A1 = z if z > A1 else A1
for z in bvals[1:]: B1 = z if z > B1 else B1
for z in avals: assert not (z > A1)
for z in bvals: assert not (z > B1)
Ap, Bp = ceil96(A1), ceil96(B1)
assert Ap >= Bp > 0
res.update(A_new=str(Ap), B_new=str(Bp), A_new_float=int(Ap.p)/int(Ap.q), B_new_float=int(Bp.p)/int(Bp.q),
           dA=int((A0 - Ap).p)/int((A0 - Ap).q), dB=int((B0 - Bp).p)/int((B0 - Bp).q), prec=ctx.prec)
print(json.dumps(res))
