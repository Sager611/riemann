#!/usr/bin/env python3
# dbnrecaudit scan_wrap.py: independent global search of shell 0 (u in [0, 1/100]) of one row for phase alignments (falsification of bb_rec's
# sup, nothing certified). Own mpmath amplitudes (falsify_dense.py formulas) and own uint64 phase data; C helper 'scan' ranks all K anchors
# by the groups b in {3,5,7,11} (table lookups) plus the 2-group's value (constant over anchors) at NOFF in-anchor offsets; the top NTOP anchors are then evaluated in full
# (all 4095 terms of Q+, 4000 points over the whole anchor, C helper 'dense'). Compares with bb's CERT sup and the cert's binding cell.
# usage: scan_wrap.py ROW NTOP [NOFF=8]   (cwd = repo root; writes Riemann/Deep/DBNRecAudit/logs/scan_r<ROW>_<NOFF>off.json)
import sys, json, subprocess, os, math
sys.argv, args = sys.argv[:1], sys.argv[1:]
sys.path.insert(0, 'Riemann/Deep/DBNRecAudit/code')
import falsify_dense as fd          # no rows given: only its definitions are loaded
mp = fd.mp
from fractions import Fraction as Fr
row, ntop = int(args[0]), int(args[1]); noff = int(args[2]) if len(args) > 2 else 8
p = fd.plan[row]; K_ = fd.K; ln = fd.ln; X = fd.X; PI = fd.PI; Dd = fd.Dd
tm, tp = (mp.mpf(Fr(v).numerator)/Fr(v).denominator for v in p['time'])
nodes = [(Fr(n['height']), Fr(n['shift']), Fr(n['weight'])) for n in p['nodes']]
q = lambda r: mp.mpf(r.numerator)/r.denominator
F = [mp.mpf(0)]*(K_ + 1)
for n in range(2, K_ + 1): F[n] = sum(q(w)*mp.exp(-q(h)*ln[n]/2)*mp.cos(q(a)*ln[n]/2) for h, a, w in nodes)
dm, dp = fd.dcoef(tm), fd.dcoef(tp)
L0 = mp.log(X/(4*PI)); Llo, Lhi = L0, L0 + mp.mpf(1)/100
lo, hi = {}, {}
for n in range(2, K_ + 1):
    aF = abs(F[n])
    if fd.isprime(n):
        base = ln[n]*n**mp.mpf(-0.5)/2; H = base*mp.exp(tp*ln[n]*(ln[n] - Llo)/4)*aF; l = base*mp.exp(tp*ln[n]*(ln[n] - Lhi)/4)*aF
    else:
        H = dp[n]*mp.exp(-(mp.mpf(1)/2 + tp*Llo/4)*ln[n])*aF/2; l = dm[n]*mp.exp(-(mp.mpf(1)/2 + tp*Lhi/4)*ln[n])*aF/2
    lo[n], hi[n] = (float(l), float(H)) if F[n] > 0 else (float(-H), float(-l))
xs = X + PI*tm/4; xe = X*mp.exp(mp.mpf(1)/100) + PI*tp/4
K = int(mp.ceil((xe - xs)/mp.mpf(Dd))) + 1
two64 = mp.mpf(2)**64
def u64(v): return int(mp.floor(mp.frac(v)*two64)) % (1 << 64)
groups = {}
for b in (2, 3, 5, 7, 11):
    e, m = 1, b
    while m <= K_: groups.setdefault(b, []).append((e, m)); e += 1; m *= b
def gval(b, th):   # th in turns
    s = 0.0
    for e, m in groups[b]:
        c = math.cos(2*math.pi*e*th); s += hi[m]*c if c >= 0 else lo[m]*c
    return s
w2 = mp.log(2)/(4*PI); th20 = float(mp.frac(xs*w2))
grid = [Dd*i/20000 for i in range(20000)]
g2 = [gval(2, th20 + d*float(w2)) for d in grid]
g2max = max(g2); high = [d for d, v in zip(grid, g2) if v >= g2max - 0.08]
d0, d1 = min(high), max(high)
offs = [d0 + (d1 - d0)*(i + 0.5)/noff for i in range(noff)]
fn = '/tmp/dbnra_keep/scan_in.txt'
with open(fn, 'w') as fo:
    fo.write('%d %d %d\n' % (len(groups) - 1, K, len(offs)))
    fo.write(' '.join('%.12g' % gval(2, th20 + o*float(w2)) for o in offs) + '\n')   # 2-group: constant over anchors (drift < 10^-6 turn)
    for b in groups:
        if b == 2: continue
        w = mp.log(b)/(4*PI)
        fo.write('%d %d\n' % (u64(xs*w), int(mp.nint(mp.frac(mp.mpf(Dd)*w)*two64)) % (1 << 64)))
        fo.write(' '.join(str(int(mp.nint(mp.frac(mp.mpf(o)*w)*two64)) % (1 << 64)) for o in offs) + '\n')
        fo.write(' '.join('%.12g' % gval(b, (i + 0.5)/4096) for i in range(4096)) + '\n')
out = subprocess.run(['/tmp/dbnra_keep/scan', fn, str(ntop)], capture_output=True, text=True, check=True).stdout.split('\n')
top = sorted(([int(a), int(o), float(v)] for a, o, v in (l.split() for l in out if l.strip())), key=lambda z: -z[2])
os.remove(fn)
res = dict(row=row, K=K, g2_high_region=[d0, d1], offsets=offs, ntop=len(top), top10=top[:10], full=[])
binding = 3244584625
res['binding_rank'] = next((i for i, z in enumerate(top) if z[0] == binding), None)
fb = '/tmp/dbnra_keep/scan_b.txt'
for k, o, v in top[:min(len(top), ntop)]:
    xk = xs + k*mp.mpf(Dd)
    with open(fb, 'w') as fo:
        fo.write('%d\n' % (K_ - 1))
        for n in range(2, K_ + 1): fo.write('%.17g %.17g %.17g %.17g\n' % (float(ln[n]/2), lo[n], hi[n], float(mp.fmod(xk*ln[n]/2, 2*PI))))
    r = subprocess.run(['/tmp/dbnra_keep/dense', fb, 'b', '0', repr(Dd), '4000'], capture_output=True, text=True, check=True).stdout.split()
    res['full'].append(dict(k=k, P=v, Qplus_max=float(r[0]), at=float(r[1])))
os.remove(fb)
res['full'].sort(key=lambda z: -z['Qplus_max'])
json.dump(res, open('Riemann/Deep/DBNRecAudit/logs/scan_r%d_%doff.json' % (row, noff), 'w'), indent=0)
print(json.dumps(dict(row=row, K=K, binding_rank=res['binding_rank'], top3=top[:3], best_full=res['full'][:5])))
