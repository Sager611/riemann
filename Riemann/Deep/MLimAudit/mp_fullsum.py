#!/usr/bin/env python3
"""mlimaudit: library-independent full recomputation of the known-part sums S+- (zeros 1..14400) and Weisz's JP values in
mpmath (floating point, not interval): theta/(2 pi) mod 1 at 33,300 bits from Weisz's gamma midpoints, then cos/sin, kernel and
sum at 80 digits. Kernel k from R(d) = sum_j b_j b_{j+d} computed here. Compares with the Arb values of lean_certify.log.
Usage: python3 mp_fullsum.py DATADIR [N]"""
import sys, os, time, json
sys.set_int_max_str_digits(0)
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from npzstream import NpzRows
t0 = time.time()
D = sys.argv[1]; N = int(sys.argv[2]) if len(sys.argv) > 2 else 14400
W = os.path.join(D, 'wbeni95-mertens-function-bounds-46f5314'); NPZ = os.path.join(D, 'zeta_zeros_14400.npz')
Ys = {w: open(os.path.join(W, 'points', 'y_%s.txt' % w)).read().strip() for w in ('plus', 'minus')}
def load_R(w):
    b = [int(l) for l in open(os.path.join(HERE, '..', 'MLimsup', 'kernel_%s.txt' % w)) if l.strip() and not l.startswith('#')]
    M = len(b)
    return [sum(b[j] * b[j + d] for j in range(M - d)) for d in range(M)] + [0, 0]
R = {w: load_R(w) for w in Ys}
HP, LP = 33300, 80
mp.mp.prec = HP
inv2pi = 1 / (2 * mp.pi)
yv = {w: mp.mpf(Ys[w]) for w in Ys}
gN = None
for i, gs in enumerate(NpzRows(NPZ, 'gamma')):
    gN = gs
mp.mp.dps = LP
TJP = mp.mpf(gN[:100]); h = mp.mpf(1355) / 100
S = {k: mp.mpf(0) for k in ('k_plus', 'k_minus', 'jp_plus', 'jp_minus')}
streams = [NpzRows(NPZ, nm) for nm in ('gamma', 'c_re', 'c_im')]
for i, (gs, cr, ci) in enumerate(zip(*streams)):
    if i >= N: break
    mp.mp.prec = HP
    g = mp.mpf(gs)
    fr = {}
    for w in yv:
        t = g * yv[w] * inv2pi
        fr[w] = t - mp.floor(t)
    mp.mp.dps = LP
    fr = {w: +fr[w] for w in fr}
    gl = mp.mpf(gs[:100]); cre = mp.mpf(cr); cim = mp.mpf(ci)
    u = gl / TJP
    fjp = (1 - u) * mp.cospi(u) + mp.sinpi(u) / mp.pi
    tt = gl / h; d = int(mp.floor(tt)); tau = tt - d
    for w in yv:
        re = cre * mp.cospi(2 * fr[w]) - cim * mp.sinpi(2 * fr[w])
        kv = (R[w][d] * (1 - tau) + R[w][d + 1] * tau) / R[w][0] if d < 1350 else mp.mpf(0)
        S['k_' + w] += 2 * kv * re
        S['jp_' + w] += 2 * fjp * re
    if (i + 1) % 2000 == 0:
        print('  %d zeros, %.0f s' % (i + 1, time.time() - t0), flush=True)
out = {k: mp.nstr(v, 50) for k, v in S.items()}
print(json.dumps(out, indent=1))
json.dump(out, open(os.path.join(HERE, 'mp_fullsum.json'), 'w'), indent=1)
print('time %.0f s' % (time.time() - t0))
