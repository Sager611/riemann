#!/usr/bin/env python3
"""mlimaudit: memory-lean, independently written recomputation of the mlimsup certificate.
Streams Weisz's zeros one row at a time (npzstream, no full np.load); per zero: gamma*y/(2 pi) at PH bits,
fractional part f, cos/sin via cos_pi/sin_pi(2f) at 320 bits (certify.py instead subtracts kk*2pi and calls cos).
Kernel autocorrelation R(d) via an fmpz_poly product (certify.py: direct double sum). Extra zeros at 128 bits
(certify.py: 96). Prints h_K(y+-) bounds, Weisz's JP values, radii, and checks the two claimed inequalities.
Usage: python3 lean_certify.py DATADIR [--nzeros n] [--noext] [--dump out.json]"""
import sys, os, time, json, hashlib
sys.set_int_max_str_digits(0)
from flint import arb, acb, acb_series, ctx, fmpz_poly
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from npzstream import NpzRows

t0 = time.time()
A = sys.argv[1:]
def opt(key, conv=str, default=None):
    if key in A:
        i = A.index(key); return conv(A[i + 1])
    return default
D = A[0]
NZ = opt('--nzeros', int, 14400)
NOEXT = '--noext' in A
DUMP = opt('--dump')
W = os.path.join(D, 'wbeni95-mertens-function-bounds-46f5314')
NPZ = os.path.join(D, 'zeta_zeros_14400.npz')
N = 14400
CLAIM = {'plus': '2.0739001417', 'minus': '-2.0764178466'}
WEISZ = {'plus': ('2.0617536514585282157830679158150887214080897256322036948605', '1.30e-59'),
         'minus': ('-2.0647155435291726514920038268907360718361372493174447898790', '5.59e-59')}

def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for ch in iter(lambda: f.read(1 << 20), b''):
            h.update(ch)
    return h.hexdigest()
sums = {}
for line in open(os.path.join(W, 'SHA256SUMS.txt')):
    if line.strip():
        hx, nm = line.split(); sums[nm.lstrip('*')] = hx
assert sha(NPZ) == sums['data/zeta_zeros_14400.npz'], 'npz sha256'
Ys = {}
for w in ('plus', 'minus'):
    p = os.path.join(W, 'points', 'y_%s.txt' % w)
    assert sha(p) == sums['points/y_%s.txt' % w], 'y sha256'
    Ys[w] = open(p).read().strip()
print('sha256 ok', flush=True)

def load_kernel(w):
    hdr, b = {}, []
    for line in open(os.path.join(HERE, '..', 'MLimsup', 'kernel_%s.txt' % w)):
        s = line.strip()
        if not s:
            continue
        if s.startswith('#'):
            if '=' in s:
                k_, v_ = s[1:].split('=', 1); hdr[k_.strip()] = v_.strip()
            continue
        b.append(int(s))
    M = int(hdr['M']); assert len(b) == M and min(b) > 0
    PQ = fmpz_poly(b) * fmpz_poly(b[::-1])          # coefficient M-1+d = sum_j b_j b_{j+d} = R(d)
    R = [int(PQ[M - 1 + d]) for d in range(M)] + [0, 0]
    for d in (0, 1, 2, 500, M - 2, M - 1):
        assert R[d] == sum(b[j] * b[j + d] for j in range(M - d))
    return {'M': M, 'hn': int(hdr['h_num']), 'hd': int(hdr['h_den']), 'R': R}
KER = {w: load_kernel(w) for w in ('plus', 'minus')}

def kval(K, g):
    """ball for k(g): linear interpolation of R(d)/R(0) at nodes d*h; 0 for g >= M h."""
    t = g * K['hd'] / K['hn']
    d = int(t.mid().floor().unique_fmpz())
    if d >= K['M']:
        assert t >= K['M']; return arb(0)
    assert t >= d and t <= d + 1, 'ball straddles a node'
    tau = t - d
    return (K['R'][d] * (1 - tau) + K['R'][d + 1] * tau) / K['R'][0]

def dball(mid, rad):
    """Weisz's ball rule: [mid - r - u, mid + r + u], u = one unit in the last printed digit."""
    s = mid.strip().lower(); m, _, e = s.partition('e')
    u = (int(e) if e else 0) - len(m.partition('.')[2])
    return arb(mid) + arb(0, (arb(rad) + arb('1e%d' % u)).upper())

# pre-pass: monotone ordinates, gamma_14400 (Weisz's JP support)
prev, last = 0.0, None
for i, gs in enumerate(NpzRows(NPZ, 'gamma')):
    v = float(gs[:25]); assert v > prev, ('not increasing', i); prev = v; last = gs
gr_last = list(NpzRows(NPZ, 'gamma_rad'))[-1]
print('pre-pass: 14400 ordinates strictly increasing; gamma_14400 = %s...' % last[:20], flush=True)

PH, PL = 33200, 320
ctx.prec = PH
inv2pi = 1 / (2 * arb.pi())
Y = {w: arb(Ys[w]) for w in Ys}
TJP = dball(last, gr_last)
del last
ctx.prec = PL
S = {k: arb(0) for k in ('k_plus', 'k_minus', 'jp_plus', 'jp_minus')}
maxrad = {'phase2f': 0.0, 'c': 0.0, 'term_k': 0.0}
picks = {0, 1, 2, 71, 999, 4999, 9999, 14398, 14399}
dump = {}
streams = [NpzRows(NPZ, nm) for nm in ('gamma', 'gamma_rad', 'c_re', 'c_re_rad', 'c_im', 'c_im_rad')]
for i, (gs, gr, crs, crr, cis, cir) in enumerate(zip(*streams)):
    if i >= NZ:
        break
    ctx.prec = PH
    g = dball(gs, gr)
    f2 = {}
    for w in Y:
        t = g * Y[w] * inv2pi                          # gamma*y/(2 pi)
        n = t.mid().floor()                           # exact integer (an arb)
        assert n.is_exact()
        f2[w] = 2 * (t - n)                           # theta/pi mod 2, rigorous
    u = g / TJP
    ctx.prec = PL
    g = +g; u = +u
    cre = dball(crs, crr); cim = dball(cis, cir)
    fjp = (1 - u) * u.cos_pi() + u.sin_pi() / arb.pi()
    for w in Y:
        a = +f2[w]
        re = cre * a.cos_pi() - cim * a.sin_pi()       # Re(c e^{i gamma y})
        kv = kval(KER[w], g)
        S['k_' + w] += 2 * kv * re
        S['jp_' + w] += 2 * fjp * re
        maxrad['phase2f'] = max(maxrad['phase2f'], float(a.rad()))
        maxrad['term_k'] = max(maxrad['term_k'], float((2 * kv * re).rad()))
        if i in picks:
            dump.setdefault(str(i + 1), {})[w] = {'theta_over_pi_mod2': a.mid().str(50), 'k': kv.str(30), 'term_k': (2 * kv * re).str(30), 'f_jp': fjp.str(30)}
    maxrad['c'] = max(maxrad['c'], float(cre.rad()), float(cim.rad()))
    if (i + 1) % 2000 == 0:
        print('  %d zeros, %.0f s' % (i + 1, time.time() - t0), flush=True)
del streams
print('known part done: %.0f s; max radii %s' % (time.time() - t0, maxrad), flush=True)

pen = {w: arb(0) for w in KER}
ext = {}
if not NOEXT and NZ == N:
    ctx.prec = 128
    T1 = {w: arb(KER[w]['hn'] * KER[w]['M']) / KER[w]['hd'] for w in KER}
    assert all(T1[w].is_exact() for w in KER) and T1['plus'] == T1['minus']
    T1v = T1['plus']
    NT = T1v.zeta_nzeros(); assert NT.is_exact(); NT = int(NT.unique_fmpz())
    g_last = +TJP
    rows, n, nused = [], N + 1, 0
    while n <= NT + 1:
        cnt = min(250, NT + 2 - n)
        for j, zz in enumerate(acb.zeta_zeros(n, cnt)):
            m = n + j; g = zz.imag
            assert zz.real == arb(0.5), 'zero not exactly on the line'
            if m == N + 1:
                assert g > g_last, 'zero 14401 not above Weisz gamma_14400'
            if m == NT + 1:
                assert g > T1v, 'zero N(T1)+1 below T1'; continue
            assert g < T1v, 'zero below N(T1) not below T1'
            zp = acb_series([zz, 1], prec=2).zeta().coeffs()[1]
            ca = 1 / (abs(zz) * abs(zp))
            for w in KER:
                pen[w] += 2 * abs(kval(KER[w], g)) * ca
            nused += 1
            rows.append('%d %s %s %s' % (m, g.mid().str(25), arb(ca.upper()).str(15), float(g.rad())))
        n += cnt
        print('  ext zeros up to %d, %.0f s' % (n - 1, time.time() - t0), flush=True)
    assert nused == NT - N
    with open(os.path.join(HERE, 'ext_zeros_128.txt'), 'w') as fo:
        fo.write('# n gamma_mid |c|_upper gamma_rad  (Arb 128 bits, acb.zeta_zeros; N(T1)=%d)\n' % NT)
        fo.write('\n'.join(rows) + '\n')
    ext = {'T1': T1v.str(10), 'N(T1)': NT, 'used': nused, 'pen_plus': pen['plus'].str(20), 'pen_minus': pen['minus'].str(20)}
ctx.prec = PL
out = {k: v.str(45) for k, v in S.items()}
out.update({'rad_' + k: float(v.rad()) for k, v in S.items()})
out['ext'] = ext
lo = S['k_plus'] - pen['plus']; hi = S['k_minus'] + pen['minus']
out['h_plus_lower'] = lo.str(30); out['h_minus_upper'] = hi.str(30)
out['h_plus_lower_endpoint'] = arb(lo.lower()).str(25); out['h_minus_upper_endpoint'] = arb(hi.upper()).str(25)
out['rad_h_plus'] = float(lo.rad()); out['rad_h_minus'] = float(hi.rad())
print(json.dumps(out, indent=1))
for w in ('plus', 'minus'):
    if NZ == N:
        wb = arb(WEISZ[w][0]) + arb(0, arb(WEISZ[w][1]).upper())
        print('JP %s vs Weisz certificate: overlap=%s, |diff mid|=%s' % (w, S['jp_' + w].overlaps(wb), abs(S['jp_' + w].mid() - wb.mid()).str(5)))
if NZ == N and not NOEXT:
    c1 = bool(lo > arb(CLAIM['plus'])); c2 = bool(hi < arb(CLAIM['minus']))
    print('h_K+(y+) > %s: %s ; h_K-(y-) < %s: %s' % (CLAIM['plus'], c1, CLAIM['minus'], c2))
    print('LEAN RECOMPUTATION: ' + ('BOTH CLAIMED INEQUALITIES HOLD' if c1 and c2 else 'CLAIM NOT REPRODUCED'))
if DUMP:
    json.dump(dump, open(DUMP, 'w'), indent=1)
print('time %.0f s' % (time.time() - t0))
