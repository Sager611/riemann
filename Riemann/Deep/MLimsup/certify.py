#!/usr/bin/env python3
"""Rigorous certificate (FLINT/Arb ball arithmetic) for the mlimsup bounds
    lim sup M(x)/sqrt(x) > CLAIM_PLUS   and   lim inf M(x)/sqrt(x) < CLAIM_MINUS.

Theorem used (Ingham 1942; Odlyzko-te Riele 1985, Sec. 2): if K >= 0 is even, K in C^2, K(y) = O((1+y^2)^-1),
k(t) = int K(y) e^{ity} dy has k(0) = 1 and k(t) = 0 for |t| >= T1, and the zeros with 0 < gamma < T1 are on the
critical line and simple, then for every real y
    liminf M(x)/sqrt(x) <= h_k(y) <= limsup M(x)/sqrt(x),   h_k(y) = 2 sum_{0<gamma<T1} k(gamma) Re(c_rho e^{i gamma y}),
c_rho = 1/(rho zeta'(rho)).  (If RH fails or a zero is multiple, both limits are infinite: the bounds are unconditional.)

Kernel (kernel_<point>.txt): g = step function on [-T1/2, T1/2], M cells of width h = h_num/h_den, integer heights b_j;
k = (g*g~)/||g||^2, i.e. k = linear interpolation of R(d)/R(0) at t = d h, R(d) = sum_j b_j b_{j+d}; k = 0 for |t| >= T1 = M h.
Then K = |g^|^2/(2 pi ||g||^2) >= 0, even, entire, O(y^-2).

Point y and the first 14400 zeros (gamma to ~9860 digits, c to 60 digits, as balls): B. Weisz, Zenodo 23651459 (v2.0.0).
Zeros 14401..N(T1): computed here with Arb (acb.zeta_zeros: rigorous isolation on the critical line, count by Turing's
method; N(T1) by acb zeta_nzeros); their phases at y are unknown, so each such term is bounded by 2|k(gamma)||c|:
    h_k(y) >= S_known - P   (plus point),   h_k(y) <= S_known + P   (minus point).
Cross-check: the script re-evaluates Weisz's Jurkat-Peyerimhoff values at the same points.

Usage: python3 certify.py [datadir] [--nzeros n] [--phases out.npz]     (datadir default .scratch/mlimsup/data)
"""
import sys, os, io, time, json, hashlib, zipfile, urllib.request
import numpy as np
from flint import arb, acb, acb_series, ctx

CLAIMS = {'plus': '2.0739001417', 'minus': '-2.0764178466'}
sys.set_int_max_str_digits(0)
t00 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
args = sys.argv[1:]
opt = {}
for key in ('--nzeros', '--phases'):
    if key in args:
        i = args.index(key); opt[key] = args[i + 1]; args = args[:i] + args[i + 2:]
D = args[0] if args else '.scratch/mlimsup/data'
URL = 'https://zenodo.org/api/records/23651459/files/wbeni95/mertens-function-bounds-v2.0.0.zip/content'
os.makedirs(D, exist_ok=True)
zp = os.path.join(D, 'weisz.zip')
if not os.path.exists(zp):
    with urllib.request.urlopen(urllib.request.Request(URL, headers={'User-Agent': 'mlimsup-certify'}), timeout=900) as r, open(zp, 'wb') as out:
        out.write(r.read())
zf = zipfile.ZipFile(zp)
root = zf.namelist()[0].split('/')[0]
def zread(name):
    return zf.read(root + '/' + name)
sums = {l.split()[1].lstrip('*'): l.split()[0] for l in zread('SHA256SUMS.txt').decode().splitlines() if l.strip()}
npz_bytes = zread('data/zeta_zeros_14400.npz')
assert hashlib.sha256(npz_bytes).hexdigest() == sums['data/zeta_zeros_14400.npz'], 'data checksum mismatch'
Y = {}
for w in ('plus', 'minus'):
    b = zread('points/y_%s.txt' % w)
    assert hashlib.sha256(b).hexdigest() == sums['points/y_%s.txt' % w]
    Y[w] = b.decode().strip()
z = np.load(io.BytesIO(npz_bytes), allow_pickle=False)
gam_s, gam_r = z['gamma'], z['gamma_rad']
cre_s, cre_r, cim_s, cim_r = z['c_re'], z['c_re_rad'], z['c_im'], z['c_im_rad']     # load each array once
del npz_bytes
N = 14400
NZ = int(opt.get('--nzeros', N))
print('Weisz data sha256 ok (zeros npz, y+, y-)', flush=True)

def load_kernel(path):
    hdr, b = {}, []
    for line in open(path):
        line = line.strip()
        if not line: continue
        if line.startswith('#'):
            if '=' in line and line[1:].split('=')[0].strip() in ('M', 'h_num', 'h_den'):
                k_, v_ = line[1:].split('=', 1); hdr[k_.strip()] = v_.strip()
            continue
        b.append(int(line))
    M = int(hdr['M']); assert len(b) == M
    R = [sum(b[j] * b[j + d] for j in range(M - d)) for d in range(M)] + [0, 0]
    return {'M': M, 'h_num': int(hdr['h_num']), 'h_den': int(hdr['h_den']), 'R': R}
KER = {w: load_kernel(os.path.join(HERE, 'kernel_%s.txt' % w)) for w in ('plus', 'minus')}

def dec_ulp(s):
    s = s.strip().lower(); mant, _, ex = s.partition('e'); frac = mant.partition('.')[2]
    return '1e%d' % ((int(ex) if ex else 0) - len(frac))
def ball(mid, rad):
    mid, rad = str(mid), str(rad)
    return arb(mid) + arb(0, arb(rad).upper()) + arb(0, arb(dec_ulp(mid)).upper())

P_LO = 320
ctx.prec = P_LO
Rb = {w: [arb(r) for r in KER[w]['R']] for w in KER}
hb = {w: arb(KER[w]['h_num']) / KER[w]['h_den'] for w in KER}
def kernel_at(w, g):
    """ball for k_w(g), g a ball > 0"""
    K = KER[w]; t = g / hb[w]
    d = int(t.mid().floor().unique_fmpz())
    if d >= K['M']:
        assert t.lower() >= K['M']
        return arb(0)
    assert t.lower() >= d and t.upper() <= d + 1, 'node ambiguity'
    tau = t - d
    return (Rb[w][d] * (1 - tau) + Rb[w][d + 1] * tau) / Rb[w][0]

idig = max(len(Y[w].lstrip('-').split('.')[0]) for w in Y)
p_hi = max(int(idig * 3.3220) + 400, int(max(len(str(s)) for s in gam_s[:50]) * 3.3220) + 64)
ctx.prec = p_hi
yb = {w: arb(Y[w]) for w in Y}
twopi = 2 * arb.pi()
T_JP = ball(gam_s[N - 1], gam_r[N - 1])
ctx.prec = P_LO
pi = arb.pi()
tot = {key: arb(0) for key in ('jp_plus', 'jp_minus', 'k_plus', 'k_minus')}
ph_float = {w: np.zeros(NZ) for w in Y}
for i in range(NZ):
    ctx.prec = p_hi
    g = ball(gam_s[i], gam_r[i])
    red = {}
    for w in Y:
        th = g * yb[w]
        qm, qe = (th.mid() / twopi.mid()).mid().man_exp()
        qm, qe = int(qm), int(qe)
        kk = arb(qm >> (-qe) if qe < 0 else qm << qe)
        assert kk.is_exact()
        red[w] = th - kk * twopi                      # rigorous ball for gamma*y - 2 pi kk
    u = g / T_JP
    ctx.prec = P_LO
    g_lo = +g; u = +u
    f = (1 - u) * (pi * u).cos() + (pi * u).sin() / pi
    cre = ball(cre_s[i], cre_r[i]); cim = ball(cim_s[i], cim_r[i])
    for w in Y:
        r_ = +red[w]
        re_part = cre * r_.cos() - cim * r_.sin()      # Re(c e^{i gamma y})
        ph_float[w][i] = float(r_.mid())
        tot['jp_' + w] += 2 * f * re_part
        tot['k_' + w] += 2 * kernel_at(w, g_lo) * re_part
    if (i + 1) % 2000 == 0:
        print('  %d zeros, %.0f s' % (i + 1, time.time() - t00), flush=True)

# zeros above gamma_14400 inside the kernel support: worst case
pen = {w: arb(0) for w in KER}
ext = {}
if NZ == N:
    ctx.prec = 96
    T1 = {w: arb(KER[w]['h_num'] * KER[w]['M']) / KER[w]['h_den'] for w in KER}
    T1max = max(T1.values(), key=lambda v: float(v))
    NT = T1max.zeta_nzeros()
    assert NT.is_exact(), 'zero count at T1 not exact'
    NT = int(NT.unique_fmpz())
    g14400 = ball(gam_s[N - 1], gam_r[N - 1])
    n = N + 1
    nmax = {w: N for w in KER}
    while n <= NT:
        zs = acb.zeta_zeros(n, min(500, NT + 1 - n))
        for j, zz in enumerate(zs):
            g = zz.imag
            if n + j == N + 1:
                assert g > g14400 and abs(zz.real - arb('0.5')) == 0
            zpr = acb_series([zz, 1], prec=2).zeta().coeffs()[1]
            cab = 1 / (abs(zz) * abs(zpr))
            for w in KER:
                if g.lower() < T1[w]:
                    pen[w] += 2 * abs(kernel_at(w, g)) * cab
                    nmax[w] = n + j
        n += len(zs)
        print('  ext zeros up to %d, %.0f s' % (n - 1, time.time() - t00), flush=True)
    for w in KER:
        ext[w] = {'T1': T1[w].str(12), 'N(T1)': int(T1[w].zeta_nzeros().unique_fmpz()), 'last_zero_used': nmax[w], 'penalty': pen[w].str(20, radius=True)}
        assert ext[w]['N(T1)'] == nmax[w], 'zero list below T1 incomplete'
ctx.prec = P_LO
hk = {'plus': tot['k_plus'] - pen['plus'], 'minus': tot['k_minus'] + pen['minus']}
res = {key: v.str(40, radius=True) for key, v in tot.items()}
res.update({'ext_' + w: ext.get(w) for w in KER})
res.update({'bound_plus (h_k(y+) >=)': hk['plus'].str(30, radius=True), 'bound_minus (h_k(y-) <=)': hk['minus'].str(30, radius=True)})
print(json.dumps(res, indent=1))
if '--phases' in opt:
    gf = np.array([float(str(s)[:30]) for s in gam_s[:NZ]])
    np.savez(opt['--phases'], theta_plus=ph_float['plus'], theta_minus=ph_float['minus'], gam=gf,
             c_re=np.array([float(s) for s in cre_s[:NZ]]), c_im=np.array([float(s) for s in cim_s[:NZ]]))
if NZ == N:
    ok_jp = bool(abs(tot['jp_plus'] - arb('2.0617536514585282157830679158150887214080897256322036948605')) < arb('1e-40')) and \
            bool(abs(tot['jp_minus'] - arb('-2.0647155435291726514920038268907360718361372493174447898790')) < arb('1e-40'))
    print('cross-check: Weisz JP values reproduced to 1e-40:', ok_jp)
    c1 = bool(hk['plus'] > arb(CLAIMS['plus'])); c2 = bool(hk['minus'] < arb(CLAIMS['minus']))
    print('claim: limsup M(x)/sqrt(x) >= h_k+(y+) > %s: %s' % (CLAIMS['plus'], c1))
    print('claim: liminf M(x)/sqrt(x) <= h_k-(y-) < %s: %s' % (CLAIMS['minus'], c2))
    print('ALL CLAIMS VERIFIED' if (c1 and c2 and ok_jp) else 'CLAIM NOT VERIFIED')
print('time %.0f s' % (time.time() - t00))
