#!/usr/bin/env python3
"""Independent spot check of Weisz's zero data (Zenodo 23651459, data/zeta_zeros_14400.npz) with FLINT/Arb.

(a) Indexing and low digits, for every checked index n: acb.zeta_zero(n) at 256 bits (the n-th zero, isolated by Arb with
    Turing's method) and c = 1/(rho zeta'(rho)); the stored gamma and c balls must agree to 1e-60 / 1e-55.
(b) Full precision, for the indices in HIGH (stored ball rule: midpoint +- (stored radius + one unit in the last digit)):
    B = box of radius r ~ 2 x (stored gamma radius) around 1/2 + i*gamma_mid; one interval Newton step
    N(B) = s0 - zeta(s0)/zeta'(B) at ~33,000 bits.
    - N(B) inside B: B contains exactly one zero; B is symmetric about Re s = 1/2, so that zero is on the line;
      with (a) it is the n-th zero.
    - Im N(B) inside the stored gamma ball: the stored gamma ball is correct.
    - 1/(B zeta'(B)) inside the stored Re c and Im c balls: the stored c balls are correct.
Selection: --high i,j,... (explicit), --heavy K (union over the two kernels kernel_plus/minus.txt of the K largest
2 k(gamma)|c| among zeros 1..14400), --random R (R further random indices, seed 20261106). --dry prints the list only.
Usage: python3 spotcheck.py [datadir] [--high list] [--heavy K] [--random R] [--dry]"""
import sys, os, io, time, zipfile, random
import numpy as np
from flint import arb, acb, acb_series, ctx
sys.set_int_max_str_digits(0)
HERE = os.path.dirname(os.path.abspath(__file__))
args = sys.argv[1:]
opt = {}
for key in ('--high', '--heavy', '--random'):
    if key in args:
        i = args.index(key); opt[key] = args[i + 1]; args = args[:i] + args[i + 2:]
dry = '--dry' in args
args = [a for a in args if a != '--dry']
D = args[0] if args else '.scratch/mlimsup/data'
t0 = time.time()
zf = zipfile.ZipFile(os.path.join(D, 'weisz.zip')); root = zf.namelist()[0].split('/')[0]
z = np.load(io.BytesIO(zf.read(root + '/data/zeta_zeros_14400.npz')), allow_pickle=False)
gam_s, gam_r = z['gamma'], z['gamma_rad']
cre_s, cre_r, cim_s, cim_r = z['c_re'], z['c_re_rad'], z['c_im'], z['c_im_rad']
N = 14400
def dec_ulp(s):
    s = s.strip().lower(); mant, _, ex = s.partition('e'); frac = mant.partition('.')[2]
    return '1e%d' % ((int(ex) if ex else 0) - len(frac))
def ball(mid, rad):
    mid, rad = str(mid), str(rad)
    return arb(mid) + arb(0, arb(rad).upper()) + arb(0, arb(dec_ulp(mid)).upper())
def inside(x, b):            # ball x certainly inside ball b
    return bool(x.lower() >= b.lower()) and bool(x.upper() <= b.upper())

high = set(int(v) for v in opt['--high'].split(',')) if '--high' in opt else set()
heavy_info = ''
if '--heavy' in opt:
    K = int(opt['--heavy'])
    gam = np.array([float(str(s)[:30]) for s in gam_s])
    absc = np.abs(np.array([float(s) for s in cre_s]) + 1j * np.array([float(s) for s in cim_s]))
    sets = []
    for w in ('plus', 'minus'):
        hdr, b = {}, []
        for line in open(os.path.join(HERE, 'kernel_%s.txt' % w)):
            line = line.strip()
            if line.startswith('#'):
                if '=' in line and line[1:].split('=')[0].strip() in ('M', 'h_num', 'h_den'):
                    k_, v_ = line[1:].split('=', 1); hdr[k_.strip()] = int(v_)
            elif line:
                b.append(float(line))
        M = hdr['M']; h = hdr['h_num'] / hdr['h_den']; b = np.array(b)
        r = np.correlate(b, b, 'full')[M - 1:]; r = np.append(r / r[0], [0.0, 0.0])
        x = gam / h; i0 = np.floor(x).astype(int); fr = x - i0
        k = np.where(i0 < M, r[np.minimum(i0, M)] * (1 - fr) + r[np.minimum(i0 + 1, M + 1)] * fr, 0.0)
        wgt = 2 * k * absc
        top = np.argsort(-wgt)[:K] + 1
        sets.append(set(int(v) for v in top))
        heavy_info += ' %s: top-%d weight %.4f of %.4f;' % (w, K, np.sort(wgt)[::-1][:K].sum(), wgt.sum())
    heavy = sets[0] | sets[1]
    high |= heavy
    heavy_info = 'heavy union %d indices (max %d, max gamma %.1f);' % (len(heavy), max(heavy), max(gam[n - 1] for n in heavy)) + heavy_info
rand = []
if '--random' in opt:
    random.seed(20261106)
    pool = [n for n in range(1, N + 1) if n not in high]
    rand = sorted(random.sample(pool, int(opt['--random'])))
    high |= set(rand)
high = sorted(high)
print('selection: %s random %s' % (heavy_info, rand), flush=True)
print('HIGH (%d): %s' % (len(high), high), flush=True)
if dry:
    sys.exit(0)

random.seed(20261106)
LOW = sorted(set([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 14399, 14400] + random.sample(range(11, 14399), 28)) | set(high))
ctx.prec = 256
bad = []
for n in LOW:
    rho = acb.zeta_zero(n)
    zp = acb_series([rho, 1], prec=2).zeta().coeffs()[1]
    c = 1 / (rho * zp)
    g = ball(str(gam_s[n - 1])[:90], '1e-80')
    cr = ball(cre_s[n - 1], cre_r[n - 1]); ci = ball(cim_s[n - 1], cim_r[n - 1])
    ok = abs(rho.imag - g) < arb('1e-60') and abs(c.real - cr) < arb('1e-55') and abs(c.imag - ci) < arb('1e-55')
    if not ok: bad.append(n)
print('(a) %d zeros at 256 bits (1-10, 14399, 14400, 28 random, all of HIGH): %s  (%.0f s)'
      % (len(LOW), 'all agree' if not bad else 'MISMATCH at %s' % bad, time.time() - t0), flush=True)
res = {'newton': 0, 'gamma_ball': 0, 'c_ball': 0}
for n in high:
    t1 = time.time()
    gs = str(gam_s[n - 1]); gr = str(gam_r[n - 1])
    prec = int(len(gs) * 3.3220) + 256
    ctx.prec = prec
    g = ball(gs, gr)
    r = g.rad() * 2 + arb('1e-%d' % (len(gs) - 8))
    s0 = acb(arb('0.5'), g.mid())
    B = acb(arb('0.5', r), arb(g.mid(), r))
    zs = acb_series([s0], prec=1).zeta().coeffs()[0]
    zpB = acb_series([B, 1], prec=2).zeta().coeffs()[1]
    Nb = s0 - zs / zpB
    in_B = bool(abs(Nb.real - arb('0.5')) < r) and bool(abs(Nb.imag - g.mid()) < r)
    in_g = inside(Nb.imag, g)
    cenc = 1 / (B * zpB)
    cr = ball(cre_s[n - 1], cre_r[n - 1]); ci = ball(cim_s[n - 1], cim_r[n - 1])
    in_c = inside(cenc.real, cr) and inside(cenc.imag, ci)
    res['newton'] += in_B; res['gamma_ball'] += in_B and in_g; res['c_ball'] += in_B and in_c
    print('(b) zero %5d (gamma %.2f): N(B) in B: %s; Im N(B) in stored gamma ball: %s; c in stored c balls: %s (%d bits, %.0f s)'
          % (n, float(g.mid()), in_B, in_g, in_c, prec, time.time() - t1), flush=True)
print('SUMMARY: %d zeros at full precision: N(B) in B %d/%d, stored gamma ball certified %d/%d, stored c balls certified %d/%d; (a) %s; total %.0f s'
      % (len(high), res['newton'], len(high), res['gamma_ball'], len(high), res['c_ball'], len(high),
         'ok' if not bad else 'MISMATCH', time.time() - t0), flush=True)
