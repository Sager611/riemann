#!/usr/bin/env python3
"""Checks behind the NO RECORD verdict of Riemann/Deep/MertensX/REPORT.md (mertensx, 2026-11-05).

Rerun (about 2 minutes plus a 91 MB download; single core, < 400 MB RAM):
  OMP_NUM_THREADS=1 uv run --no-project --with python-flint==0.9.0 --with numpy \
      python3 Riemann/Deep/MertensX/check.py [datadir]
datadir defaults to .scratch/mertensx/data (delete it afterwards).

1. N(74000) = 98625 (Arb), so h_StR of Kim-Nguyen sums the first 98,625 zeros.
2. Spot check of the Kim-Nguyen zero data (Zenodo 10775723) against Arb at 40 indices:
   gamma_n (acb.zeta_zero), alpha_n = 1/|rho zeta'(rho)|, psi_n = arg(rho zeta'(rho)).
3. Ball-arithmetic value of h_StR(y) = 2 sum_{gamma<74000} alpha e^{-3e-9 gamma^2} cos(gamma y - psi)
   at the Kim-Nguyen point y = 19571878850562201959.215107 (from their data).
4. Every published candidate below the current record exp(3.485e18): Kim-Nguyen's enumeration
   output (140-119.OUT) and Bi-Hou's sub-record list; h_StR at the point and after local
   maximisation in y (float64 phases from exact integer reduction mod 2 pi).
"""
import sys, os, time, csv, math, urllib.request
from decimal import Decimal, getcontext, ROUND_FLOOR
import numpy as np
from flint import arb, acb, acb_series, ctx

t0 = time.time()
D = sys.argv[1] if len(sys.argv) > 1 else '.scratch/mertensx/data'
KN = 'https://zenodo.org/api/records/10775723/files/%s/content'
BH = 'https://raw.githubusercontent.com/bzyzz/mertens_bound/HEAD/search/results/subrecord_peaks_summary.csv'
SRC = {'zeroes.txt': KN % 'zeroes.txt', 'alpha.txt': KN % 'alpha.txt', 'beta.txt': KN % 'beta.txt',
       '140-119.OUT': KN % '140-119.OUT', 'bh_peaks.csv': BH}
os.makedirs(D, exist_ok=True)
for f, u in SRC.items():
    p = os.path.join(D, f)
    if not os.path.exists(p):
        req = urllib.request.Request(u, headers={'User-Agent': 'mertensx-check'})
        with urllib.request.urlopen(req, timeout=600) as r, open(p, 'wb') as out:
            out.write(r.read())

NZ = 98625
K = 3e-9
RECORD = 3.484971332733178124e18          # Bi-Hou / Gallego: x < exp(3 484 971 332 733 178 124)
KN_Y = '19571878850562201959.215107'


def head(f, n):
    out = []
    with open(os.path.join(D, f)) as fh:
        for i, line in enumerate(fh):
            if i >= n:
                break
            out.append(line.strip())
    return out


G = head('zeroes.txt', NZ + 1)
A = head('alpha.txt', NZ)
P = head('beta.txt', NZ)

# 1. zero count
ctx.prec = 200
n74 = arb(74000).zeta_nzeros()
print('1. Arb N(74000) =', n74, '| data gamma_98625 = %.6f, gamma_98626 = %.6f' % (float(G[NZ - 1]), float(G[NZ])))
assert n74 == NZ and float(G[NZ - 1]) < 74000 < float(G[NZ])
G = G[:NZ]

# 2. spot check against Arb
idx = list(range(1, 21)) + [int(round(x)) for x in np.linspace(21, NZ, 20)]
dg = da = dp = 0.0
for n in idx:
    z = acb.zeta_zero(n)
    w = z * acb_series([z, 1], prec=2).zeta().coeffs()[1]
    dg = max(dg, abs(float((z.imag - arb(G[n - 1])).mid())))
    da = max(da, abs(float((1 / abs(w) - arb(A[n - 1])).mid())) / float(A[n - 1]))
    dp = max(dp, abs(float((w.arg() - arb(P[n - 1])).mid())))
print('2. spot check of %d zeros vs Arb: max |d gamma| = %.1e, max rel |d alpha| = %.1e, max |d psi| = %.1e'
      % (len(idx), dg, da, dp))

# 3. ball-arithmetic h_StR at the Kim-Nguyen point
ka = arb('3e-9')
y = arb(KN_Y)
s = arb(0)
for g, a, p in zip(G, A, P):
    g = arb(g)
    s += arb(a) * (-(ka * g * g)).exp() * (g * y - arb(p)).cos()
print('3. h_StR(%s) = %s  (Arb, 200 bits, Kim-Nguyen data)' % (KN_Y, (2 * s).str(12)))

# 4. candidates below the current record, local maximisation in y
getcontext().prec = 400
SCALE = 45


def fix(sv, e):
    return int((Decimal(sv).scaleb(e)).to_integral_value(rounding=ROUND_FLOOR))


Gi = [fix(g, SCALE) for g in G]
Pi = [fix(p, SCALE) for p in P]
TP = fix((2 * arb.pi()).str(120, radius=False), SCALE + 20)       # 2 pi * 10^65
gam = np.array([float(g) for g in G])
wgt = np.array([float(a) for a in A]) * np.exp(-K * gam * gam)


def phases(ys):
    Y = fix(ys, 20)
    s20 = 10 ** 20
    return np.array([float((gi * Y - pi * s20) % TP) for gi, pi in zip(Gi, Pi)]) / float(TP) * 2 * math.pi


def hval(ph, d):
    return 2.0 * np.dot(wgt, np.cos(ph + gam * d))


def polish(ys, W=1e-3, step=1e-6):
    ph = phases(ys)
    h0 = hval(ph, 0.0)
    sg = 1.0 if h0 > 0 else -1.0
    ds = np.arange(-W, W + step / 2, step)
    hs = np.empty(len(ds))
    for i in range(0, len(ds), 64):
        blk = ds[i:i + 64]
        hs[i:i + 64] = 2.0 * (np.cos(ph[None, :] + np.outer(blk, gam)) @ wgt)
    jb = int(np.argmax(sg * hs))
    best = (ds[jb], hs[jb])
    for j in np.argsort(-sg * hs)[:8]:                     # Newton from the best grid points
        d = ds[j]
        for _ in range(30):
            arg = ph + gam * d
            h1 = -2.0 * np.dot(wgt * gam, np.sin(arg))
            h2 = -2.0 * np.dot(wgt * gam * gam, np.cos(arg))
            if h2 * sg >= 0:
                break
            dn = d - h1 / h2
            if abs(dn - d) < 1e-15:
                d = dn
                break
            d = dn
        hv = hval(ph, d)
        if abs(d) <= 2 * W and sg * hv > sg * best[1]:
            best = (d, hv)
    return h0, best[1], best[0]


cands = []
clusters = {}
for line in open(os.path.join(D, '140-119.OUT')):
    p = line.split()
    if len(p) >= 10:
        yv = float(p[9])
        key = round(yv)
        if yv < RECORD and (key not in clusters or abs(float(p[8])) > abs(clusters[key][0])):
            clusters[key] = (float(p[8]), p[9])
for key, (h, ys) in sorted(clusters.items()):
    cands.append(('KN 140-119.OUT', ys, None))
for r in csv.DictReader(open(os.path.join(D, 'bh_peaks.csv'))):
    if float(r['y_exact_base']) < RECORD:
        cands.append(('Bi-Hou list', r['y_exact_base'], float(r['abs_h_gauss_k1.5e-10_600k'])))

print('4. published candidates y < %.6e (record), k = 3e-9, 98,625 zeros:' % RECORD)
print('   %-15s %-32s %10s %10s %10s %10s' % ('source', 'y', 'h(y)', 'polished', 'dy', '|K2| 1.5e-10'))
for src, ys, k15 in [('KN record', KN_Y, None)] + cands:
    h0, hp, d = polish(ys)
    print('   %-15s %-32s %+10.6f %+10.6f %+10.2e %10s' % (src, ys[:32], h0, hp, d,
          '%.6f' % k15 if k15 else '-'))
print('done in %.0f s' % (time.time() - t0))
