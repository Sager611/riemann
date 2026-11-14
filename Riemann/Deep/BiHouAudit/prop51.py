# Recheck Bi-Hou Prop 5.1: min over cells I_nu covering [1.3e6, 4.3e6] of N_nu = sum_{mu != nu, dist(gamma_mu, I_nu) < 1} log dist,
# using Gallego's certified zeros (gamma radius 5e-26). Reads only the gamma columns of the .crho files.
import sys, os, json, lzma, struct, math, time
import numpy as np
from decimal import Decimal, getcontext
getcontext().prec = 60
DATA = sys.argv[1]
FILES = ['cc_02228959.crho', 'cc_02600452.crho', 'cc_02971945.crho', 'cc_03343438.crho', 'cc_03714927.crho',
         'cc_04800350.crho', 'cc_05885773.crho', 'cc_06971196.crho', 'cc_08056619.crho']
NLO, NHI = 2325500, 8600000          # index window (N(1.3e6) = 2325572; N(4.3e6) ~ 8.51e6)
t0 = time.time()
gam = np.zeros(NHI - NLO + 1); have = np.zeros(NHI - NLO + 1, dtype=bool); exact = {}
for fn in FILES:
    with open(os.path.join(DATA, fn), 'rb') as f:
        assert f.read(6) == b'CRHO1\n'
        L = struct.unpack('<Q', f.read(8))[0]; h = json.loads(f.read(L).decode())
        cols = {}
        for c in h['columnas']:
            if c['nombre'].startswith('gamma'):
                a = np.frombuffer(lzma.decompress(f.read(c['bytes'])), dtype=np.dtype(c['dtype']))
                cols[c['nombre']] = a.reshape(-1, c['ancho']) if c.get('ancho') else a
            else:
                f.seek(c['bytes'], 1)
    n0 = h['n0']; cnt = h['filas']; assert np.all(cols['gamma_flags'] == 2)
    lf = cols['gamma_lfrac']; mant = cols['gamma_mant']
    for i in range(cnt):
        n = n0 + i
        if n < NLO or n > NHI or have[n - NLO]: continue
        G = int.from_bytes(mant[i].tobytes(), 'little'); gam[n - NLO] = G / 10 ** int(lf[i]); have[n - NLO] = True
    print(fn, n0, cnt, '%.0fs' % (time.time() - t0), flush=True)
    del cols, mant, lf
assert have.all(), 'missing rows'
assert np.all(np.diff(gam) > 0), 'not increasing'
theta = (gam[:-1] + gam[1:]) / 2            # theta[j] between gam[j], gam[j+1]
# cells nu: I = [theta[nu-1], theta[nu]] (array index j for gamma), covering [1.3e6, 4.3e6]
j_lo = int(np.searchsorted(theta, 1.3e6)); j_hi = int(np.searchsorted(theta, 4.3e6)) + 1
assert theta[j_lo - 1] <= 1.3e6 and theta[j_hi - 1] >= 4.3e6
J = 10; cand = []; ncell = j_hi - j_lo; CH = 400000
for st in range(j_lo, j_hi, CH):
    jj = np.arange(st, min(st + CH, j_hi)); lo = theta[jj - 1]; hi = theta[jj]; Nn = np.zeros(len(jj))
    for d in range(1, J + 1):
        dl = lo - gam[jj - d]; dr = gam[jj + d] - hi
        Nn += np.where(dl < 1, np.log(np.maximum(dl, 1e-300)), 0.0) + np.where(dr < 1, np.log(np.maximum(dr, 1e-300)), 0.0)
        if d == J: assert dl.min() >= 1 and dr.min() >= 1, 'J too small'
    o = np.argsort(Nn)[:5]
    cand += [(float(Nn[q]), int(NLO + jj[q]), float(gam[jj[q]]), float(lo[q]), float(hi[q])) for q in o]
    del jj, lo, hi, Nn, dl, dr
cand.sort()
Nmin, nu, gnu, lo0, hi0 = cand[0]
print('cells', ncell, ' nu range', NLO + j_lo, NLO + j_hi - 1, ' min N_nu = %.6f at nu = %d (gamma_nu = %.5f, cell [%.5f, %.5f])' % (Nmin, nu, gnu, lo0, hi0))
for c in cand[:5]: print('  N = %.6f  nu = %d  gamma = %.5f' % (c[0], c[1], c[2]))
print('M = 14 - min N_nu = %.6f   (Bi-Hou: N_nu >= -8.76707, M = 22.7671; 6 185 380 cells; nu = 5 042 995)  %.0fs' % (14 - Nmin, time.time() - t0))
