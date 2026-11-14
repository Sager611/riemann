# K2(v,T0,k) at Bi-Hou's (v, T0, k) from Gallego's certified zero data (Zenodo 23700212, .crho files).
# usage: python3 k2_gallego.py DATA_DIR OUT_JSON ; needs crho_format.py (Gallego's reader) on sys.path.
# Phases gamma*v mod 2pi in exact integer arithmetic; terms in float64, summed with math.fsum.
import sys, math, json, time, random, os
from decimal import Decimal, getcontext
from crho_format import CrhoFile
getcontext().prec = 120
DATA, OUT = sys.argv[1], sys.argv[2]
FILES = ['cc_00000001.crho', 'cc_00240000.crho', 'cc_00305000.crho', 'cc_00371494.crho', 'cc_00742987.crho',
         'cc_01114480.crho', 'cc_01485973.crho', 'cc_01857466.crho', 'cc_02228959.crho']
V_STR = '3484971332733132396.5383341489440461'          # Bi-Hou evaluation point (exact decimal)
Y_GAL = '3484971332733132396.538348536063730810969705971170'  # Gallego's point
V_INT, V_SC = 34849713327331323965383341489440461, 16     # v = V_INT / 10^16
K = 1.5e-10; T0 = 1300000
P = 30
PI_DEC = Decimal('3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955058223172535940812848111745028410270193852110555964462294895493038196')
_twopi = {}
def twopi_int(a):
    if a not in _twopi: _twopi[a] = int((2 * PI_DEC * (Decimal(10) ** (a + V_SC + P))).to_integral_value(rounding='ROUND_FLOOR'))
    return _twopi[a]
def dec_val(flags, lint, lfrac, e, mant):
    x = mant / 10.0 ** lfrac if lfrac < 300 else 0.0
    if flags & 4: x *= 10.0 ** e
    return -x if flags & 1 else x
class Acc:
    """Exact-ish running sum: math.fsum over chunks, then over chunk sums (bounded memory)."""
    def __init__(self): self.buf = []; self.parts = []
    def append(self, x):
        self.buf.append(x)
        if len(self.buf) >= 20000: self.parts.append(math.fsum(self.buf)); self.buf = []
    def total(self): return math.fsum(self.parts + [math.fsum(self.buf)])
dvs = {'m1e-6': -1e-6, 'p1e-6': 1e-6, 'gallego_y': float(Decimal(Y_GAL) - Decimal(V_STR)), 'm1e-5': -1e-5, 'p1e-5': 1e-5}
terms = Acc(); terms_d = {k: Acc() for k in dvs}; dterms = Acc(); absw = Acc(); errc = Acc()
parts = {}; t_start = time.time(); last_n = 0; last_g = 0.0; n_below = 0; g_last_below = None; g_first_above = None
top = []; sample_rows = []; rng = random.Random(20261114)
crit = {1, 2, 98625, 98626, 600000, 600001, 2325571, 2325572, 2325573}
for fn in FILES:
    cf = CrhoFile(os.path.join(DATA, fn)); c = cf.col; rad_g = float(cf.radius_gamma_txt); t1 = time.time()
    for i in range(cf.count):
        n = cf.n0 + i
        if n <= last_n: continue
        assert n == last_n + 1, (fn, n, last_n)
        gf, gl, ga, ge = int(c['gamma_flags'][i]), int(c['gamma_lint'][i]), int(c['gamma_lfrac'][i]), int(c['gamma_exp'][i])
        assert gf == 2, (n, gf)
        G = int.from_bytes(c['gamma_mant'][i].tobytes(), 'little')
        g = G / 10.0 ** ga
        assert g > last_g, (n, g, last_g)
        if g >= T0:
            g_first_above = (n, cf.txt('gamma', i)); break
        n_below = n; g_last_below = (n, cf.txt('gamma', i))
        re = dec_val(int(c['re_flags'][i]), int(c['re_lint'][i]), int(c['re_lfrac'][i]), int(c['re_exp'][i]), int.from_bytes(c['re_mant'][i].tobytes(), 'little'))
        im = dec_val(int(c['im_flags'][i]), int(c['im_lint'][i]), int(c['im_lfrac'][i]), int(c['im_exp'][i]), int.from_bytes(c['im_mant'][i].tobytes(), 'little'))
        rc = float(c['radio_c'][i])
        tp = twopi_int(ga)
        r = (G * V_INT * 10 ** P) % tp
        th = r / 10.0 ** (ga + V_SC + P)
        w2 = 2.0 * math.exp(-K * g * g)
        cs, sn = math.cos(th), math.sin(th)
        t = w2 * (re * cs - im * sn); terms.append(t)
        dterms.append(-w2 * g * (re * sn + im * cs))
        for kk, dv in dvs.items():
            th2 = th + g * dv; terms_d[kk].append(w2 * (re * math.cos(th2) - im * math.sin(th2)))
        a = w2 * math.hypot(re, im); absw.append(a); errc.append(w2 * math.sqrt(2.0) * rc)
        if len(top) < 40 or a > top[-1][0]:
            top.append((a, n, g, t)); top.sort(reverse=True); top = top[:40]
        if n in crit or rng.random() < 120.0 / 2325572:
            sample_rows.append({'n': n, 'gamma': cf.txt('gamma', i), 're': cf.txt('re', i), 'im': cf.txt('im', i), 'radio_c': rc, 'radio_gamma': cf.radius_gamma_txt})
        for lab, cut in (('n<=98625', 98625), ('n<=600000', 600000)):
            if n == cut: parts[lab] = (terms.total(), absw.total())
        last_n, last_g = n, g
    else:
        print(fn, 'rows', cf.count, 'done to n', last_n, '%.1fs' % (time.time() - t1), flush=True); del cf; continue
    print(fn, 'stopped at T0, n', last_n, flush=True); break
# rows straddling T0
sample_rows.append({'n': g_first_above[0], 'gamma': g_first_above[1], 'note': 'first above T0'})
K2 = terms.total(); S_abs = absw.total(); vrad = float(Decimal(V_STR)) * rad_g
res = {'k': K, 'T0': T0, 'v': V_STR, 'N_below_T0': n_below, 'gamma_last_below': g_last_below, 'gamma_first_above': g_first_above,
       'K2_float': repr(K2), 'K2_main_n<=600000': repr(parts['n<=600000'][0]), 'K2_tail_n>600000': repr(K2 - parts['n<=600000'][0]),
       'abs_tail_n>600000': repr(S_abs - parts['n<=600000'][1]), 'K2_n<=98625': repr(parts['n<=98625'][0]),
       'S_abs=sum 2w|c|': repr(S_abs), 'rigorous_err_phase=S_abs*v*rad_gamma': repr(S_abs * vrad), 'rigorous_err_c': repr(errc.total()),
       'dK2/dv': repr(dterms.total()), 'K2_at_offsets': {kk: repr(v.total()) for kk, v in terms_d.items()},
       'top_weights': [(n, g, a, t) for a, n, g, t in top[:15]], 'seconds': time.time() - t_start}
json.dump(res, open(OUT, 'w'), indent=1)
json.dump(sample_rows, open(OUT.replace('.json', '_sample.json'), 'w'), indent=0)
print(json.dumps({k: v for k, v in res.items() if k != 'top_weights'}, indent=1))
