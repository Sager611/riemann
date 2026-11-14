# Rescore Bi-Hou's 25 sub-record candidates (search/results/subrecord_peaks_summary.csv) at k=1.5e-10 with Gallego's zeros
# n <= 742986 (gamma < 4.6e5; omitted tail < 1e-12). One pass: K2, K2', K2'' at y_base; quadratic polish K2max = K2 - K2'^2/(2K2'').
import sys, math, csv, os, json
from decimal import Decimal, getcontext
from crho_format import CrhoFile
getcontext().prec = 120
DATA, CSV, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
FILES = ['cc_00000001.crho', 'cc_00240000.crho', 'cc_00305000.crho', 'cc_00371494.crho']
K = 1.5e-10; P = 30
PI = Decimal('3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955058223172535940812848111745028410270193852110555964462294895493038196')
cands = [r for r in csv.DictReader(open(CSV))]
ys = [r['y_exact_base'] for r in cands] + ['3484971332733132396.5383341489440461']   # + Bi-Hou's certified v
Ys = []
for y in ys:
    ip, fp = y.split('.'); Ys.append((int(ip + fp), len(fp)))
tw = {}
def twopi(sc):
    if sc not in tw: tw[sc] = int((2 * PI * Decimal(10) ** sc).to_integral_value(rounding='ROUND_FLOOR'))
    return tw[sc]
def dec_val(flags, lfrac, e, mant):
    x = mant / 10.0 ** lfrac
    if flags & 4: x *= 10.0 ** e
    return -x if flags & 1 else x
class Acc:
    def __init__(self): self.buf = []; self.parts = []
    def add(self, x):
        self.buf.append(x)
        if len(self.buf) >= 20000: self.parts.append(math.fsum(self.buf)); self.buf = []
    def tot(self): return math.fsum(self.parts + [math.fsum(self.buf)])
S0 = [Acc() for _ in ys]; S1 = [Acc() for _ in ys]; S2 = [Acc() for _ in ys]
last = 0
for fn in FILES:
    cf = CrhoFile(os.path.join(DATA, fn)); c = cf.col
    for i in range(cf.count):
        n = cf.n0 + i
        if n <= last: continue
        last = n
        ga = int(c['gamma_lfrac'][i]); G = int.from_bytes(c['gamma_mant'][i].tobytes(), 'little'); g = G / 10.0 ** ga
        re = dec_val(int(c['re_flags'][i]), int(c['re_lfrac'][i]), int(c['re_exp'][i]), int.from_bytes(c['re_mant'][i].tobytes(), 'little'))
        im = dec_val(int(c['im_flags'][i]), int(c['im_lfrac'][i]), int(c['im_exp'][i]), int.from_bytes(c['im_mant'][i].tobytes(), 'little'))
        w2 = 2.0 * math.exp(-K * g * g)
        for j, (Y, f) in enumerate(Ys):
            sc = ga + f + P; th = ((G * Y * 10 ** P) % twopi(sc)) / 10.0 ** sc
            cs, sn = math.cos(th), math.sin(th)
            S0[j].add(w2 * (re * cs - im * sn)); S1[j].add(-w2 * g * (re * sn + im * cs)); S2[j].add(-w2 * g * g * (re * cs - im * sn))
    del cf
res = []
for j, y in enumerate(ys):
    k0, k1, k2 = S0[j].tot(), S1[j].tot(), S2[j].tot()
    dy = -k1 / k2 if k2 != 0 else 0.0
    kmax = k0 - k1 * k1 / (2 * k2) if k2 != 0 else k0
    bh = cands[j]['h_gauss_k1.5e-10_600k'] if j < len(cands) else 'certified v'
    res.append({'y': y, 'K2_at_y': k0, 'dy_polish': dy, 'K2_polished': kmax, 'BiHou_csv_k1.5e-10_600k': bh})
    print('%-36s K2(y)=%+.9f  dy=%+.2e  K2max=%+.9f  BH=%s' % (y, k0, dy, kmax, bh))
json.dump(res, open(OUT, 'w'), indent=1)
