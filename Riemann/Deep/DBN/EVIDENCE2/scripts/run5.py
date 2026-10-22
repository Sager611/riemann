# run5.py -- production run: new paired-probe floors for every (cell, height, shift).
# usage: python run5.py SHARD NSHARDS OUT.jsonl [cells]
import sys, json, time, math, os, signal
from fractions import Fraction as Fr
import numpy as np
from flint import arb, fmpq  # pyright: ignore[reportMissingImports]
import lib5 as L, cert5 as C
SHARD, NSH = int(sys.argv[1]), int(sys.argv[2]); OUT = sys.argv[3]
SRC = os.path.join(os.environ.get('DBN4V_WORK', os.path.join(L._repo(), '.scratch', 'dbn4v')), 'arxiv/anc/certificate/inputs/refined') + '/'
NAMES = ['SOURCE_OLD_ENHANCED.json', 'SOURCE_RETUNED.json']
SH = ['9/20', '19/40', '49/100', '1/2', '51/100', '21/40', '107/200', '113/200', '23/40', '3/5']
SUBR = [(Fr(0), Fr(1, 250)), (Fr(1, 250), Fr(2, 125)), (Fr(2, 125), Fr(3, 50)), (Fr(3, 50), Fr(1, 5)),
        (Fr(1, 5), Fr(3, 5)), (Fr(3, 5), Fr(8, 5)), (Fr(8, 5), None)]
KSUB = [32, 32, 32, 16, 16, 16, 16]
plans = []; old = []
for nm in NAMES:
    try:
        with open(SRC + nm) as _f: d = json.load(_f)
    except (OSError, ValueError) as _e:
        raise SystemExit('cannot read source %s: %s' % (nm, _e))
    plans.append([Fr(r['height']) for r in d['plan']['rows']])
    old.append([{s: Fr(d['rows'][i]['r_' + s.replace('/', '_')]) for s in SH} for i in range(608)])
combos = []
for i in range(608):
    hs = sorted(set([plans[0][i], plans[1][i]]))
    for h in hs: combos.append((i, h))
if len(sys.argv) > 4: combos = [c for c in combos if c[0] in set(int(x) for x in sys.argv[4].split(','))]
mine = combos[SHARD::NSH]
done = set()
if os.path.exists(OUT):
    for line in open(OUT):
        r = json.loads(line); done.add((r['cell'], r['h']))
def dummyS(Gm): return {'Clo': np.ones(Gm.nb), 'Chi': np.ones(Gm.nb), 'Wc': Gm.W0, 'Ws': Gm.W0}
def rat_floor(x):
    lo = x.lower(); m, e = lo.man_exp(); f = Fr(int(m))*(Fr(2)**int(e))
    return Fr(math.floor(f*10**12), 10**12)
fo = open(OUT, 'a')
for (cell, h) in mine:
    if (cell, str(h)) in done: continue
    t0 = time.time()
    C.clear_cache()
    Psi = L.PsiTable(cell)
    geos = []
    for (s1, s2) in SUBR:
        R = L.Rect(h, cell, s1, s2); Gm = C.Geometry(R, Psi)
        Ct = C.one_point_fast(Gm, dummyS(Gm)); E1 = C.euler_one(R.rmax)
        f1 = C.floor_value(R, Gm, Ct, E1)
        geos.append((R, Gm, Ct, E1, f1))
    rec = {'cell': cell, 'h': str(h), 'shifts': {}}
    best_old = max(old[p][cell][s] for p in (0, 1) for s in SH if plans[p][cell] == h)
    hint = None
    for sh in SH:
        a = Fr(sh)*h
        R0 = geos[0][0]
        Ep = C.euler_pair(R0.rmin, R0.rmax, a)       # rho_max decreases with L1: valid for all subranges
        subs = []; target = None
        order = list(range(len(SUBR)))
        for si in order:
            R, Gm, Ct, E1, f1 = geos[si]
            f1p = C.floor_value(R, Gm, Ct, Ep)        # c = Ctop with the paired Euler term (valid)
            if target is not None and float(f1p.lower()) >= target:
                subs.append([str(SUBR[si][0]), str(SUBR[si][1]), Ct, Ct, 0, float(f1p.lower())]); continue
            S = Gm.shift_data(a)
            if target is None:
                res = C.paired_fast(Gm, S, Ct, K=KSUB[si], c_hint=hint)
                c = res[0] if res else Ct; K = KSUB[si] if res else 0
                hint = c if res else None
                fv = C.floor_value(R, Gm, c, Ep)
                target = float(fv.lower())
                subs.append([str(SUBR[si][0]), str(SUBR[si][1]), Ct, c, K, target])
                continue
            pe = C.phys_err(R, Gm.Aplus, Gm.Bplus)
            cn = 2*(R.L1/4 - arb(fmpq(1, 10**6)) - pe - Ep - arb(target))
            cneed = float(cn.lower())
            cneed = float(np.nextafter(cneed, -np.inf))
            ok = None
            if cneed > 0.3*Ct:
                ok = C.verify_cert(Gm, S, cneed, Ct, KSUB[si])
            if ok is not None:
                fv = C.floor_value(R, Gm, cneed, Ep)
                subs.append([str(SUBR[si][0]), str(SUBR[si][1]), Ct, cneed, KSUB[si], float(fv.lower())])
            else:
                res = C.paired_fast(Gm, S, Ct, K=32)
                c = res[0] if res else Ct; K = 32 if res else 0
                fv = C.floor_value(R, Gm, c, Ep)
                v = float(fv.lower()); target = min(target, v)
                subs.append([str(SUBR[si][0]), str(SUBR[si][1]), Ct, c, K, v])
        fl = min(x[5] for x in subs)
        flr = Fr(math.floor(Fr(fl)*10**12), 10**12)
        rec['shifts'][sh] = {'floor': str(flr), 'subs': subs, 'Ep': float(Ep.upper())}
    rec['seconds'] = round(time.time() - t0, 2)
    fo.write(json.dumps(rec) + '\n'); fo.flush()
    print(cell, str(h), 'min new', min(float(Fr(v['floor'])) for v in rec['shifts'].values()), 'old best', float(best_old), rec['seconds'], 's', flush=True)
fo.close()
