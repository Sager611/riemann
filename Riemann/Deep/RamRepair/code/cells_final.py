# cells_final.py -- certified log-density bound per cell (Padhi Thm 6.1 / Lemma 6.2 with our corrections; NOTE.md sections 3.3-3.4):
#   log N(1 - k/1e4, e^l) <= log(nt d4 d5 kap2 C/(w (c1^2 - d4 c5r))) + 2 x delta_k l - log(delta_k l)
# with the cell parts of cells_mine.cell_parts (out/cell_parts.json) and d4 from d4_mine.d4_bound (own Barban-Vehov bound,
# own constants out/d4_consts.json, k1/k0 from out/kconst_cert.json).  Cells with c1 <= 0 or c1^2 - d4 c5r < 0.2 are dropped.
# Every parameter is used as the same exact double everywhere (arb(float(.))).  Output: out/cells_final.json {key: {...}}.
import sys, os, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from flint import arb, ctx
import d4_mine as D4
ctx.prec = 100
P = json.load(open(os.path.join(HERE, 'out', 'cell_parts.json')))
OUT = os.path.join(HERE, 'out', 'cells_final.json')
done = json.load(open(OUT)) if os.path.exists(OUT) else {}
t0 = time.time(); n = 0
keys = sorted(P, key=lambda s: (float(s.split(',')[1]), int(s.split(',')[0])))
for key in keys:
    if key in done: continue
    part = P[key]
    if part is None or 'error' in part:
        done[key] = None; continue
    u, t, w, x = [arb(float(part[nm])) for nm in ('u', 't', 'w', 'x')]
    L = arb(part['l']); de = arb(part['k'])/10**4
    C, kap2, d5, c1, c5r = [arb(part[nm]) for nm in ('C', 'kap2', 'd5', 'c1', 'c5r')]
    d4 = D4.d4_bound(u, t, x, L)
    den = c1*c1 - d4*c5r
    if not (c1 > 0 and den > arb('0.2')):
        done[key] = {'logN_up': None, 'd4': d4.str(12), 'den': den.str(8)}
    else:
        lg = (part['nt']*d4*d5*kap2*C/(w*den)).log() + 2*x*de*L - (de*L).log()
        done[key] = {'logN_up': repr(float(lg.upper())), 'd4': d4.str(12), 'den': den.str(8)}
    n += 1
    if n % 50 == 0:
        json.dump(done, open(OUT, 'w')); print(n, key, done[key], '%.0fs' % (time.time() - t0), flush=True)
json.dump(done, open(OUT, 'w'))
print('done', n, 'new cells', len(done), 'total', '%.0fs' % (time.time() - t0), flush=True)

