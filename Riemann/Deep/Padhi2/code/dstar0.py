"""padhi2 sensitivity: recompute Padhi's certified Lemma D cells (certify/arb/cells_arb.cell, his stored parameters) with the
two-zero repulsion switched off (NO_REP=1 -> d* = 0 in Lemma pair), for the heights lk given. Output: logs/<out>.json in
the cells_arb format. usage: dstar0.py <lk,lk,...> <out-name> [kmax]  (cells in increasing k; k > kmax keep the published value)"""
import os, sys, json, pickle, glob, time
SCALE = os.environ.get('DSTAR_SCALE')          # if set: d* multiplied by SCALE instead of switched off
if SCALE is None: os.environ['NO_REP'] = '1'
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '1'
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..'))
REPO = os.environ.get('PADHI_REPO', os.path.join(ROOT, '.scratch', 'padhi2', 'repo', 'ramanujans-inequality-main'))
sys.path.insert(0, os.path.join(REPO, 'certify', 'arb'))
import cells_arb as CA
if SCALE is not None:
    _orig = CA.D._dstar
    CA.D._dstar = lambda L, de: _orig(L, de)*float(SCALE)
lks = [float(x) for x in sys.argv[1].split(',')]
KMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 300
outf = os.path.join(os.path.dirname(__file__), '..', 'logs', sys.argv[2] + '.json')
tab = {}
for f in sorted(glob.glob(os.path.join(CA.CODE, 'lemmaD_tableGH0.01_*.pkl'))): tab.update(pickle.load(open(f, 'rb')))
pub = {}
for f in glob.glob(os.path.join(CA.HERE, 'cells_arb_*.json')): pub.update(json.load(open(f)))
out = {}; t0 = time.time(); mx = None
for key in sorted(k for k in tab if tab[k] is not None):
    k, lk = key[0], key[1]
    if float(lk) not in lks or k > KMAX: continue
    fv, p = tab[key]
    if fv > 1e8: continue
    nm = '%d,%s' % (k, lk)
    lo, val = CA.cell(k, lk, p)
    v = CA.fup(lo); out[nm] = {'logN_up': repr(v)}
    d = v - float(pub[nm]['logN_up'])
    if mx is None or d > mx[0]: mx = (d, nm)
    if len(out) % 20 == 0:
        json.dump(out, open(outf, 'w')); print(len(out), 'cells %.0fs, max increase of log N+ so far %.4f at %s' % (time.time() - t0, mx[0], mx[1]), flush=True)
json.dump(out, open(outf, 'w'))
print('done', len(out), 'cells %.0fs; max increase of log N+ vs published %.4f at %s' % (time.time() - t0, mx[0], mx[1]))
