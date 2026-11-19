"""padhi2: shared setup. Imports Padhi's certified scan (certify/arb/scan_arb.py) from the upstream repo copy
(env PADHI_REPO, default .scratch/padhi2/repo/ramanujans-inequality-main) and defines the N3 variants."""
import os, sys, math
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..'))
REPO = os.environ.get('PADHI_REPO', os.path.join(ROOT, '.scratch', 'padhi2', 'repo', 'ramanujans-inequality-main'))
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '1'
sys.path.insert(0, os.path.join(REPO, 'certify', 'arb'))
import scan_arb as S                     # loads Padhi's certified cells (cells_arb_*.json), lambda table, KLN, Bellotti rows
from flint import arb
import json as _json
if os.environ.get('CELLS_OVERRIDE'):           # e.g. the NO_REP cells of dstar0.py
    for _k, _v in _json.load(open(os.environ['CELLS_OVERRIDE'])).items(): S.CELL[_k] = arb(_v['logN_up'])

# Bellotti Table 2 rows 1-13 (T0 < T <= T1 for log T1 = 29..41): printed C and the corrected C of BellottiFix/logs/table.md
C_PRINTED = [37341.72, 37564.18, 40014.59, 42840.01, 46015.02, 49270.71, 52249.66, 55524.80, 59715.72, 64228.93, 68446.65, 71907.32, 76267.80]
C_BOTH = [66414.87, 66849.95, 71240.48, 76305.96, 82030.96, 87885.92, 93268.04, 99202.68, 106734.91, 114908.23, 122531.22, 128782.01, 136683.00]
C_ALL = [76414.26, 75333.28, 75998.39, 79273.09, 84129.86, 89588.55, 94790.68, 100687.03, 108282.25, 116564.75, 124303.64, 130645.33, 138682.14]
UNIF_RATIO = 2.046                       # largest 'all fixes' ratio (row 1); the uniform 1.62e11 bound was not recomputed by BellottiFix
ORIG = dict(ROWS=list(S.ROWS), ROWLOG=list(S.ROWLOG), C1=S.LOG_C1, C2=S.LOG_C2)

def set_variant(name):
    """published | both | all | none.  'both'/'all': Table 2 C replaced; uniform C scaled by UNIF_RATIO; e^29 bound = row 1.
    'none': N3 removed entirely (every Bellotti candidate dropped) - valid for any correction of Bellotti's constants."""
    assert [float(r[3].mid()) for r in ORIG['ROWS']] == C_PRINTED and len(C_BOTH) == len(C_ALL) == 13
    if name == 'published':
        S.ROWS[:] = ORIG['ROWS']; S.ROWLOG[:] = ORIG['ROWLOG']; S.LOG_C1, S.LOG_C2 = ORIG['C1'], ORIG['C2']; return
    if name == 'kln2':                    # sensitivity: Kadiri-Lumley-Ng A and B doubled
        S.KLN[:] = [(s0, A*2, B*2) for s0, A, B in S.KLN]; S._LT.clear(); return
    if name == 'nokln':                   # sensitivity: Kadiri-Lumley-Ng removed
        S.KLN_S = []; S._LT.clear(); return
    if name == 'none':
        S.ROWS[:] = []; S.ROWLOG[:] = []; S.LOG_C1 = S.LOG_C2 = arb(10**6); return
    Cs = {'both': C_BOTH, 'all': C_ALL}[name]
    rows = []
    for (a, b, al, C, B), Cn in zip(ORIG['ROWS'], Cs):
        rows.append((a, b, al, arb(str(Cn)), B))
    S.ROWS[:] = rows; S.ROWLOG[:] = [r[3].log() for r in rows]
    S.LOG_C1 = (arb('1.62e11')*arb(str(UNIF_RATIO))).log(); S.LOG_C2 = arb(str(Cs[0])).log()
