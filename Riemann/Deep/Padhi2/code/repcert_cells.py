"""padhi2 follow-up: rerun Padhi's Thm 5.2 certification (code/rep/rep_cert.certify, ball arithmetic, his code unchanged)
at chosen cells (lk, k) with the d* stored in dstar2c_all.json. usage: repcert_cells.py lk,k [lk,k ...] [--depth N]"""
import os, sys, json, time
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[v] = '1'
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..'))
REPO = os.environ.get('PADHI_REPO', os.path.join(ROOT, '.scratch', 'padhi2', 'repo', 'ramanujans-inequality-main'))
REP = os.path.join(REPO, 'code', 'rep'); os.chdir(REP); sys.path.insert(0, REP)
import rep_cert as RC
RC.lamt = {float(k): v for k, v in json.load(open('lamstar2c.json')).items()}
DA = json.load(open('dstar2c_all.json')); DC = json.load(open('dcert_29_30_31_32_33_34_35_36_37.json'))
args = [a for a in sys.argv[1:] if not a.startswith('--')]
depth = int(sys.argv[sys.argv.index('--depth') + 1]) if '--depth' in sys.argv else 40
scale = float(sys.argv[sys.argv.index('--scale') + 1]) if '--scale' in sys.argv else 1.0
for a in args:
    if ',' not in a: continue
    lk, k = [int(x) for x in a.split(',')]
    d = DA['%d,%d' % (lk, k)]*scale; t0 = time.time()
    ok, nb = RC.certify(lk, k*1e-4*lk, d, maxdepth=depth)
    print('cell (%d,%d): d* = %.6f (directly certified in dcert: %s) -> certify ok=%s boxes=%d  %.0fs'
          % (lk, k, d, ('%d,%d' % (lk, k)) in DC, ok, nb, time.time() - t0), flush=True)
