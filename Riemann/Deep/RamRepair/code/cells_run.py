# cells_run.py -- compute cell_parts for every stored parameter cell (Padhi's tables only supply (u, tau, w, x)); output parts JSON
import sys, os, json, time, glob, pickle
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import cells_mine as CM
R = os.path.join(HERE, 'repo', 'ramanujans-inequality-main', 'code')
def load_params():
    tab = {}
    for f in sorted(glob.glob(os.path.join(R, 'lemmaD_tableGH0.01_*.pkl'))): tab.update(pickle.load(open(f, 'rb')))
    out = {}
    for key, v in tab.items():
        if v is None: continue
        fv, p = v
        if fv > 1e8: continue
        out['%d,%s' % (key[0], key[1])] = (key[0], float(key[1]), [float(x) for x in p])
    return out
def job(item):
    nm, (k, l, p) = item
    try: return nm, CM.cell_parts(k, l, p)
    except Exception as e: return nm, {'error': repr(e)}
if __name__ == "__main__":
    outf = sys.argv[1] if len(sys.argv) > 1 else 'out/cell_parts.json'
    P = load_params()
    try: done = json.load(open(outf))
    except Exception: done = {}
    todo = sorted((nm, v) for nm, v in P.items() if nm not in done)
    print('cells', len(P), 'todo', len(todo), flush=True)
    t0 = time.time(); n = 0
    with Pool(3) as pool:
        for nm, r in pool.imap_unordered(job, todo, chunksize=4):
            done[nm] = r; n += 1
            if n % 50 == 0:
                json.dump(done, open(outf + '.tmp', 'w')); os.replace(outf + '.tmp', outf)
                print(n, nm, '%.0fs' % (time.time() - t0), flush=True)
    json.dump(done, open(outf + '.tmp', 'w')); os.replace(outf + '.tmp', outf)
    print('done', n, '%.0fs' % (time.time() - t0), flush=True)
