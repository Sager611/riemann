import os
BASE = os.path.abspath(os.environ.get('RAMAUDIT', '.scratch/ramaudit')) + '/'
import sys, os, time, json
R=BASE + 'repo/ramanujans-inequality-HEAD'
sys.path.insert(0, R+'/certify/arb'); sys.path.insert(0, R+'/certify')
os.chdir(R+'/certify/arb')
import scan_arb as SA
from flint import arb
V = sys.argv[1]; part, nparts = int(sys.argv[2]), int(sys.argv[3]); noN3 = (sys.argv[4] == 'noN3')
MOD = {}
for f in ['cellmod_A.json', 'cellmod_B.json', 'cellmod_C.json', 'cellmod_D.json']:
    MOD.update(json.load(open(BASE + f)))
missing = [k for k in SA.CELL if k not in MOD]
print('cells in scan:', len(SA.CELL), 'recomputed:', len(MOD), 'missing:', len(missing), flush=True)
worse = [(k, float(v[V]) - float(v['old'])) for k, v in MOD.items() if v.get(V) is not None and float(v[V]) > float(v['old'])]
print('cells where', V, '> old:', len(worse), sorted(worse, key=lambda z: -z[1])[:5], flush=True)
for k, v in MOD.items():
    if v.get(V) is not None: SA.CELL[k] = arb(v[V])
    else: SA.CELL.pop(k, None)        # inadmissible after correction: drop the cell bound
if noN3:
    SA.ROWS = []; SA.LOG_C1 = arb(10)**9; SA.LOG_C2 = arb(10)**9
ivs = []
for fn in ['scan_arb_0.txt', 'scan_arb_1.txt']:
    for line in open(fn):
        p = line.split()
        ivs.append((float(p[0]), float(p[1]), int(p[3].split('=')[1]), float(p[4].split('=')[1]), float(p[2])))
ivs.sort()
out = BASE + 'scan_%s_%s_%d.txt' % (V, sys.argv[4], part)
t0 = time.time(); worst = (0, None); nworse = 0
with open(out, 'w') as fo:
    for i, (a, b, c, f, q0) in enumerate(ivs):
        if i % nparts != part: continue
        q, parts = SA.q_arb(a, b, 3.0e12, c, f); qv = SA.fup_safe(q)
        fo.write('%r %r %r  c=%d f=%.2f\n' % (a, b, qv, c, f)); fo.flush()
        if a >= 57 and qv > worst[0]: worst = (qv, a, b)
        if qv > q0 + 1e-12: nworse += 1
print('done', V, sys.argv[4], 'part', part, 'worst q (L1>=57):', worst, 'intervals with q > published:', nworse, '%.0fs' % (time.time() - t0), flush=True)

