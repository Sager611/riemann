import os
BASE = os.path.abspath(os.environ.get('RAMAUDIT', '.scratch/ramaudit')) + '/'
import sys, os, time, json
R=BASE + 'repo/ramanujans-inequality-HEAD'
sys.path.insert(0, R+'/certify/arb'); sys.path.insert(0, R+'/certify')
os.chdir(R+'/certify/arb')
import scan_arb as SA
from flint import arb
BASE = dict(SA.CELL)
MOD = {}
for f in ['cellmod_A.json', 'cellmod_B.json']:
    MOD.update(json.load(open(BASE + f)))
def setvar(V):
    SA.CELL.clear(); SA.CELL.update(BASE)
    if V == 'orig': return
    for k, v in MOD.items():
        if v.get(V) is not None: SA.CELL[k] = arb(v[V])
# their chosen params
chosen = {}
for fn in ['scan_arb_0.txt', 'scan_arb_1.txt']:
    for line in open(fn):
        p = line.split()
        chosen[(float(p[0]), float(p[1]))] = (int(p[3].split('=')[1]), float(p[4].split('=')[1]))
mode = sys.argv[1]
H = 3.0e12
t0 = time.time()
if mode == 'fixed':
    ivs = sorted(k for k in chosen if 797 <= k[0] < 1000)
    for V in ['orig', 'far2', 'rep']:
        setvar(V); worst = (0, None)
        for (a, b) in ivs:
            c, f = chosen[(a, b)]
            q, parts = SA.q_arb(a, b, H, c, f)
            qv = SA.fup(q)
            if qv > worst[0]: worst = (qv, (a, b, c, f, parts['tail']))
        print(V, 'worst q over [797,1000) with their (c,f):', worst, '%.0fs' % (time.time()-t0), flush=True)
elif mode == 'search':
    V = sys.argv[2]; a = float(sys.argv[3]); b = float(sys.argv[4])
    setvar(V); best = (9, None)
    for c in range(50, 111, 10):
        for f in [0.40 + 0.05*i for i in range(9)]:
            q, parts = SA.q_arb(a, b, H, c, f); qv = SA.fup(q)
            if qv < best[0]: best = (qv, c, f, parts['tail'])
    c0, f0 = best[1], best[2]
    for c in range(max(30, c0-8), min(110, c0+8)+1, 2):
        for f in [f0 - 0.04 + 0.01*i for i in range(9)]:
            if not (0.3 <= f <= 0.98): continue
            q, parts = SA.q_arb(a, b, H, c, f); qv = SA.fup(q)
            if qv < best[0]: best = (qv, c, f, parts['tail'])
    print(V, 'interval', a, b, 'best q', best, '%.0fs' % (time.time()-t0), flush=True)

