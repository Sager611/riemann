# check5.py -- re-verify a random sample of the certified floors that Planat's producer consumed (reproducibility check).
# For each sampled (cell, height, shift): rebuild every subrange from scratch, re-verify the one-point top Ctop (mu = 0),
# re-verify the paired certificate at the stored c (mu re-optimised, then interval check), recompute the floor.
import sys, json, gzip, random, time, os
from fractions import Fraction as F
import numpy as np
import lib5 as L, cert5 as C
E2 = os.environ.get('DBN5_E2', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
NS = int(sys.argv[1]) if len(sys.argv) > 1 else 48
certs = {(r['cell'], r['h']): r for r in json.load(gzip.open(os.path.join(E2, 'certs_dbn5.json.gz'), 'rt'))}
fl = {(r['cell'], r['h']): r['floor'] for r in json.load(open(os.path.join(E2, 'floors_dbn5.json')))['floors']}
H = json.load(gzip.open(os.path.join(E2, 'classical_dbn5.json.gz'), 'rt'))
usedf = {(int(rec['cell']), F(rec['height']), F(rec['floor'])) for rec in H['records']}
cons = [(cell, h, s, v) for (cell, h), fv in fl.items() for s, v in fv.items() if (cell, F(h), F(v)) in usedf]
random.seed(20261021)
sample = random.sample(cons, min(NS, len(cons)))
print('consumed floors that are ours:', len(cons), '; sample', len(sample), flush=True)
ok = 0
for (cell, h, s, v) in sorted(sample):
    t0 = time.time(); Psi = L.PsiTable(cell); C.clear_cache(); good = True; vals = []
    subs = certs[(cell, h)]['shifts'][s]['subs']
    a = F(s)*F(h)
    R0 = L.Rect(F(h), cell, F(subs[0][0]), None if subs[0][1] == 'None' else F(subs[0][1]))
    Ep = C.euler_pair(R0.rmin, R0.rmax, a)
    for (s1, s2, Ct, c, K, fsub) in subs:
        R = L.Rect(F(h), cell, F(s1), None if s2 == 'None' else F(s2)); G = C.Geometry(R, Psi)
        Sd = {'Clo': np.ones(G.nb), 'Chi': np.ones(G.nb), 'Wc': G.W0, 'Ws': G.W0}
        if C.psi_ia(G, Sd, Ct, Ct, 0.0) < 0: good = False
        if K > 0 and C.verify_cert(G, G.shift_data(a), c, Ct, K) is None: good = False
        fvv = C.floor_value(R, G, c, Ep); vals.append(float(fvv.lower()))
        if float(fvv.lower()) < fsub - 1e-12: good = False
    if min(vals) < float(F(v)): good = False
    ok += good
    print(cell, h, s, 'floor', float(F(v)), 'recomputed min', round(min(vals), 9), 'OK' if good else 'FAIL', round(time.time() - t0, 1), 's', flush=True)
print('re-verified', ok, 'of', len(sample), flush=True)
