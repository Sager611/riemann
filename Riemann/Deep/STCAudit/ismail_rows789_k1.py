
# Ismail rows 7-9 (Patel-Yang sub-Weyl input) re-evaluated with k1 = 66.7 (as published), 68.54 (what Patel-Yang's own
# parameters give, REFUTATIONS.md) and 54.6 (our sub-Weyl constant, discoveries/sub-weyl-constant.md); Ismail's sbound.py.
import sys, json
sys.path.insert(0, '.scratch/stcaudit/ism')
import mpmath as mp
from sbound import Params, c1_of_Q0
ROWS = [(1e40, 1.0660270, 1.215723, 0.066026), (1e60, 1.0677880, 1.135577, 0.067787), (1e100, 1.0428930, 1.085787, 0.042892)]
Q0 = 1e9; out = []
for k1 in ['66.7', '68.54', '54.6']:
    for Topt, c, r, e in ROWS:
        P = Params(c, r, e, n=4, k=(mp.mpf(k1), mp.mpf(27)/164, 0), c1=c1_of_Q0(Q0), c2=1, Q0=Q0, Q1=1.5, Q2=1.5, Q3=Q0, Q10=Q0, Q11=Q0, dps=30)
        o = P.constants()
        rec = dict(k1=k1, Topt=Topt, admissible=bool(P.admissible()), C1=float(o['C1']), C2=float(o['C2']), C3T=float(o['C3T']), C2p=float(o['C2p']), C3Tp=float(o['C3Tp']))
        out.append(rec); print(rec); sys.stdout.flush()
json.dump(out, open('Riemann/Deep/STCAudit/ismail_rows789_k1.json', 'w'), indent=1)
