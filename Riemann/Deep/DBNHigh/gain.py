#!/usr/bin/env python3
# dbnhigh gain.py: endpoint of the theta141 float chain with per-row source gains dS_i interpolated from the LP sweeps.
import json, math, sys
sys.argv = [sys.argv[0], sys.argv[1], '/dev/null']
src = open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), 'chain.py')).read()
src = src[:src.index('t0 = time.time()')]
exec(src)
def interp(pts):
    ks = sorted(pts); out = []
    for i in range(67):
        if i <= ks[0]: out.append(pts[ks[0]]); continue
        if i >= ks[-1]: out.append(pts[ks[-1]]); continue
        for a, b in zip(ks, ks[1:]):
            if a <= i <= b: out.append(pts[a] + (pts[b] - pts[a])*(i - a)/(b - a)); break
    return out
E0 = run()[0]
cases = {
 'A_jets_gain_abs_v_model': {0: 0.0004, 10: 0.0026, 20: 0.0047, 30: 0.0030, 40: 0.0058, 50: 0.0049, 60: 0.0186},
 'B_jets_gain_envelope_model': {0: 0.00034, 20: 0.00043, 40: 0.0018, 60: 0.0069},
 'C_generous_upper': None,
 'D_jets_gain_rich_abs_v': {0: 0.0003, 5: 0.0039, 10: 0.0026, 15: 0.0042, 20: 0.0047, 25: 0.0061, 30: 0.0030, 40: 0.0058, 50: 0.0049, 60: 0.0186, 65: 0.0237},
 'E_twice_D': {0: 0.0006, 5: 0.0078, 10: 0.0052, 15: 0.0084, 20: 0.0094, 25: 0.0122, 30: 0.0060, 40: 0.0116, 50: 0.0098, 60: 0.0372, 65: 0.0474},
}
res = {'base': E0}
for k, pts in cases.items():
    dS = [0.006 if i < 50 else (0.02 if i < 60 else 0.03) for i in range(67)] if pts is None else interp(pts)
    E = run(dS=dS)[0]; res[k] = dict(endpoint=E, gain=E0 - E)
    print(k, repr(E), 'gain %.3e' % (E0 - E), flush=True)
json.dump(res, open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), 'logs', 'gain.json'), 'w'), indent=1)
