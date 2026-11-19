# Compare the Lean certificate lists certB/certB1 (TheoremB.lean) with the JSON etas as rounded by verify_dual.py.
import json, re
from fractions import Fraction as Fr
src = open('Riemann/Deep/ExploreDPP/Lean/TheoremB.lean').read()
def lean_list(name):
    m = re.search(r'def ' + name + r' : List ℤ := \[(.*?)\]', src, re.S)
    return [int(x) for x in m.group(1).replace('\n', ' ').split(',') if x.strip()]
D = 10**12
for name, js in [('certB', 'Riemann/Deep/ExploreDPP/logs/dual_t47_k7X26.json'),
                 ('certB1', 'Riemann/Deep/ExploreDPP/logs/cert_t1k3_X40_lam0.595.json')]:
    d = json.load(open(js))
    L = lean_list(name)
    lam = Fr(str(d['lam'])); k0 = d['k0']; eta = d['eta']
    R = []
    for k, e in enumerate(eta, start=1):
        if k < k0: q = -D
        else:
            q = round(e * D)
            if q < -D: q = -D
        R.append(q)
    print(name, 'keys', sorted(d.keys()), 'theta', d.get('theta'), 'lam', lam, 'k0', k0, 'h', lam / k0,
          'len lean', len(L), 'len json', len(R), 'equal', L == R,
          'min', min(L), 'max', max(L), 'n(-1)', sum(1 for x in L if x == -D), 'last k*h', float(len(L) * lam / k0))
