# e1arb: float optimum (opt_<nu>.json, entry with m) -> exact rational parameter file for cert.py.
# R rounded to 1e-6, x_j to 1e-12, r (= float lam) to 1e-9, N = 12 (P = odd Taylor polynomial of sinh, degree 23, at r x, over its value at r).
# usage: mkpar.py opt.json m theta nu out.json
import json, sys
from fractions import Fraction as F
opt, m, th, nu, out = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4], sys.argv[5]
d = [e for e in json.load(open(opt)) if e['m'] == m][0]
rnd = lambda v, e: str(F(round(v*10**e), 10**e))
json.dump({'theta': th, 'nu': nu, 'R': rnd(d['R'], 6), 'r': rnd(d['lam'], 9), 'N': 12, 'x': [rnd(x, 12) for x in d['x']],
           'float_kappa': d['kappa'], 'float_nu': d['nu']}, open(out, 'w'))
