# prolcount.py: M_low(t, tau) = #{even m : delta_m(t) <= tau} (exact rungs, ACLGFals prolate.py) at t^2 = c and c*sigma^2
import sys, json, math
sys.path.insert(0, '/Users/adriansager/playground/riemann/Riemann/Deep/ACLGFals/code')
import mpmath as mp
from prolate import rungs
out = []
for c in [float(x) for x in sys.argv[1].split(',')]:
    nmax = int(4*c) + 40
    r, dps, K = rungs(mp.mpf(c), nmax, 80)
    d = {n: float(v['delta']) for n, v in r.items()}
    row = {'c': c, 'nmax': nmax, 'K': K}
    for tau in (1e-2, 1e-5, 1e-10, 1e-20):
        row['le_%g' % tau] = sum(1 for n, v in d.items() if v <= tau)
        row['le_%g_mod0' % tau] = sum(1 for n, v in d.items() if v <= tau and n % 4 == 0)
    row['max_n_used'] = max(d)
    row['delta_last'] = d[max(d)]
    out.append(row); print(json.dumps(row), flush=True)
