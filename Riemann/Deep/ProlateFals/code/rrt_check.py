# compare with the certified eigenvalues of Ramis, Richard-Jung, Thomann, C. R. Math. 363 (2025) 1065-1081, sec. 4.2.3 (Lambda = sqrt2)
from decimal import Decimal as D, getcontext
import sys
getcontext().prec = 30
ev = {int(l.split()[0]): D(l.split()[1]) for l in open(sys.argv[1]) if not l.startswith('#')}
rrt = {2: D('-39.3832165744668224'), 148: D('-9104.3331714128495040')}
for lab, mu in rrt.items():
    k = lab//2; E = ev[k]; m = -(E*E)
    print(f'RRT mu_-{lab} = {mu}; ours (even #{k}) -E^2 = {m:.16f}; diff {m-mu:.2e}, relative {(m-mu)/mu:.1e}')
