import json, numpy as np, sys
sys.argv = ['x']
import importlib.util
spec = importlib.util.spec_from_file_location('rm', 'rmodel_lib.py'); rm = importlib.util.module_from_spec(spec); spec.loader.exec_module(rm)
g = json.load(open('gideal.json')); ts = np.array(g['ts']); hs = np.array(g['hs']); G = np.array(g['G'])
def Gi(t, h):
    # bilinear interpolation in t (G decreasing in t: use the grid value at the lower t node = conservative) and log-linear in h
    i = max(0, np.searchsorted(ts, t, side='right') - 1)
    row = np.log(G[i]); return float(np.exp(np.interp(h, hs, row)))
for add in [0.0, 0.1]:
    rm.Gmodel = lambda t, h, fac, add2, add=add: Gi(t, h) + add
    print('derivative, G = Gideal + %.2f:' % add, rm.run(True), flush=True)
