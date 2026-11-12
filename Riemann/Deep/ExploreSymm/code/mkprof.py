import numpy as np, sys
d = np.load(sys.argv[1]); xs = d['xs']; lf = d['logV'].astype(float) + d['logP'].astype(float)
h = np.array([0.5*(np.log(np.mean(np.exp(2*(l - l.max())))) + 2*l.max()) for l in lf])
m = lf.mean(axis=1)
np.savez(sys.argv[2], xs=xs, h=h, m=m)
print(np.round(xs, 3)); print(np.round(h, 4)); print(np.round(m, 4))
