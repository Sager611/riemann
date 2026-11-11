
import numpy as np
d = np.load("t2_grid.npz"); ts = d["ts"]; Z = d["Z"]; Zq = Z[0]; ZG = Z[1] + Z[2]
r = Zq/ZG; hh = ts[1]-ts[0]
pos = []
for i in range(1, len(ts)-1):
    if ZG[i-1]*ZG[i] <= 0 or ZG[i]*ZG[i+1] <= 0: continue
    a, b, c_ = r[i-1], r[i], r[i+1]
    if b > a and b > c_ and b > 0:
        den = a - 2*b + c_; x = 0.5*(a-c_)/den; pos.append((ts[i]+x*hh, b - 0.25*(a-c_)*x))
pos.sort(); best = 1e9; rec = []
for t, v in pos:
    if v < best: best = v; rec.append((round(float(t),3), round(float(v),5), round(float(np.arctan(v)),5)))
print("c<0 births (positive local maxima of r): %d; records (t*, r*, |c*|):" % len(pos)); [print("  ", x) for x in rec]
