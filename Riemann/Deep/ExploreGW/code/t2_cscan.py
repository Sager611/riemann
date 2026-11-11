
import numpy as np
d = np.load("t2_grid.npz"); ts = d["ts"]; Z = d["Z"]; Zq = Z[0]; ZG = Z[1] + Z[2]
m = ts <= 790.0; ts, Zq, ZG = ts[m], Zq[m], ZG[m]
def nsc(z): return int(np.sum(np.signbit(z[1:]) != np.signbit(z[:-1])))
N0 = nsc(Zq)
cs = np.round(np.arange(-1.570, 1.5705, 0.0005), 4)
rows = []
for c in cs:
    z = np.cos(c)*Zq + np.sin(c)*ZG
    dz1 = z[1:-1] - z[:-2]; dz2 = z[1:-1] - z[2:]
    ws = ((dz1 > 0) & (dz2 > 0) & (z[1:-1] < 0)) | ((dz1 < 0) & (dz2 < 0) & (z[1:-1] > 0))
    nws = int(ws.sum()); first = float(ts[1:-1][ws][0]) if nws else None
    rows.append((c, (N0 - nsc(z))/2, nws, first))
clean = [r[0] for r in rows if r[1] == 0 and r[2] == 0]
print("c-grid step 0.0005 on [-1.570, 1.570], T = 790: c with deficit 0 and no wrong-sign extremum: %d values, range [%.4f, %.4f]" % (len(clean), min(clean), max(clean)))
print("all c outside that range have a wrong-sign extremum below 790:", all(r[2] > 0 for r in rows if r[0] < min(clean) or r[0] > max(clean)))
gaps = [r for r in rows if (r[0] < min(clean) or r[0] > max(clean)) and r[2] == 0]
print("exceptions:", gaps[:10])
for c in (-0.0025, -0.003, -0.01, -0.05, -0.1, -0.5, -1.0, 0.0135, 0.02, 0.1, 0.5, 1.0):
    r = min(rows, key=lambda x: abs(x[0]-c)); print("  c=%+.4f deficit %.1f, wrong-sign extrema %d, lowest at %s" % r)
