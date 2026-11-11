
import numpy as np, json
d = np.load("t2_grid.npz"); ts = d["ts"]; Z = d["Z"]; Zq = Z[0]; ZG = Z[1] + Z[2]
# zeros of Zq near 604 (linear interpolation)
i0, i1 = np.searchsorted(ts, 600), np.searchsorted(ts, 608)
zs = []
for i in range(i0, i1):
    if Zq[i] * Zq[i+1] < 0:
        zs.append(ts[i] - Zq[i]*(ts[i+1]-ts[i])/(Zq[i+1]-Zq[i]))
zs = np.array(zs)
print("zeros of Z_12 in [600,608]:", np.round(zs, 4).tolist())
print("gaps:", np.round(np.diff(zs), 4).tolist(), " mean spacing 2pi/log(13 t/2pi) at 604 =", round(2*np.pi/np.log(13*604/(2*np.pi)), 4))
k = np.searchsorted(ts, 604.0)
w = slice(k-40, k+40)
j = np.argmin(np.abs(Zq[w])) 
print("max |Z_12| between the close pair and Z_A there:", float(np.max(np.abs(Zq[np.searchsorted(ts, zs[np.argmin(np.diff(zs))]) : np.searchsorted(ts, zs[np.argmin(np.diff(zs))+1])]))), float(ZG[k]))
crit = json.load(open("t2_crit.json")); a = np.array(sorted(x[1] for x in crit["negmin"]))
print("negative minima |r*| sorted (first 15):", np.round(a[:15], 4).tolist())
for eps in (0.02, 0.05, 0.1, 0.2, 0.5, 1.0):
    print("  #(|r*| < %.2f) = %d" % (eps, int(np.sum(a < eps))))
# Lehmer-type check: normalized gaps of Z_12 zeros below 800, smallest ones and their heights
zz = []
for i in range(len(ts)-1):
    if Zq[i]*Zq[i+1] < 0: zz.append(ts[i] - Zq[i]*(ts[i+1]-ts[i])/(Zq[i+1]-Zq[i]))
zz = np.array(zz); g = np.diff(zz); mid = (zz[1:]+zz[:-1])/2
gn = g * np.log(13*mid/(2*np.pi))/(2*np.pi)
o = np.argsort(gn)[:6]
print("smallest normalized gaps of L(s,chi_13) zeros below 800 (height, gap/mean):", [(round(float(mid[i]),3), round(float(gn[i]),4)) for i in o])
