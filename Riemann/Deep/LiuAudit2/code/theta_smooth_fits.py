"""liuaudit2 T13b: coherent-growth fits for the smoothed full theta sums (theta_smooth_q12.py), with look-elsewhere
correction, and the degree-2-prime diagnosis. Recovers each squarefree c's prime factors by replaying the deterministic
DFS of theta_bias_q12.py on the saved per-prime norms (verified against the saved norms, classes and Lambda values).
Statistics per (k, b): weighted least squares S_{k,b}(D_j) ~ C D_j^theta over 8 disjoint dyadic windows, random-walk
variance sigma_j^2 = sum |b|^2 W^2 (exact); z^2 = |C_hat|^2/Var(C_hat) ~ Exp(1) under the null; family-wise p over 81 classes.
theta = 5/6: Patterson pole (constant-term power |y|^{1/3}); theta = 1/2: coherent lower-order bias (no pole).
Run: uv run --no-project --with numpy python3 theta_smooth_fits.py data.npz
"""
import sys, math
import numpy as np
d = np.load(sys.argv[1])
N = d['N'].astype(np.int64); B = d['b'].astype(int); G = d['g']; L = d['lam']
pq = d['pq'].astype(np.int64); pb = d['pb'].astype(int); plam = d['plam']
XMAX = int(N.max()); nb = 81
# ---- replay DFS (same order as theta_bias_q12.py) to recover factorizations ----
deg2 = np.array([int(round(math.isqrt(int(q))))**2 == q for q in pq])
sys.setrecursionlimit(10000)
fac = [()]
def dfs(start, Nn, members):
    for j in range(start, len(pq)):
        q = int(pq[j])
        if Nn*q > XMAX: break
        fac.append(tuple(members + [j]))
        if j + 1 < len(pq) and Nn*q*int(pq[j+1]) <= XMAX:
            dfs(j+1, Nn*q, members + [j])
dfs(0, 1, [])
assert len(fac) == len(N), (len(fac), len(N))
ok = all(int(np.prod([pq[i] for i in f])) == N[t] if f else N[t] == 1 for t, f in enumerate(fac[:5000]))
lam_ok = max(abs(np.prod([plam[i] for i in f]) - L[t]) if f else abs(L[t]-1) for t, f in enumerate(fac[:5000]))
b_ok = all(((sum(np.array([[(pb[i] // 3**s) % 3 for s in range(4)] for i in f]).reshape(-1, 4)) % 3) if f else np.zeros(4, int)).dot([1, 3, 9, 27]) == B[t] for t, f in enumerate(fac[:3000]))
has2 = np.array([any(deg2[i] for i in f) for f in fac])
print("replayed %d squarefree ideals; norms match: %s; Lambda match err %.1e; classes match: %s; with a degree-2 prime factor: %d" % (
      len(fac), ok, lam_ok, b_ok, int(has2.sum())))
import os
XCAP = int(float(os.environ.get('XCAP', XMAX)))
if XCAP < XMAX:
    keep = N <= XCAP
    N, B, G, L, has2 = N[keep], B[keep], G[keep], L[keep], has2[keep]
    XMAX = XCAP
    print("capped at X = %d: %d ideals" % (XCAP, len(N)))
def W(y):
    out = np.zeros_like(y, dtype=float); m = (y > 0.5) & (y < 1.0); t = (y[m] - 0.75)/0.25
    out[m] = np.exp(-1.0/(1.0 - t*t)); return out
ncube = [(1.0, 1.0+0j, False)] + [(float(pq[i]), plam[i], bool(deg2[i])) for i in range(len(pq)) if pq[i]**3 <= XMAX]
NW = int(sys.argv[2]) if len(sys.argv) > 2 else 8
Ds = XMAX/2.0**np.arange(0, NW)[::-1]
ks = [-3, -2, -1, 0, 1, 2, 3]
def sums(k, mask_c=None, drop_deg2_cube=False):
    S = np.zeros((len(Ds), nb), dtype=complex); V = np.zeros((len(Ds), nb))
    base = G*L**k
    sel = np.ones(len(N), bool) if mask_c is None else mask_c
    for (qn, ln, d2) in ncube:
        if drop_deg2_cube and d2: continue
        coef = math.sqrt(qn)*ln**(3*k)
        for j, D in enumerate(Ds):
            w = W(N*qn**3/D); m = (w > 0) & sel
            np.add.at(S[j], B[m], coef*base[m]*w[m]); np.add.at(V[j], B[m], qn*w[m]**2)
    return S, V
def fit(S, V, theta):
    w = Ds**theta
    ok = V > 0
    num = np.where(ok, S*(w[:, None]/np.where(ok, V, 1)), 0).sum(0)
    den = np.where(ok, (w[:, None]**2)/np.where(ok, V, 1), 0).sum(0)
    C = num/den; z2 = np.abs(C)**2*den
    return C, z2
def report(tag, mask_c=None, drop=False):
    print(tag)
    print("   k   R(Dmax) | theta=5/6: max z2 (class)  #z2>9  mean z2  p_fw     | theta=1/2: max z2 (class)  #z2>9  mean z2  p_fw")
    for k in ks:
        S, V = sums(k, mask_c, drop)
        R = float(np.sum(np.abs(S[-1])**2)/np.sum(V[-1]))
        out = []
        for th in (5/6, 0.5):
            C, z2 = fit(S, V, th); zm = float(z2.max())
            out.append("%8.2f (%2d)  %5d  %7.2f  %.2g" % (zm, int(np.argmax(z2)), int((z2 > 9).sum()), float(z2.mean()), 1 - (1 - math.exp(-zm))**81))
        print("  %+d  %7.3f | %s | %s" % (k, R, out[0], out[1]))
report("(A) all coefficients (Liu's full sum: c squarefree, n any; 81 classes; %d dyadic windows %d..%d)" % (NW, int(Ds[0]), XMAX))
report("(B) without degree-2 primes (c with a degree-2 factor and cube parts n of degree 2 removed)", ~has2, True)
print()
print("(C) class 0, k = -1, per window (all coefficients): |S|/sigma, arg/2pi; same with degree-2 terms removed")
for tag, msk, dr in (("all", None, False), ("no deg-2", ~has2, True)):
    S, V = sums(-1, msk, dr)
    print("   %-9s " % tag + "  ".join("%d:%.2f/%.2f" % (D, abs(S[j, 0])/math.sqrt(V[j, 0]) if V[j, 0] > 0 else float('nan'), np.angle(S[j, 0])/(2*math.pi)) for j, D in enumerate(Ds)))
S1, V1 = sums(-1, has2, False)
print("   only c with a degree-2 factor: " + "  ".join("%d:%.2f/%.2f" % (D, abs(S1[j, 0])/math.sqrt(V1[j, 0]) if V1[j, 0] > 0 else float('nan'), np.angle(S1[j, 0])/(2*math.pi)) for j, D in enumerate(Ds)))

print()
print("(D) growth and phase per (k, b): slope of log(|S|/sigma) against log D (pole 5/6: +1/3; random walk or O(sqrt D) bias: 0),")
print("    phase concentration Rbar = |mean_j S_j/|S_j|| (fixed phase: 1; random phases: about 1/sqrt(n_windows)); empirical null = the 5 other nonzero twists")
nw = len(Ds); lD = np.log(Ds)
rows = {}
for k in ks:
    S, V = sums(k)
    ok = (V > 0).all(0)
    rat = np.where(V > 0, np.abs(S)/np.sqrt(np.where(V > 0, V, 1)), np.nan)
    slope = np.array([np.polyfit(lD, np.log(rat[:, b]), 1)[0] if ok[b] else np.nan for b in range(nb)])
    Rbar = np.abs(np.nanmean(np.where(np.abs(S) > 0, S/np.where(np.abs(S) > 0, np.abs(S), 1), np.nan), axis=0))
    C, z2 = fit(S, V, 5/6)
    rows[k] = (slope, Rbar, z2)
    print("  k=%+d: mean slope %.3f (sd %.3f); class 0 slope %.3f, Rbar %.2f, z2 %.2f; max Rbar %.2f (class %d); classes with Rbar > 0.8: %d" % (
        k, float(np.nanmean(slope)), float(np.nanstd(slope)), float(slope[0]), float(Rbar[0]), float(z2[0]), float(np.nanmax(Rbar)), int(np.nanargmax(Rbar)), int((Rbar > 0.8).sum())))
null_z = np.concatenate([rows[k][2] for k in ks if k not in (0, -1)])
null_R = np.concatenate([rows[k][1] for k in ks if k not in (0, -1)])
for lab, arr, val in (("z2(5/6), class 0, k=-1", null_z, rows[-1][2][0]), ("Rbar, class 0, k=-1", null_R, rows[-1][1][0]),
                      ("max_b z2(5/6), k=-1", np.array([rows[k][2].max() for k in ks if k not in (0, -1)]), rows[-1][2].max())):
    print("  empirical null for %s: value %.3f; fraction of null class statistics >= value: %.4f (n = %d)" % (lab, val, float(np.mean(arr >= val)), len(arr)))
