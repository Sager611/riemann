"""liuaudit2 T12b: control family of angular twists on the saved data (same 81 classes, X <= 3e5).
For k in -3..3: sums of Lambda^k gamma_2 and of Lambda^k alone; class energy E = mean|S_b|^2/mean n_b, number of 3-sigma classes,
class 0 |S_0|/sqrt(n_0). Under the theory (constant terms only in angular type (0,0)), only k = 0 with gamma_2 should show main terms.
Note: on primes Lambda = -gamma_2^3, so Lambda^k gamma_2 = (-1)^k gamma_2^(3k+1) on primes.
Run: uv run --no-project --with numpy python3 theta_bias_twists.py data.npz
"""
import sys, math
import numpy as np
d = np.load(sys.argv[1])
N = d['N']; b = d['b'].astype(int); g = d['g']; lam = d['lam']
import os
_cap = float(os.environ.get('XCAP', 1e30)); _k = N <= _cap
N, b, g, lam = N[_k], b[_k], g[_k], lam[_k]
Xmax = int(N.max()); nb = 81
def stats(v, X):
    m = N <= X
    S = np.zeros(nb, dtype=complex); np.add.at(S, b[m], v[m]); n = np.bincount(b[m], minlength=nb).astype(float)
    return S, n
print("%-18s %8s %8s %8s %10s %12s" % ("coefficient", "E(3e4)", "E(1e5)", "E(3e5)", "#3sigma", "class0 z"))
for k in range(-3, 4):
    for lab, v in (("Lambda^%+d g2" % k, lam**k*g), ("Lambda^%+d" % k, lam**k)):
        if lab == "Lambda^+0": continue
        Es = []
        for X in (30000, 100000, Xmax):
            S, n = stats(v, X); Es.append(float(np.mean(np.abs(S)**2)/np.mean(n)))
        S, n = stats(v, Xmax); z = np.abs(S)**2/n
        print("%-18s %8.2f %8.2f %8.2f %10d %12.2f" % (lab, Es[0], Es[1], Es[2], int((z > 9).sum()), math.sqrt(z[0])))
