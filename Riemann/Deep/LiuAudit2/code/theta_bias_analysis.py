"""liuaudit2 T12 analysis of the saved squarefree-ideal data of theta_bias_q12.py (npz).
Per B_S^(3) class b (81 classes): partial sums S_b(X) of gamma_2, conj(Lambda) gamma_2 (= Liu's a_0), Lambda gamma_2.
Statistics: (1) normalized class energies E(X) = mean_b |S_b(X)|^2 / mean_b n_b(X) (random walk: ~1);
(2) coherence C(X1,X2) = Re sum_b S_b(X1) conj S_b(X2) / sqrt(sum|S_b(X1)|^2 sum|S_b(X2)|^2) against the random-walk value sqrt(n(X1)/n(X2));
(3) the main-term fit S_b(X) ~ c_b X^theta (theta = 5/6), fraction of class energy explained at X_max;
(4) class 0 (canonical: generator a local cube at every place of S) trajectories.
Run: uv run --no-project --with numpy python3 theta_bias_analysis.py data.npz
"""
import sys, math
import numpy as np
d = np.load(sys.argv[1])
N = d['N']; b = d['b'].astype(int); g = d['g']; lam = d['lam']
import os
_cap = float(os.environ.get('XCAP', 1e30)); _k = N <= _cap
N, b, g, lam = N[_k], b[_k], g[_k], lam[_k]
types = {'g2': g, 'Lbar g2 (a_0)': np.conj(lam)*g, 'L g2': lam*g}
Xmax = int(N.max()); order = np.argsort(N, kind='stable'); N = N[order]; b = b[order]
for k in types: types[k] = types[k][order]
CH = np.unique(np.round(np.logspace(3, math.log10(Xmax), 12)).astype(np.int64))
nb = 81
def class_sums(vals, X):
    m = N <= X
    S = np.zeros(nb, dtype=complex); np.add.at(S, b[m], vals[m])
    n = np.bincount(b[m], minlength=nb).astype(float)
    return S, n
print("squarefree exterior ideals: %d, X_max = %d" % (len(N), Xmax))
S0, n0 = class_sums(types['g2'], Xmax)
print("class sizes at X_max: mean %.1f, min %d, max %d, class 0: %d" % (n0.mean(), n0.min(), n0.max(), n0[0]))
print()
print("(1) normalized class energy E(X) = mean_b |S_b|^2 / mean_b n_b  (random walk ~ 1)")
print("%9s " % "X" + " ".join("%16s" % k for k in types))
E = {k: [] for k in types}
for X in CH:
    row = []
    for k, v in types.items():
        S, n = class_sums(v, X); e = float(np.mean(np.abs(S)**2)/np.mean(n)); E[k].append(e); row.append(e)
    print("%9d " % X + " ".join("%16.2f" % e for e in row))
print()
print("(2) coherence between X1 = X_max/10 and X_max (random walk ~ sqrt(n(X1)/n(X2)))")
X1 = Xmax//10
for k, v in types.items():
    Sa, na = class_sums(v, X1); Sb, nb_ = class_sums(v, Xmax)
    c = float(np.real(np.sum(Sa*np.conj(Sb)))/math.sqrt(np.sum(np.abs(Sa)**2)*np.sum(np.abs(Sb)**2)))
    print("  %-16s C = %.3f   random-walk value %.3f" % (k, c, math.sqrt(na.sum()/nb_.sum())))
print()
print("(3) main-term fit S_b(X) = c_b X^(5/6) over the 12 checkpoints, fraction of sum_b sum_X |S_b(X)|^2 explained")
for k, v in types.items():
    M = np.array([class_sums(v, X)[0] for X in CH])          # (12, 81)
    w = CH.astype(float)**(5/6)
    cb = (M*w[:, None]).sum(0)/(w**2).sum()
    resid = M - np.outer(w, cb)
    frac = 1 - float(np.sum(np.abs(resid)**2)/np.sum(np.abs(M)**2))
    print("  %-16s explained fraction %.3f; |c_b| rms %.4f; top classes %s" % (k, frac, float(np.sqrt(np.mean(np.abs(cb)**2))),
          [(int(i), round(float(abs(cb[i])), 4)) for i in np.argsort(-np.abs(cb))[:4]]))
print()
print("(4) class 0 trajectory |S_0(X)| / sqrt(n_0(X)) and arg S_0(X)/(2 pi)")
for k, v in types.items():
    out = []
    for X in CH[::2]:
        S, n = class_sums(v, X); out.append("%d:%.2f/%.2f" % (X, abs(S[0])/math.sqrt(n[0]), np.angle(S[0])/(2*math.pi)))
    print("  %-16s %s" % (k, "  ".join(out)))
print()
print("(5) at X_max: classes with |S_b|^2/n_b > 9 (3-sigma for a random walk):")
for k, v in types.items():
    S, n = class_sums(v, Xmax); z = np.abs(S)**2/n
    print("  %-16s %d classes; max %.1f at class %d; expected count under random walk %.2f" % (k, int((z > 9).sum()), float(z.max()), int(np.argmax(z)), 81*math.exp(-9)))
