"""liuaudit2 T13: Liu's own smoothed theta sums (Lemma 5.1, l.2630-2640; Prop 5.2, l.2899-2903) over F = Q(zeta12),
with the full coefficient set a(c n^3) = gamma_2(c) vartheta(c n^3) / (sqrt(q_c) q_n), c squarefree, n arbitrary (no
squarefree restriction on the ideal), vartheta = Lambda^k (k = -1 is Liu's conj(Lambda); with eta = G_0 his
vartheta = conj(Lambda) conj(G_0) eta equals conj(Lambda); characters of B_S^(3) are in his finite group Upsilon_S).
Normalization b = q^{1/2} a: b(c n^3) = q_n^{1/2} gamma_2(c) Lambda(c)^k Lambda(n)^{3k}.
Smoothed sums S_{k,b}(D) = sum_{a in class b} b(a) W(q_a/D), W a C_c^infty bump on (1/2,1) (Liu: W in C_c^infty(I)).
Predictions: Liu (entire D_m, l.2892, Prop 5.2 with rho_nu = |nu|_S >= 1): k = -1 sums are O(D^{1/2-A}), i.e. far below the
random-walk size sigma = (sum |b|^2 W^2)^{1/2}. Patterson pole (k = 0, constant term |y|^{1/3}): S = C_b W^(5/6) D^{5/6} + O(D^{-A}).
Input: npz of theta_bias_q12.py (squarefree c <= 3e5 with gamma_2, Lambda, class; per-prime arrays).
Run: uv run --no-project --with numpy python3 theta_smooth_q12.py data.npz
"""
import sys, math
import numpy as np
d = np.load(sys.argv[1])
N = d['N'].astype(float); B = d['b'].astype(int); G = d['g']; L = d['lam']
pq = d['pq'].astype(float); plam = d['plam']
Xmax = N.max(); nb = 81
def W(y):
    out = np.zeros_like(y)
    m = (y > 0.5) & (y < 1.0)
    t = (y[m] - 0.75)/0.25
    out[m] = np.exp(-1.0/(1.0 - t*t))
    return out
# Mellin transform W^(s) = int W(y) y^(s-1) dy at s = 5/6 (main-term weight)
yy = np.linspace(0.5, 1.0, 200001); fy = W(yy)*yy**(5/6 - 1); Wh = float(np.sum((fy[1:] + fy[:-1])*np.diff(yy))/2)
# cube cofactors n: all exterior ideals with q_n^3 <= Xmax (here primes only; products exceed the bound)
ncube = [(1.0, 1.0+0j)] + [(pq[i], plam[i]) for i in range(len(pq)) if pq[i]**3 <= Xmax]
print("cube cofactors n with q_n^3 <= %.0f: %d (norms %s)" % (Xmax, len(ncube), sorted(set(int(q) for q, _ in ncube))))
Ds = Xmax/2.0**np.arange(0, 9)[::-1]          # 9 scales, factor 2 apart: windows [D/2, D] disjoint
Ds = Ds[Ds >= 1000]
ks = [-3, -2, -1, 0, 1, 2, 3]
res = {}
for k in ks:
    S = np.zeros((len(Ds), nb), dtype=complex); V = np.zeros((len(Ds), nb))
    base = G*L**k
    for (qn, ln) in ncube:
        coef = math.sqrt(qn)*ln**(3*k)
        for j, D in enumerate(Ds):
            w = W(N*qn**3/D)
            m = w > 0
            np.add.at(S[j], B[m], coef*base[m]*w[m])
            np.add.at(V[j], B[m], (qn)*w[m]**2)
    res[k] = (S, V)
print("W^(5/6) = %.6f; scales D = %s" % (Wh, [int(x) for x in Ds]))
print()
print("(1) energy ratio R_k(D) = sum_b |S_{k,b}(D)|^2 / sum_b sigma_{k,b}(D)^2  (random walk ~ 1; Liu k=-1: << 1; pole: grows like D^{2/3})")
print("%9s " % "D" + " ".join("%9s" % ("k=%+d" % k) for k in ks))
for j, D in enumerate(Ds):
    print("%9d " % D + " ".join("%9.3f" % (np.sum(np.abs(res[k][0][j])**2)/np.sum(res[k][1][j])) for k in ks))
print()
print("(2) k = 0 (positive control): fit S_{0,b}(D) = C_b W^(5/6) D^(5/6) over the scales; relative residual energy")
S0 = res[0][0]; w56 = Wh*Ds**(5/6)
Cb = (S0*w56[:, None]).sum(0)/(w56**2).sum()
resid = S0 - np.outer(w56, Cb)
print("    fraction of energy explained by C_b D^(5/6): %.4f ; residual energy / random-walk energy: %.3f" % (
      1 - np.sum(np.abs(resid)**2)/np.sum(np.abs(S0)**2), np.sum(np.abs(resid)**2)/np.sum(res[0][1])))
for th in (0.5, 2/3, 5/6, 1.0):
    w = Wh*Ds**th; C = (S0*w[:, None]).sum(0)/(w**2).sum(); r_ = S0 - np.outer(w, C)
    print("    theta = %.4f: explained %.4f" % (th, 1 - np.sum(np.abs(r_)**2)/np.sum(np.abs(S0)**2)))
print("    |C_b|: rms %.4f, max %.4f (class %d); total sum over classes: |sum_b C_b| = %.4f" % (
      float(np.sqrt(np.mean(np.abs(Cb)**2))), float(np.abs(Cb).max()), int(np.argmax(np.abs(Cb))), float(abs(Cb.sum()))))
print()
print("(3) Patterson-type coherent fit for every k: C_{k,b} from S_{k,b}(D) ~ C W^ D^(5/6), z_{k,b}^2 = |C_hat|^2 / Var(C_hat) under the")
print("    random-walk null (Var from sigma^2; complex Gaussian: P(z^2 > t) = e^{-t}); look-elsewhere over 81 classes: p_fw = 1-(1-e^{-zmax^2})^81")
for k in ks:
    S, V = res[k]
    w = Wh*Ds**(5/6)
    num = (S*(w[:, None]/V)).sum(0); den = ((w[:, None]**2)/V).sum(0)
    Chat = num/den; var = 1.0/den
    z2 = np.abs(Chat)**2/var
    zmax = float(z2.max()); arg = int(np.argmax(z2))
    pfw = 1 - (1 - math.exp(-zmax))**81
    print("    k=%+d: max z^2 = %7.2f (class %2d), classes with z^2 > 9: %2d, mean z^2 = %6.2f, family-wise p = %.3g" % (
          k, zmax, arg, int((z2 > 9).sum()), float(z2.mean()), pfw))
print()
print("(4) class 0 and class 9 trajectories, k = -1 (Liu's a_0) and k = 0:  |S|/sigma and arg(S)/2pi per scale")
for k in (-1, 0, 1):
    for b in (0, 9):
        S, V = res[k]
        print("    k=%+d b=%d: " % (k, b) + "  ".join("%d:%.2f/%.2f" % (D, abs(S[j, b])/math.sqrt(V[j, b]), np.angle(S[j, b])/(2*math.pi)) for j, D in enumerate(Ds)))
