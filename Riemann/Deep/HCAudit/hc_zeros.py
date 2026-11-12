# hcaudit falsification on actual zeta zeros: Lamzouri's E_lambda (l.405, fixed L = log T), Z_0, and every step of (iii')
# with the cert768 f (EPS = 1/1000), A, beta and c_cert. Usage: EPS=1/1000 python3 hc_zeros.py cert768.json zeta_zeros.npy -0.3478616532
import sys, json, os, math
from fractions import Fraction as Fr
from flint import fmpq_poly, fmpq, arb, arb_poly, ctx
import numpy as np
ctx.prec = 128
d = json.load(open(sys.argv[1])); g_all = np.sort(np.load(sys.argv[2]).astype(float)); c = float(sys.argv[3])
lam = Fr(d['lam']); D = d['D']; s = Fr(d['s']); U = Fr(d['U']); dg = Fr(d['dg'])
B = [Fr(b) for b in d['B']]; B[0] += Fr(os.environ.get('EPS', '0')); A = [float(Fr(a)) for a in d['A']]; beta = float(Fr(d['beta']))
q2 = lambda x: fmpq(x.numerator, x.denominator)
H = [fmpq_poly([1]), fmpq_poly([0, 2])]
for n in range(1, 2*D): H.append(fmpq_poly([0, 2])*H[n] - 2*n*H[n-1])
Q = sum((q2(B[n]*(-1)**n)*H[2*n] for n in range(D+1)), fmpq_poly([0]))
QA = arb_poly([arb(Q[i]) for i in range(Q.degree()+1)]); sA = arb(q2(s)); k_ = (2*arb.pi()).sqrt()*sA
def fh(u): v = k_*arb(float(u)); return float((sA*QA(v)*(-(v*v)/2).exp()).mid())
lamf, Uf, dgf = float(lam), float(U), float(dg); nodes = np.array([lamf + k*dgf for k in range(len(A))])
Af = lambda x: float(np.interp(x, nodes, A, right=A[-1]))
fh0 = fh(0.0)
print('zeros: %d, gamma in [%.3f, %.3f], sorted, min diff %.3e' % (len(g_all), g_all[0], g_all[-1], np.diff(g_all).min()))
# T between the last two zeros so that every zero <= T has its right neighbour in the file
T = 0.5*(g_all[-2] + g_all[-1]); L = math.log(T); Z = g_all[g_all <= T]; NT = len(Z)
thr = 2*math.pi*lamf/L
gap = np.diff(g_all)[:NT]              # gap[n] = gamma_{n+1} - gamma_n (0-based n), right neighbour exists
left = np.r_[np.inf, gap[:NT-1]]       # gamma_n - gamma_{n-1}
inE = (np.minimum(gap, left) <= thr); inE[0] = False   # n >= 2 in Lamzouri's 1-based indexing
Z0 = Z[~inE]; N0 = len(Z0); x = Z0*L/(2*math.pi); g = np.diff(x)
print('T = %.4f, L = %.5f, N(T) = %d (RvM main term %.1f), TL/2pi = %.1f, #E = %d (%.4f N), N0 = %d' % (T, L, NT,
      T/(2*math.pi)*math.log(T/(2*math.pi)) - T/(2*math.pi) + 7/8, T*L/(2*math.pi), inE.sum(), inE.sum()/NT, N0))
print('min Z_0 gap (unfolded, fixed L) = %.6f > lambda = %.3f : %s' % (g.min(), lamf, g.min() > lamf))
# (1) split of S(f) by index distance, with Montgomery weight w and with w = 1
K = 1
while (x[K:] - x[:-K]).min() < 14: K += 1
S = {}; S1 = {}; maxpos3 = -1e9; maxw = 0.0
for k in range(1, K+1):
    u = x[k:] - x[:-k]; vals = np.array([fh(t) for t in u]); dgam = (Z0[k:] - Z0[:-k]); w = 4/(4 + dgam**2)
    S[k] = 2*(vals*w).sum(); S1[k] = 2*vals.sum()
    if k >= 3: maxpos3 = max(maxpos3, vals.max()); assert u.min() > 3*lamf
    if k <= 2: maxw = max(maxw, np.abs(vals*(w - 1)).max())
Sw = N0*fh0 + sum(S.values()); Sdrop = N0*fh0 + S[1] + S[2]
print('K = %d index distances cover |u| < 14; S_w = %.4f, S_w/N(T) = %.5f, S_w/N0 = %.5f' % (K, Sw, Sw/NT, Sw/N0))
print('(1) k>=3 terms: max f^ = %.3e (<= 0 required), their sum %.4f <= 0: %s; min distance at k=3: %.5f vs U = %.3f' % (maxpos3,
      sum(S[k] for k in range(3, K+1)), sum(S[k] for k in range(3, K+1)) <= 0, (x[3:] - x[:-3]).min(), Uf))
print('(2) weight change on k=1,2: |S_w - S_1| = %.3e <= bound 4 pi^2 M_f N0/L^2 with max single term %.3e' % (abs(S[1] + S[2] - S1[1] - S1[2]), maxw))
# (3) (H) at every consecutive pair, telescoping, beta term
Phi = np.array([fh(g[i]) + fh(g[i] + g[i+1]) + Af(g[i]) - Af(g[i+1]) + beta*(1 - g[i]) for i in range(N0 - 2)])
print('(3) max_i Phi(g_i, g_{i+1}) = %.6f <= c_cert = %.6f : %s ; argmax gaps (%.4f, %.4f)' % (Phi.max(), c, Phi.max() <= c,
      g[Phi.argmax()], g[Phi.argmax() + 1]))
lhs = sum(fh(g[i]) + fh(g[i] + g[i+1]) for i in range(N0 - 2))
rhs_sum = sum(c - Af(g[i]) + Af(g[i+1]) - beta*(1 - g[i]) for i in range(N0 - 2))
rhs_closed = (N0 - 2)*c + Af(g[N0-2]) - Af(g[0]) + beta*(x[N0-2] - x[0] - (N0 - 2))
print('    sum (H)-LHS = %.4f <= sum RHS = %.4f ; telescoped closed form %.4f (identity error %.2e)' % (lhs, rhs_sum, rhs_closed, rhs_sum - rhs_closed))
bound = N0*fh0 + 2*fh(g[N0-2]) + 2*rhs_closed
print('(4) chain: S_w = %.4f <= S_drop = %.4f + weight err <= bound = %.4f : %s' % (Sw, Sdrop, bound, Sw <= bound + 4*N0*maxw))
print('    beta-term 2 beta (x_{N0-1} - x_1 - N0) = %.2f ; #E_lambda = %d ; TL/2pi - N(T) = %.1f ; 3N/L = %.1f' % (2*beta*(x[N0-2] - x[0] - N0), inE.sum(), T*L/(2*math.pi) - NT, 3*NT/L))
print('    N0 (f^(0) + 2c) = %.2f ; per N(T): S_w/N = %.4f, hard-core ceiling f^(0)+2c = %.4f, J = 0.30684' % (N0*(fh0 + 2*c), Sw/NT, fh0 + 2*c))
