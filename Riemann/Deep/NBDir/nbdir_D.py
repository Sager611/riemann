"""nbdir: continuous Nyman-Beurling distance for a real primitive Dirichlet character.

D_chi(lam)^2 = dist^2_{L^2(0,inf)}(1_(0,1), span{A_chi(theta/x): lam<=theta<=1}), A_chi(y)=sum_{n<=y} chi(n).
Grid theta_j = e^{-j h}, j=0..L/h (L = log 1/lam). With f_j = A(theta_j/x)/sqrt(theta_j):
  Gram_ij = Ghat(|u_i-u_j|), Ghat(v) = e^{-v/2} K(e^v), K(r) = int_1^inf A(rz)A(z) z^-2 dz (r>=1),
  <1_(0,1), f_j> = e^{-u_j/2} L(1,chi), so D_h^2 = 1 - L(1)^2 e^T Gram^-1 e, e_j = e^{-u_j/2}.
K(r) = r sum_n chi(n) W(rn), W(a) = int_a^inf A y^-2 dy = A(a)/a + T(floor(a)+1), T(M) = sum_{m>=M} chi(m)/m,
     = c0 L(1) + r^-1 sum_n chi(n) V(rn)/n^2,  V(a) = a(A(a)-c0) + a^2 T(floor(a)+1)  (c0 = mean of A = L(0,chi)).
Usage: python3 nbdir_D.py <chi: m4|5|8> <N_terms> <h_min> <L_max> [check]
"""
import sys, time, json
import numpy as np
from scipy.linalg import cho_factor, cho_solve, eigvalsh, solve_toeplitz
import mpmath as mp

CHI = {"m4": (4, [0, 1, 0, -1], 1), "5": (5, [0, 1, -1, -1, 1], 0), "8": (8, [0, 1, 0, -1, 0, -1, 0, 1], 0)}
name = sys.argv[1]; NT = int(sys.argv[2]); hmin = float(sys.argv[3]); Lmax = float(sys.argv[4])
check = len(sys.argv) > 5 and sys.argv[5] == "check"
q, chil, par = CHI[name]
chiv = np.array(chil, float)
Atab = np.cumsum(chiv)            # A(m) for m = 0..q-1 (chi(0)=0) ; periodic mod q
c0 = Atab.mean()
mp.mp.dps = 30
L1 = float(mp.dirichlet(1, chil)); L0 = float(mp.dirichlet(0, chil))
assert abs(c0 - L0) < 1e-12, (c0, L0)

# T(M) exact table for M <= M1, asymptotic (digamma series, log1p) beyond
M1 = 2000
def T_asym(M):
    M = np.asarray(M, float); out = np.zeros_like(M)
    for j in range(q):
        x = M + j
        c = chiv[(M.astype(np.int64) + j) % q]
        out += c * (np.log1p(j / M) - q / (2 * x) - q**2 / (12 * x**2) + q**4 / (120 * x**4) - q**6 / (252 * x**6))
    return -out / q
Ttab = np.zeros(M1 + 2)
Ttab[M1 + 1] = T_asym(np.array([M1 + 1.0]))[0]
for M in range(M1, 0, -1):
    Ttab[M] = Ttab[M + 1] + chiv[M % q] / M
def Tfun(M):  # M int64 array >= 1
    out = np.empty(M.shape)
    sm = M <= M1 + 1
    out[sm] = Ttab[M[sm]]
    out[~sm] = T_asym(M[~sm])
    return out

n = np.arange(1, NT + 1, dtype=np.float64)
chin = chiv[np.arange(1, NT + 1) % q]
w2 = chin / n**2
def K_of_r(r):
    a = r * n
    fl = np.floor(a).astype(np.int64)
    V = a * (Atab[fl % q] - c0) + a * a * Tfun(fl + 1)
    return c0 * L1 + np.dot(w2, V) / r
# exact K(1) = (1/q) sum_b A(b)^2 [psi((b+1)/q) - psi(b/q)]
K1 = float(sum(mp.mpf(Atab[b])**2 * (mp.digamma(mp.mpf(b + 1) / q) - mp.digamma(mp.mpf(b) / q)) for b in range(1, q)) / q)

if check:  # brute force K(r) on [1,Z] by exact piecewise integration + c0^2/Z tail
    for r in [1.0, 1.37, 2.7182818, 7.389056]:
        Z = 200000
        bp = np.unique(np.concatenate([np.arange(1, Z + 1, dtype=float), np.arange(np.ceil(r), r * Z + 1) / r]))
        bp = bp[(bp >= 1) & (bp <= Z)]
        mid = 0.5 * (bp[:-1] + bp[1:])
        Az = Atab[np.floor(mid).astype(np.int64) % q]; Arz = Atab[np.floor(r * mid).astype(np.int64) % q]
        brute = np.sum(Az * Arz * (1 / bp[:-1] - 1 / bp[1:])) + c0**2 / Z
        print(f"check r={r}: formula {K_of_r(r):.15f} brute {brute:.15f} exactK1 {K1:.15f}")
    sys.exit(0)

t0 = time.time()
kmax = int(round(Lmax / hmin))
G = np.empty(kmax + 1)
G[0] = K1
for k in range(1, kmax + 1):
    v = k * hmin
    G[k] = np.exp(-v / 2) * K_of_r(np.exp(v))
print(f"chi={name} q={q} a={par} c0={c0} L(1)={L1:.15f} K(1)={K1:.15f} NT={NT} hmin={hmin} Lmax={Lmax} Ghat time {time.time()-t0:.1f}s", flush=True)
np.save(f"/Users/adriansager/playground/riemann/.scratch/nbdir/G_{name}_{NT}_{hmin}_{Lmax}.npy", G)

res = []
hm = 1
mults = []
while hmin * hm <= 0.1 + 1e-12:
    mults.append(hm); hm *= 2
for h_mult in mults[::-1]:
    h = hmin * h_mult
    for L in [2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 20, 24]:
        if L > Lmax + 1e-9: continue
        m = int(round(L / h)) + 1
        col = G[np.arange(m) * h_mult]
        e = np.exp(-np.arange(m) * h / 2)
        try:
            x = solve_toeplitz(col, e)
            D2 = 1 - L1**2 * e @ x
            if m <= 700:
                Tm = col[np.abs(np.arange(m)[:, None] - np.arange(m)[None, :])]
                xc = cho_solve(cho_factor(Tm, lower=True), e)
                assert abs((1 - L1**2 * e @ xc) - D2) < 1e-9, "toeplitz vs cholesky"
                ev = eigvalsh(Tm)
            else:
                ev = np.array([np.nan, np.nan])
            res.append(dict(h=h, L=L, m=m, D2=D2, LD2=L * D2, lmin=float(ev[0]), lmax=float(ev[-1]), cnorm=float(np.linalg.norm(L1 * x))))
            print(f"h={h:.4f} L={L:5.1f} m={m:5d} D2={D2:.10f} L*D2={L*D2:.8f} eig[{ev[0]:.2e},{ev[-1]:.2e}] |c|={np.linalg.norm(L1*x):.3e}", flush=True)
        except Exception as ex:
            print(f"h={h} L={L} failed {ex}", flush=True)
json.dump(res, open(f"/Users/adriansager/playground/riemann/.scratch/nbdir/D_{name}_{NT}_{hmin}_{Lmax}.json", "w"))
print(f"total {time.time()-t0:.1f}s")
