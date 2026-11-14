# E1Audit: end-to-end check of SL (12)+(14)+E1 Step 1 (x-space Bessel form, factor -pi, K-branch 2K_d) with complex shifts.
# primal  S = sum_n sigma_{a,b}(n) e(A n/q) F(n),  F(x) = A(xp/q), A = smooth Gaussian bump (complex), sigma_{a,b}(n) = sum_{uv=n} u^-a v^-b
# claim   S = q^{-1+a-b} zeta(1-a+b) Ft(1-a) + q^{-1-a+b} zeta(1+a-b) Ft(1-b)
#           + sum_n sigma_{-a,-b}(n) (2/p) (q/p)^{-z/2} n^{-z/2} [cos(pi d/2) e(Abar n/q) I_K(n) + e(-Abar n/q) I_osc(n)],
#         I_osc = -pi int A(xi) xi^{-z/2} (sin(pi d/2) J_d + cos(pi d/2) Y_d)(X) dxi,  I_K = 2 int A xi^{-z/2} K_d(X) dxi,  X = 4 pi sqrt(n xi/(pq)),
#         z = a+b, d = a-b; with a = alpha, b = -beta this is E1's z = alpha-beta, d = alpha+beta.
import sys, time
from mpmath import mp, mpf, mpc, exp, pi, sqrt, besselj, bessely, besselk, zeta, quad, sin, cos, linspace
mp.dps = 22
def run(alpha, beta, p, q, xic, w, nd, eps=-1):
    a, b = alpha, -beta; z, d = a + b, a - b
    An = eps * p; Ab = pow(An % q, -1, q) if q > 1 else 0
    A = lambda xi: exp(-((xi - xic) / w) ** 2 + 0.3j * (xi - xic) / w)
    lo, hi = max(mpf('0.5'), xic - 9 * w), xic + 9 * w
    cuts = linspace(lo, hi, 25)
    sig = lambda n, s1, s2: sum(mpf(u) ** (-s1) * mpf(n // u) ** (-s2) for u in range(1, n + 1) if n % u == 0)
    e = lambda x: exp(2j * pi * x)
    S = mpc(0)
    for n in range(1, int(hi * q / p) + 2):
        S += sig(n, a, b) * e(mpf(An * n) / q) * A(mpf(n) * p / q)
    At = lambda s: quad(lambda xi: A(xi) * xi ** (s - 1), cuts)
    Ft = lambda s: (mpf(q) / p) ** s * At(s)
    Res = mpf(q) ** (-1 + a - b) * zeta(1 - a + b) * Ft(1 - a) + mpf(q) ** (-1 - a + b) * zeta(1 + a - b) * Ft(1 - b)
    Dsum = mpc(0); last = []
    for n in range(1, nd + 1):
        X = lambda xi: 4 * pi * sqrt(n * xi / (p * q))
        Io = -pi * quad(lambda xi: A(xi) * xi ** (-z / 2) * (sin(pi * d / 2) * besselj(d, X(xi)) + cos(pi * d / 2) * bessely(d, X(xi))), cuts)
        IK = 2 * quad(lambda xi: A(xi) * xi ** (-z / 2) * besselk(d, X(xi)), cuts)
        c = sig(n, -a, -b) * (mpf(2) / p) * (mpf(q) / p) ** (-z / 2) * mpf(n) ** (-z / 2)
        term = c * (cos(pi * d / 2) * e(mpf(Ab * n) / q) * IK + e(-mpf(Ab * n) / q) * Io)
        Dsum += term; last.append(abs(term))
    tot = Res + Dsum
    print(f"alpha={alpha} beta={beta} (p,q)=({p},{q}) xi_c={xic} w={w} nd={nd}: |S|={mp.nstr(abs(S),8)} |Res|={mp.nstr(abs(Res),8)} "
          f"|dual|={mp.nstr(abs(Dsum),8)} relerr={mp.nstr(abs(S-tot)/abs(S),3)} last|term|={mp.nstr(last[-1],3)} K-branch|n=1|={mp.nstr(abs(IK),3) if nd==1 else ''}", flush=True)
t0 = time.time()
cases = [(mpc('0.03','0.02'), mpc('0.05','-0.01'), 1, 3, 60, 6, 24),
         (mpc('0.03','0.02'), mpc('0.05','-0.01'), 2, 7, 8, mpf('1.2'), int(sys.argv[1]) if len(sys.argv) > 1 else 40),
         (mpc('0.001','0.0005'), mpc('-0.0005','0.0003'), 1, 3, 60, 6, 24),
         (mpc('0.0','0.0'), mpc('0.0','0.0'), 1, 3, 60, 6, 24)]
for cs in cases[:int(sys.argv[2]) if len(sys.argv) > 2 else 4]:
    if cs[0] == 0 and cs[1] == 0: continue  # a = b needs the combined residue; skipped (covered by the near-coalescent case)
    run(*cs); print(f"  elapsed {time.time()-t0:.1f}s", flush=True)
