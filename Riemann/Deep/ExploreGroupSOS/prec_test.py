# prec_test.py: Ozawa-style "numerical SOS + exact residual absorption" on the even Galerkin matrix Q_N of the Weil window form.
# For working precision b: round Q to b bits, Cholesky of (Q_b - tau I) in b-bit arithmetic, exact residual E = Q - tau I - L L^T,
# success iff ||E||_F < tau (then Q >= tau - ||E|| > 0 on the Galerkin space). tau = lambda_min/2.
import sys, mpmath as mp, time
from flint import arb, ctx
sys.path.insert(0, '/Users/adriansager/playground/riemann/Riemann/Deep/ChekAudit/evidence')
import galerkin_fast as gf
def run(c, N, hp_bits):
    ctx.prec = hp_bits; mp.mp.prec = hp_bits
    M = gf.build(c, N)
    Q = mp.matrix(N+1, N+1)
    for i in range(N+1):
        for j in range(N+1): Q[i, j] = mp.mpf(M[i, j].mid().str(int(hp_bits*0.31)+10, radius=False))
    lam, hist, x = gf.lam_min_inv(M, 4)
    lam = mp.mpf(lam.mid().str(30, radius=False))
    nrm = max(abs(v) for v in mp.eigsy(Q, eigvals_only=True)) if N <= 40 else None
    tau = lam/2
    res = []
    for b in [53, 64, 80, 96, 112, 128, 160, 192, 224, 256, 288, 320]:
        if b >= hp_bits - 20: break
        mp.mp.prec = b
        A = mp.matrix(N+1, N+1)
        for i in range(N+1):
            for j in range(N+1): A[i, j] = +Q[i, j] - (tau if i == j else 0)
        ok = True
        try:
            Lc = mp.cholesky(A)
        except Exception as e:
            ok = False
        mp.mp.prec = hp_bits
        if ok:
            E2 = mp.mpf(0)
            for i in range(N+1):
                for j in range(N+1):
                    s = sum(Lc[i, k]*Lc[j, k] for k in range(min(i, j)+1))
                    e = Q[i, j] - (tau if i == j else 0) - s
                    E2 += e*e
            ef = mp.sqrt(E2)
            ok2 = ef < tau
            res.append((b, 'chol ok', mp.nstr(ef, 3), 'CERT' if ok2 else 'fail'))
            if ok2: break
        else:
            res.append((b, 'chol FAILED (not PD at b bits)'))
    return lam, res
for c, N, hp in [(5, 60, 400), (10, 60, 500), (16, 80, 700)]:
    t0 = time.time()
    lam, res = run(c, N, hp)
    print(f'c={c} N={N} lam_min(Q_N)={mp.nstr(lam, 4)} ({time.time()-t0:.0f}s):', '; '.join(str(r) for r in res), flush=True)
