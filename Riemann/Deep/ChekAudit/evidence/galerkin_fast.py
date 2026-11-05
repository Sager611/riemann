# Fast closed-form version of galerkin_x (same x-space derivation): the 1D integrals
#   Cm_j = int_0^L 2varpi (1-cos w u), S_j = int_0^L 2varpi sin(w u), Um_j = int_0^L 2varpi u cos(w u)
# are written with digamma/trigamma at 1/4 + i w/2 and exponentially convergent tails over [L, inf).
# Validated against galerkin_x (pure quadrature) at small c.  lambda_min by inverse iteration.
# Usage: galerkin_fast.py c N bits [iters]
import sys, time, random
from flint import arb, acb, arb_mat, ctx
from galerkin_x import vm

def build(c, N):
    L = arb(c).log(); pi = arb.pi(); cc = arb(c)
    om = [2*pi*j/L for j in range(N+1)]
    Nj = [(1/L).sqrt()] + [(2/L).sqrt()]*N
    # tail coefficients
    tl = []
    n = 0
    while True:
        an = arb(2*n) + arb(1)/2; t = cc**(-an); tl.append((an, t))
        if (t < arb(2)**(-(ctx.prec+40))) and n > 2: break
        n += 1
    tau0 = sum((2*t/an for an, t in tl), arb(0))
    hp0 = acb(arb(1)/4, 0).digamma().real - pi.log()
    base = hp0 + tau0
    S = [arb(0)]*(N+1); Cm = [arb(0)]*(N+1); Um = [arb(0)]*(N+1)
    for j in range(N+1):
        w = om[j]; z = acb(arb(1)/4, w/2)
        ps = z.digamma(); p1 = z.polygamma(1)
        T = -2*sum((t*an/(an*an + w*w) for an, t in tl), arb(0))
        Cm[j] = (ps.real - pi.log()) - hp0 - tau0 - T
        S[j] = ps.imag - 2*w*sum((t/(an*an + w*w) for an, t in tl), arb(0))
        um = arb(0)
        for an, t in tl:
            s = acb(an, -w)
            um += t*(L/s + 1/(s*s)).real
        Um[j] = p1.real/2 - 2*um
    M = arb_mat(N+1, N+1)
    for j in range(N+1):
        for k in range(N+1):
            if j == k:
                v = base + Um[0]/L if j == 0 else base + Cm[j] + Um[j]/L + S[j]/(L*om[j])
            else:
                sg = 1 if (j+k) % 2 == 0 else -1
                v = -(Nj[j]*Nj[k]/2)*sg*(S[k]-S[j])/(om[j]-om[k]) + (Nj[j]*Nj[k]/2)*sg*(S[k]+S[j])/(2*pi*(j+k)/L)
            M[j, k] = v
    # atoms: -sum_n 2 w_n g_jk(ln n); g_jk linear in sin(w_j t), so sum over n first
    P = [arb(0)]*(N+1); Dg = [arb(0)]*(N+1)
    n = 2
    while n < c:
        lam = vm(n)
        if lam is not None:
            t = arb(n).log(); wn = lam/arb(n).sqrt()
            for j in range(N+1):
                sj, cj = (om[j]*t).sin(), (om[j]*t).cos()
                P[j] += wn*sj
                Dg[j] += 2*wn*((1 - t/L) if j == 0 else ((L - t)*cj - sj/om[j])/L)
        n += 1
    for j in range(N+1):
        for k in range(N+1):
            if j == k: M[j, k] = M[j, k] - Dg[j]
            else:
                sg = 1 if (j+k) % 2 == 0 else -1
                M[j, k] = M[j, k] - Nj[j]*Nj[k]*sg*((P[k]-P[j])/(om[j]-om[k]) - (P[k]+P[j])/(2*pi*(j+k)/L))
    sh = (L/4).sinh()
    a = [Nj[j]*(1 if j % 2 == 0 else -1)*sh/(om[j]**2 + arb(1)/4) for j in range(N+1)]
    for j in range(N+1):
        for k in range(N+1): M[j, k] = M[j, k] + 2*a[j]*a[k]
    return M

def lam_min_inv(M, iters=4, seed=0):
    n = M.nrows()
    R = arb_mat(n, n)
    for i in range(n):
        for j in range(n): R[i, j] = M[i, j].mid()
    random.seed(seed)
    x = arb_mat(n, 1, [arb(random.uniform(-1, 1)) for _ in range(n)])
    lam = None; hist = []
    for it in range(iters):
        y = R.solve(x, algorithm="approx")
        nrm = sum((y[i, 0]**2 for i in range(n)), arb(0)).sqrt()
        x = arb_mat(n, 1, [(y[i, 0]/nrm).mid() for i in range(n)])
        Mx = R*x
        lam = sum((x[i, 0]*Mx[i, 0] for i in range(n)), arb(0))
        res = sum(((Mx[i, 0]-lam*x[i, 0])**2 for i in range(n)), arb(0)).sqrt()
        hist.append((lam.mid(), res.mid()))
    return lam, hist, x

if __name__ == "__main__":
    c = int(sys.argv[1]); N = int(sys.argv[2]); ctx.prec = int(sys.argv[3])
    iters = int(sys.argv[4]) if len(sys.argv) > 4 else 4
    t0 = time.time()
    M = build(c, N); t1 = time.time()
    lam, hist, x = lam_min_inv(M, iters)
    print(f"c={c} N={N} prec={ctx.prec} lam_min(Q_c,N)={lam.mid().str(10, radius=False)} build={t1-t0:.1f}s eig={time.time()-t1:.1f}s", flush=True)
    for (l, r) in hist: print(f"    iter lam={l.str(10, radius=False)} resid={r.str(3, radius=False)}")
