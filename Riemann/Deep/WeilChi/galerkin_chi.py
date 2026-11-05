# galerkin_chi.py: Rayleigh-Ritz UPPER bounds for lambda_min(Q_chi) on even f in L^2(-L/2, L/2), L = ln c, in the Neumann cosine basis
# phi_j (j = 0..N), from the x-space form of Q_chi (closed-form integrals). Port of ChekAudit/evidence/galerkin_fast.py:
# z = 1/4 + a/2 + i w/2 instead of 1/4 + i w/2, tail exponents b_n = 2n + 1/2 + a, constant log(q/pi), prime weights chi(n) Lambda(n)/sqrt n,
# no pole term.  Independent of the lattice identity (x-space quadrature in closed form).  Usage: galerkin_chi.py CHAR c N bits [iters]
import sys, time, random
from flint import arb, acb, arb_mat, ctx
CH = {'m4': ([0,1,0,-1], 4, 1), 'm3': ([0,1,-1], 3, 1), 'p5': ([0,1,-1,-1,1], 5, 0), 'p8': ([0,1,0,-1,0,-1,0,1], 8, 0)}
def vm(n):
    for p in range(2, n+1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return arb(p).log() if m == 1 else None
    return None
def build(key, c, N):
    tab, q, a = CH[key]
    L = arb(c).log(); pi = arb.pi(); cc = arb(c)
    om = [2*pi*j/L for j in range(N+1)]
    Nj = [(1/L).sqrt()] + [(2/L).sqrt()]*N
    tl = []; n = 0
    while True:
        bn = arb(2*n) + arb(1)/2 + a; t = cc**(-bn); tl.append((bn, t))
        if (t < arb(2)**(-(ctx.prec+40))) and n > 2: break
        n += 1
    z0 = arb(1)/4 + arb(a)/2
    tau0 = sum((2*t/bn for bn, t in tl), arb(0))
    hp0 = acb(z0, 0).digamma().real + (arb(q)/pi).log()
    base = hp0 + tau0
    S = [arb(0)]*(N+1); Cm = [arb(0)]*(N+1); Um = [arb(0)]*(N+1)
    for j in range(N+1):
        w = om[j]; z = acb(z0, w/2)
        ps = z.digamma(); p1 = z.polygamma(1)
        T = -2*sum((t*bn/(bn*bn + w*w) for bn, t in tl), arb(0))
        Cm[j] = ps.real - acb(z0, 0).digamma().real - tau0 - T
        S[j] = ps.imag - 2*w*sum((t/(bn*bn + w*w) for bn, t in tl), arb(0))
        um = arb(0)
        for bn, t in tl:
            s = acb(bn, -w); um += t*(L/s + 1/(s*s)).real
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
    P = [arb(0)]*(N+1); Dg = [arb(0)]*(N+1)
    n = 2
    while n < c:
        lam = vm(n); x = tab[n % q]
        if lam is not None and x != 0:
            t = arb(n).log(); wn = x*lam/arb(n).sqrt()
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
    return M
def lam_min_inv(M, iters=5, seed=0):
    n = M.nrows(); R = arb_mat(n, n)
    for i in range(n):
        for j in range(n): R[i, j] = M[i, j].mid()
    random.seed(seed)
    x = arb_mat(n, 1, [arb(random.uniform(-1, 1)) for _ in range(n)]); lam = None
    for it in range(iters):
        y = R.solve(x, algorithm="approx")
        nrm = sum((y[i, 0]**2 for i in range(n)), arb(0)).sqrt()
        x = arb_mat(n, 1, [(y[i, 0]/nrm).mid() for i in range(n)])
        Mx = R*x; lam = sum((x[i, 0]*Mx[i, 0] for i in range(n)), arb(0))
    return lam
if __name__ == "__main__":
    key = sys.argv[1]; c = int(sys.argv[2]); ctx.prec = int(sys.argv[4])
    for N in [int(v) for v in sys.argv[3].split(',')]:
        t0 = time.time(); M = build(key, c, N); lam = lam_min_inv(M)
        print(f"{key} c={c} N={N} prec={ctx.prec} lambda_min(Q_N, even) = {lam.mid().str(10, radius=False)}  ({time.time()-t0:.1f}s)", flush=True)
