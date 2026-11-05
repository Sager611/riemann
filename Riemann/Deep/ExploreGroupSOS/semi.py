# semi.py: even-sector Galerkin lambda_min of semilocal Weil forms Q_S on windows [-L/2, L/2] (L real),
# S = set of finite primes kept (all their powers n < e^L). Copy of ChekAudit galerkin_fast.build with L real + atom filter.
import sys, math, time
from flint import arb, acb, arb_mat, ctx
sys.path.insert(0, '/Users/adriansager/playground/riemann/Riemann/Deep/ChekAudit/evidence')
import galerkin_fast as gf

def ppow(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return p if m == 1 else None
    return None

def build_L(L, S, N):
    pi = arb.pi(); cL = L.exp()
    om = [2*pi*j/L for j in range(N+1)]
    Nj = [(1/L).sqrt()] + [(2/L).sqrt()]*N
    tl = []; n = 0
    while True:
        an = arb(2*n) + arb(1)/2; t = (-an*L).exp(); tl.append((an, t))
        if (t < arb(2)**(-(ctx.prec+40))) and n > 2: break
        n += 1
    tau0 = sum((2*t/an for an, t in tl), arb(0))
    hp0 = acb(arb(1)/4, 0).digamma().real - pi.log()
    base = hp0 + tau0
    Sv = [arb(0)]*(N+1); Cm = [arb(0)]*(N+1); Um = [arb(0)]*(N+1)
    for j in range(N+1):
        w = om[j]; z = acb(arb(1)/4, w/2)
        ps = z.digamma(); p1 = z.polygamma(1)
        T = -2*sum((t*an/(an*an + w*w) for an, t in tl), arb(0))
        Cm[j] = (ps.real - pi.log()) - hp0 - tau0 - T
        Sv[j] = ps.imag - 2*w*sum((t/(an*an + w*w) for an, t in tl), arb(0))
        um = arb(0)
        for an, t in tl:
            s = acb(an, -w); um += t*(L/s + 1/(s*s)).real
        Um[j] = p1.real/2 - 2*um
    M = arb_mat(N+1, N+1)
    for j in range(N+1):
        for k in range(N+1):
            if j == k:
                v = base + Um[0]/L if j == 0 else base + Cm[j] + Um[j]/L + Sv[j]/(L*om[j])
            else:
                sg = 1 if (j+k) % 2 == 0 else -1
                v = -(Nj[j]*Nj[k]/2)*sg*(Sv[k]-Sv[j])/(om[j]-om[k]) + (Nj[j]*Nj[k]/2)*sg*(Sv[k]+Sv[j])/(2*pi*(j+k)/L)
            M[j, k] = v
    P = [arb(0)]*(N+1); Dg = [arb(0)]*(N+1)
    n = 2; cf = float(cL.mid())
    while n < cf:
        p = ppow(n)
        if p is not None and (S is None or p in S):
            t = arb(n).log(); wn = arb(p).log()/arb(n).sqrt()
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

def lmin(L, S, N=60):
    M = build_L(L, S, N)
    # symmetric eigenvalues via approx: use inverse iteration from gf (works if lam_min is isolated); also detect negativity
    R = arb_mat(N+1, N+1)
    for i in range(N+1):
        for j in range(N+1): R[i, j] = M[i, j].mid()
    import numpy as np
    A = np.array([[float(R[i, j]) for j in range(N+1)] for i in range(N+1)])
    ev = np.linalg.eigvalsh(A)
    return ev[0], ev[1]

if __name__ == '__main__':
    ctx.prec = 200
    mode = sys.argv[1]
    Ss = {'inf': set(), 'inf2': {2}, 'inf23': {2, 3}, 'inf235': {2, 3, 5}, 'all': None}
    if mode == 'scan':
        for name in ['inf', 'inf2', 'inf23', 'inf235', 'all']:
            row = []
            for a in [0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.70, 0.80, 0.90, 1.0, 1.2]:
                L = arb(2*a)  # window half-width a -> L = 2a
                e0, e1 = lmin(L, Ss[name], 40)
                row.append(f'{a}:{e0:.2e}')
            print(name, ' '.join(row), flush=True)
