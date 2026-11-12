# Independent x-space Galerkin matrix of the Weil form Q_c in the cosine basis phi_j = N_j cos(2 pi j x/L) on [-L/2, L/2].
# Archimedean term by Weil's regularisation (paper eq. (2.1), = PROOF.md a(t) form), integrals by composite Gauss in arb.
# No lattice identity, no digamma, no Cauchy formulas.  Usage: galerkin_x.py c N bits [panels nodes]
import sys, time
from flint import arb, arb_mat, acb_mat, ctx

def vm(n):
    for p in range(2, n+1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return arb(p).log() if m == 1 else None
    return None

def build(c, N, panels, nodes):
    L = arb(c).log(); pi = arb.pi()
    om = [2*pi*j/L for j in range(N+1)]
    Nj = [(1/L).sqrt()] + [(2/L).sqrt()]*N
    gl = [arb.legendre_p_root(nodes, k, weight=True) for k in range(nodes)]
    S = [arb(0)]*(N+1); Cm = [arb(0)]*(N+1); Um = [arb(0)]*(N+1); KL = arb(0)
    for p in range(panels):
        a0 = L*p/panels; a1 = L*(p+1)/panels; h = (a1-a0)/2; m0 = (a1+a0)/2
        for (x, w) in gl:
            u = m0 + h*x; ww = w*h
            den = -((-2*u).expm1())          # 1 - e^{-2u}
            two = 2*(-u/2).exp()/den          # 2 varpi(u)
            KL += ww*(2*(-2*u).exp()/den - two)
            for j in range(N+1):
                t = om[j]*u
                s, co = t.sin(), t.cos()
                S[j] += ww*two*s
                Cm[j] += ww*two*2*(t/2).sin()**2
                Um[j] += ww*two*u*co
    C0 = -arb.const_euler() - pi.log()
    base = C0 - (-((-2*L).expm1())).log() + KL
    M = arb_mat(N+1, N+1)
    for j in range(N+1):
        for k in range(N+1):
            if j == k:
                if j == 0:
                    v = base + Um[0]*0  # placeholder
                    # int 2varpi u/L = Um[0]/L  (cos(0)=1)
                    v = base + Um[0]/L
                else:
                    v = base + Cm[j] + Um[j]/L + S[j]/(L*om[j])
            else:
                sg = 1 if (j+k) % 2 == 0 else -1
                v = -(Nj[j]*Nj[k]/2)*sg*((S[k]-S[j])/om[j-k] if j > k else (S[k]-S[j])/(-om[k-j])) 
                v += (Nj[j]*Nj[k]/2)*sg*(S[k]+S[j])/(2*pi*(j+k)/L)
            M[j, k] = v
    # atoms: -sum w_n * 2 g_jk(ln n)
    n = 2
    while n < c:
        lam = vm(n)
        if lam is not None:
            t = arb(n).log(); wn = lam/arb(n).sqrt()
            sn = [(om[j]*t).sin() for j in range(N+1)]
            cs = [(om[j]*t).cos() for j in range(N+1)]
            for j in range(N+1):
                for k in range(N+1):
                    sg = 1 if (j+k) % 2 == 0 else -1
                    if j == k: T1 = (L - t)*cs[k]
                    else: T1 = sg*(sn[k]-sn[j])/(om[j]-om[k])
                    if j + k == 0: T2 = L - t
                    else: T2 = -sg*(sn[k]+sn[j])/(2*pi*(j+k)/L)
                    g = Nj[j]*Nj[k]/2*(T1+T2)
                    M[j, k] = M[j, k] - 2*wn*g
        n += 1
    # pole: 2 a_j a_k, a_j = <phi_j, cosh(x/2)>
    sh = (L/4).sinh()
    a = [Nj[j]*(1 if j % 2 == 0 else -1)*sh/(om[j]**2 + arb(1)/4) for j in range(N+1)]
    for j in range(N+1):
        for k in range(N+1):
            M[j, k] = M[j, k] + 2*a[j]*a[k]
    return M

def lam_min(M):
    n = M.nrows(); R = arb_mat(n, n)
    for i in range(n):
        for j in range(n): R[i, j] = M[i, j].mid()
    ev = acb_mat(R).eig(algorithm="approx")
    return sorted([e.real.mid() for e in ev], key=lambda z: float(z))

if __name__ == "__main__":
    c = int(sys.argv[1]); N = int(sys.argv[2]); ctx.prec = int(sys.argv[3])
    panels = int(sys.argv[4]) if len(sys.argv) > 4 else 40
    nodes = int(sys.argv[5]) if len(sys.argv) > 5 else 60
    t0 = time.time()
    M = build(c, N, panels, nodes)
    ev = lam_min(M)
    print(f"c={c} N={N} (N+1 funcs) prec={ctx.prec} panels={panels}x{nodes} lam_min={ev[0].str(12, radius=False)} next={ev[1].str(5, radius=False)} t={time.time()-t0:.1f}s", flush=True)
    Mn = arb_mat(N, N)
    for i in range(N):
        for j in range(N): Mn[i, j] = M[i, j]
    ev2 = lam_min(Mn)
    print(f"   (first N funcs) lam_min={ev2[0].str(12, radius=False)}", flush=True)
