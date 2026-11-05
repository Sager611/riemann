# Test of the lattice identity (Theorem 3.1) on random test functions.
# Direct side: x-space Weil form (eq. 2.1) via galerkin_x.build (even) / build_odd (odd), no lattice.
# Lattice side: (1/2L) sum_k mu_L(r_k)|f^(r_k)|^2 + pole, symbol from cert.Sym (digamma + tail series + atoms).
import sys, random
from flint import arb, arb_mat, ctx
import galerkin_x as GX
import cert as CE

def build_odd(c, N, panels, nodes):
    L = arb(c).log(); pi = arb.pi()
    om = [2*pi*j/L for j in range(N+1)]
    Nj = (2/L).sqrt()
    gl = [arb.legendre_p_root(nodes, k, weight=True) for k in range(nodes)]
    S = [arb(0)]*(N+1); Cm = [arb(0)]*(N+1); Um = [arb(0)]*(N+1); KL = arb(0)
    for p in range(panels):
        a0 = L*p/panels; a1 = L*(p+1)/panels; h = (a1-a0)/2; m0 = (a1+a0)/2
        for (x, w) in gl:
            u = m0 + h*x; ww = w*h
            den = -((-2*u).expm1()); two = 2*(-u/2).exp()/den
            KL += ww*(2*(-2*u).exp()/den - two)
            for j in range(1, N+1):
                t = om[j]*u
                S[j] += ww*two*t.sin(); Cm[j] += ww*two*2*(t/2).sin()**2; Um[j] += ww*two*u*t.cos()
    C0 = -arb.const_euler() - pi.log()
    base = C0 - (-((-2*L).expm1())).log() + KL
    M = arb_mat(N, N)   # index j-1
    for j in range(1, N+1):
        for k in range(1, N+1):
            sg = 1 if (j+k) % 2 == 0 else -1
            if j == k: v = base + Cm[j] + Um[j]/L - S[j]/(L*om[j])
            else: v = -(Nj*Nj/2)*sg*((S[k]-S[j])/(om[j]-om[k]) + (S[k]+S[j])/(2*pi*(j+k)/L))
            M[j-1, k-1] = v
    n = 2
    while n < c:
        lam = GX.vm(n)
        if lam is not None:
            t = arb(n).log(); wn = lam/arb(n).sqrt()
            sn = [(om[j]*t).sin() for j in range(N+1)]; cs = [(om[j]*t).cos() for j in range(N+1)]
            for j in range(1, N+1):
                for k in range(1, N+1):
                    sg = 1 if (j+k) % 2 == 0 else -1
                    T1 = (L - t)*cs[k] if j == k else sg*(sn[k]-sn[j])/(om[j]-om[k])
                    T2 = -sg*(sn[k]+sn[j])/(2*pi*(j+k)/L)
                    M[j-1, k-1] = M[j-1, k-1] - 2*wn*(Nj*Nj/2)*(T1 - T2)
        n += 1
    sh = (L/4).sinh()
    b = [Nj*(-2)*om[j]*(1 if j % 2 == 0 else -1)*sh/(om[j]**2 + arb(1)/4) for j in range(1, N+1)]
    for j in range(N):
        for k in range(N): M[j, k] = M[j, k] - 2*b[j]*b[k]
    return M

def project(xi, vecs):
    # Gram-Schmidt: remove components along vecs (list of lists)
    import copy
    basis = []
    for v in vecs:
        w = list(v)
        for q in basis:
            d = sum((a*b for a, b in zip(w, q)), arb(0))
            w = [a - d*b for a, b in zip(w, q)]
        nrm = sum((a*a for a in w), arb(0)).sqrt()
        basis.append([a/nrm for a in w])
    x = list(xi)
    for q in basis:
        d = sum((a*b for a, b in zip(x, q)), arb(0))
        x = [a - d*b for a, b in zip(x, q)]
    return x

def run(c, N, K, sector, moments, seed, full=False):
    random.seed(seed)
    S = CE.Sym(c); L = S.L; pi = S.pi
    if sector == 'even':
        idx = list(range(0, N+1)); Nj = lambda j: (1/L).sqrt() if j == 0 else (2/L).sqrt()
        xi = [arb(random.uniform(-1, 1)) for _ in idx]
        vecs = [[Nj(j)*(1 if j % 2 == 0 else -1)*arb(4*j*j)**m for j in idx] for m in range(moments)]
        M = GX.build(c, N, 50, 70)
    else:
        idx = list(range(1, N+1)); Nj = lambda j: (2/L).sqrt()
        xi = [arb(random.uniform(-1, 1)) for _ in idx]
        vecs = [[Nj(j)*(1 if j % 2 == 0 else -1)*j*arb(4*j*j)**m for j in idx] for m in range(moments)]
        M = build_odd(c, N, 50, 70)
    if moments: xi = project(xi, vecs)
    n = len(idx)
    direct = arb(0)
    for a in range(n):
        for b in range(n): direct += xi[a]*xi[b]*M[a, b]
    norm2 = sum((x*x for x in xi), arb(0))
    # lattice
    def s_(m): return L/2 if m == 0 else L*(pi*m/2).sin()/(pi*m)
    def fhat(k):  # real number with |f^(r_k)| = |fhat|
        tot = arb(0)
        for a, j in enumerate(idx):
            if sector == 'even': tot += xi[a]*Nj(j)*(s_(2*j-k) + s_(2*j+k))
            else: tot += xi[a]*Nj(j)*(s_(2*j-k) - s_(2*j+k))
        return tot
    lat = arb(0); latfull = arb(0)
    for k in range(0, K+1):
        fh = fhat(k)
        if fh == 0 and k % 2 == 0: continue
        mu = S.mu(k)
        wgt = 1 if k == 0 else 2
        lat += wgt*mu*fh*fh
        latfull += wgt*(mu + S.T(k))*fh*fh
    sh = (L/4).sinh(); ch = (L/4).cosh()
    if sector == 'even':
        aa = sum((xi[a]*Nj(j)*(1 if j % 2 == 0 else -1)*sh/((2*pi*j/L)**2 + arb(1)/4) for a, j in enumerate(idx)), arb(0))
        pole = 2*aa*aa
    else:
        bb = sum((xi[a]*Nj(j)*(-2)*(2*pi*j/L)*(1 if j % 2 == 0 else -1)*sh/((2*pi*j/L)**2 + arb(1)/4) for a, j in enumerate(idx)), arb(0))
        pole = -2*bb*bb
    lat = lat/(2*L) + pole; latfull = latfull/(2*L) + pole
    # tail estimate for odd k > K (only matters without moment conditions)
    fb = fhat(K+1 if (K+1) % 2 == 1 else K+2)
    print(f"c={c} {sector} N={N} moments={moments} K={K} seed={seed}: direct={direct.str(25, radius=False)} lattice={lat.str(25, radius=False)} "
          f"rel.diff={((lat-direct)/direct).str(3, radius=False)} | full-symbol rel.diff={((latfull-direct)/direct).str(3, radius=False)} | |fhat(r_K+1)|={abs(fb).str(3, radius=False)} ||xi||^2={norm2.str(5, radius=False)}", flush=True)

if __name__ == "__main__":
    ctx.prec = 200
    c = int(sys.argv[1])
    for (sector, mom, K, seed) in [('even', 3, 1500, 1), ('even', 3, 1500, 2), ('odd', 3, 1500, 3), ('odd', 3, 1500, 4), ('even', 0, 6000, 5), ('odd', 0, 6000, 6)]:
        run(c, 16, K, sector, mom, seed)
