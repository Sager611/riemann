# pyright: reportMissingImports=false
# python-flint lives in a local venv (see README.md), not the system interpreter.

"""Rigorous lower bound for the Weil functional on supp f in [-a,a], a = L/2 (window L).

Chain (all steps rigorous; arb interval arithmetic throughout):
 (1) Q(f) = (1/2pi) int |fhat|^2 S + 2 Re(A conj B),  S(r) = Phi(r) - p(r),
     Phi(r) = Re psi(1/4+ir/2) - log pi,  p(r) = 2 sum_{n<e^L} Lambda(n) n^{-1/2} cos(r log n),
     A = int f e^{x/2}, B = int f e^{-x/2}           [explicit formula, h = |fhat|^2]
 (2) Phi increasing on r>=0 (Re psi(x+iy) = psi(x) + sum_k y^2/((k+x)((k+x)^2+y^2))), p <= P
     => S(r) >= eps := Phi(R) - P for |r| >= R.
 (3) Q(f) >= <M f,f>,  M = eps I + C,  C = K + Pi,
     K = band-limited operator, symbol (S - eps) 1_{[-R,R]};  Pi = A(x)B + B(x)A (poles).
 (4) Legendre ON basis phi_j, j<N:  lam_min(M) >= min(eps, lam_min(eps I + C_NN)) - ||E||,
     ||E|| <= 2 ||P_perp C||_HS <= 2 sqrt(sum_{j>=N} (sup|w| eps_j + 2|A_j| ||A||)^2),
     eps_j^2 = (1/2pi) int_{-R}^{R} |phihat_j|^2 <= (2/pi)(aR)^{2j+1}/((2j+1)!!)^2.
 (5) C_NN by composite Gauss-Legendre in r (unit panels, n nodes), error matrix Delta bounded by
     the Bernstein-ellipse estimate (64/15) M rho^{-2n}/(rho^2-1) * half-length per panel.
 (6) eps I + C_NN - (||E|| + ||Delta||_F + delta) I positive definite (interval block-Schur LDL)
     => Q(f) >= delta ||f||^2 for every complex f supported in [-a,a].
"""
import sys, time, math
from flint import arb, arb_mat, acb, ctx
import weil

DIAG = True
def phi(r):
    return acb(arb(1)/4, r/2).digamma().real - arb.pi().log()

def sbessel(nmax, a, r, bprec):
    old = ctx.prec; ctx.prec = bprec
    x = a*r                     # exact: a, r are exact points with short mantissas
    assert x.rad() == 0
    c = (arb.pi()/(2*x)).sqrt()
    out = [None]*(nmax+2)
    out[nmax+1] = c*x.bessel_j(arb(nmax+1)+arb(1)/2)
    out[nmax] = c*x.bessel_j(arb(nmax)+arb(1)/2)
    for n in range(nmax, 0, -1):
        out[n-1] = (2*n+1)/x*out[n] - out[n+1]
    ctx.prec = old
    # flush huge-exponent-range tiny values to a zero-centred ball (rigorous; avoids inf in
    # arb_mat block multiplication when exponents span thousands of decades)
    tiny = arb(2)**-400
    return [v + 0 if abs(v) > tiny else arb(0, abs(v).upper()) for v in out[:nmax]]

def dfact(m):  # (2m+1)!!
    r = arb(1)
    for k in range(1, 2*m+2, 2): r *= k
    return r

def chol_mid(G, leaf=32):
    """Approximate Cholesky (midpoints only, radii stripped). Returns lower-triangular point
    matrix L or None if a pivot is not positive. Rigor comes from the residual check, not here."""
    n = G.nrows()
    if n <= leaf:
        A = [[G[i, j].mid() for j in range(n)] for i in range(n)]
        Lm = arb_mat(n, n)
        Lr = [[arb(0)]*n for _ in range(n)]
        for j in range(n):
            s = A[j][j] - sum((Lr[j][k]**2 for k in range(j)), arb(0))
            s = s.mid()
            if not (s > 0): return None
            Lr[j][j] = s.sqrt().mid()
            for i in range(j+1, n):
                t = A[i][j] - sum((Lr[i][k]*Lr[j][k] for k in range(j)), arb(0))
                Lr[i][j] = (t/Lr[j][j]).mid()
        for i in range(n):
            for j in range(i+1): Lm[i, j] = Lr[i][j]
        return Lm
    k = n//2
    def blk(M, r0, r1, c0, c1):
        X = arb_mat(r1-r0, c1-c0)
        for i in range(r0, r1):
            for j in range(c0, c1): X[i-r0, j-c0] = M[i, j]
        return X
    L11 = chol_mid(blk(G, 0, k, 0, k), leaf)
    if L11 is None: return None
    X = L11.solve(blk(G, 0, k, k, n)).mid()        # L11 X = B,  L21 = X^T
    S = (blk(G, k, n, k, n) - X.transpose()*X).mid()
    L22 = chol_mid(S, leaf)
    if L22 is None: return None
    Lm = arb_mat(n, n)
    for i in range(k):
        for j in range(i+1): Lm[i, j] = L11[i, j]
    for i in range(n-k):
        for j in range(k): Lm[k+i, j] = X[j, i]
        for j in range(i+1): Lm[k+i, k+j] = L22[i, j]
    return Lm

def certified_shift(H):
    """Returns a rigorous upper bound b with lambda_min(H) >= -b, or None if Cholesky fails.
    H = Hmid + Rad,  Hmid = L L^T + Res  (L point matrix => L L^T >= 0)."""
    n = H.nrows()
    Hm = H.mid()
    rad2 = arb(0)
    for i in range(n):
        for j in range(n): rad2 += arb(H[i, j].rad())**2
    Lm = chol_mid(Hm)
    if Lm is None: return None
    Res = Hm - Lm*Lm.transpose()
    res2 = arb(0)
    for i in range(n):
        for j in range(n): res2 += arb(abs(Res[i, j]).upper())**2
    return arb((rad2.sqrt() + res2.sqrt()).upper())

def is_pd(G, leaf=40):
    n = G.nrows()
    if n <= leaf:
        A = [[G[i, j] for j in range(n)] for i in range(n)]
        for k in range(n):
            if not (A[k][k] > 0): return False
            for i in range(k+1, n):
                f = A[i][k]/A[k][k]
                for j in range(k+1, i+1):
                    A[i][j] -= f*A[j][k]
                    A[j][i] = A[i][j]
        return True
    k = n//2
    def blk(r0, r1, c0, c1):
        M = arb_mat(r1-r0, c1-c0)
        for i in range(r0, r1):
            for j in range(c0, c1): M[i-r0, j-c0] = G[i, j]
        return M
    A = blk(0, k, 0, k); B = blk(0, k, k, n); D = blk(k, n, k, n)
    if not is_pd(A, leaf): return False
    S = D - B.transpose()*A.solve(B)
    for i in range(S.nrows()):          # symmetrize (intersection-free: take union by averaging)
        for j in range(i):
            v = (S[i, j] + S[j, i])/2; S[i, j] = v; S[j, i] = v
    return is_pd(S, leaf)

def certify(L, R, N, n_per, deltas, b_e=arb(2)/5, chunk=400, log=print):
    a = L/2; x_max = a*R
    import envelope
    pr = envelope.primes_upto(L)
    P = sum((2*c for _, c in pr), arb(0))
    eps = phi(R) - P
    assert eps > 0, "R too small"
    nmax_n = max(n for n, _ in pr)
    bprec = ctx.prec + int(float(x_max)*math.log2(3)) + 128
    log(f"L={L.str(8, radius=False)} a={float(a):.4f} R={float(R)} primes={[n for n,_ in pr]} P={P.str(8,radius=False)} eps={eps.str(8,radius=False)} N={N} n/panel={n_per} prec={ctx.prec} bprec={bprec}")
    # --- (4) tail bound
    phi0 = phi(arb(0))
    supw = max(abs(phi0).upper(), phi(R).upper()) + P.upper() + eps.upper()
    supw = arb(supw)
    normA = (2*a.sinh()).sqrt()
    terms = []
    df = dfact(N)
    for j in range(N, N + 400):
        if j > N: df *= (2*j+1)
        ej = ((2/arb.pi())*x_max**(2*j+1)).sqrt()/df
        Aj = (2*(2*j+1)*a).sqrt()*(a/2)**j*(a/2).exp()/df
        terms.append((supw*ej + 2*Aj*normA)**2)
    rat = terms[-1]/terms[-2]
    assert rat < 1
    tail2 = sum(terms, arb(0)) + terms[-1]*rat/(1 - rat)
    Ebound = 2*tail2.sqrt()
    log(f"  tail: ||E|| <= {Ebound.str(5, radius=False)}  (sup|w| <= {supw.str(5,radius=False)})")
    # --- (5) quadrature error bound
    rho = (b_e + (b_e**2 + arb(1)/4).sqrt())*2
    smax = arb(1)/4 + (R + 2)/2
    Mpsi = (smax + 1).log() + arb.pi()/2 + 1 + 1/(arb(1)/4 - b_e/2)
    Mw = Mpsi + arb.pi().log() + eps + (b_e*arb(nmax_n).log()).cosh()*P
    Mell = Mw*(2*a*b_e).exp()
    npan = int(math.ceil(float(R)))
    errI = npan*(arb(64)/15)*Mell*(arb(1)/2)*rho**(-2*n_per)/(rho**2 - 1)
    Dbound = (2*a/arb.pi())*arb(N)**2*errI
    log(f"  quadrature: rho={rho.str(5,radius=False)} ||Delta||_F <= {Dbound.str(5, radius=False)}")
    # --- (5) assemble C_NN
    gl = weil.gauss(n_per)
    nodes = []; maxrad = 0.0
    for k in range(npan):
        lo = arb(k); hi = arb(k+1) if k+1 < float(R) else R
        hh = (hi - lo)/2; md = (hi + lo)/2
        for (x, w) in gl:
            rr = md + hh*x
            nodes.append((arb(rr.mid()), w*hh))     # exact point node; displacement bounded below
            maxrad = max(maxrad, float(rr.rad()))
    # perturbed-node error: |f'| <= Mell/b_e on the real line (Cauchy), total weight = R
    errP = R*Mell/b_e*arb(maxrad)*arb(N)
    Dbound += (2*a/arb.pi())*arb(N)**2*errP
    log(f"  node perturbation: max node radius {maxrad:.2e}; ||Delta||_F <= {Dbound.str(5, radius=False)}")
    ne, no = (N+1)//2, N//2
    K = [arb_mat(ne, ne), arb_mat(no, no)]
    t0 = time.time()
    for c0 in range(0, len(nodes), chunk):
        ch = nodes[c0:c0+chunk]; m = len(ch)
        Be, Bo, We, Wo = arb_mat(ne, m), arb_mat(no, m), arb_mat(ne, m), arb_mat(no, m)
        for q, (r, w) in enumerate(ch):
            wt = phi(r) - eps
            for n, c in pr: wt -= 2*c*(r*arb(n).log()).cos()
            wt *= w
            js = sbessel(N, a, r, bprec)
            for j in range(N):
                s = js[j] if (j//2) % 2 == 0 else -js[j]
                if j % 2 == 0: Be[j//2, q] = s; We[j//2, q] = s*wt
                else: Bo[j//2, q] = s; Wo[j//2, q] = s*wt
        K[0] += We*Be.transpose(); K[1] += Wo*Bo.transpose()
        if (c0//chunk) % 10 == 0: log(f"    nodes {c0+m}/{len(nodes)} {time.time()-t0:.0f}s")
    Av = [(2*(2*j+1)*a).sqrt()*(arb.pi()/a).sqrt()*(a/2).bessel_i(arb(j)+arb(1)/2) for j in range(N)]
    global LAST_TAU
    LAST_TAU = Ebound + Dbound
    res = {}
    for p in (0, 1):
        idx = list(range(p, N, 2)); m = len(idx)
        G = arb_mat(m, m)
        for i, j in enumerate(idx):
            for k2, l in enumerate(idx):
                G[i, k2] = (2*a/arb.pi())*(arb(2*j+1)*(2*l+1)).sqrt()*K[p][i, k2] + (2 if p == 0 else -2)*Av[j]*Av[l]
            G[i, i] += eps
        rad = max(float(G[i, k2].rad()) for i in range(m) for k2 in range(m))
        if DIAG: log(f"  diag {'even' if p==0 else 'odd'}: max radius {rad:.3e}, approx lam_min(G) = {weil.lam_min(G.mid()).str(6, radius=False)}")
        ok = None
        for d in deltas:
            tau = Ebound + Dbound + arb(d)
            H = arb_mat(m, m)
            for i in range(m):
                for k2 in range(m): H[i, k2] = G[i, k2]
                H[i, i] -= tau
            t1 = time.time()
            b = certified_shift(H)
            if b is not None and b < arb(d)/1000:
                ok = d*(1 - 1e-3); log(f"  {'even' if p==0 else 'odd'}: lam_min >= delta={d}(1-1e-3)  [residual+radius bound {b.str(3,radius=False)}, max entry radius {rad:.2e}, {time.time()-t1:.0f}s]"); break
            else:
                log(f"  {'even' if p==0 else 'odd'}: PD NOT certified at delta={d} ({time.time()-t1:.0f}s)")
        res[p] = ok
    r0, r1 = res[0], res[1]
    if r0 is None or r1 is None or not (eps > min(r0, r1) + Ebound):
        log(f"RESULT L={L.str(8,radius=False)}: not certified")
        return False
    log(f"RESULT L={L.str(8,radius=False)}: CERTIFIED Q(f) >= {min(r0, r1)} ||f||^2")
    return True

if __name__ == "__main__":
    ctx.prec = int(sys.argv[1])
    L = weil.parseL(sys.argv[2]); R = arb(sys.argv[3]); N = int(sys.argv[4]); n_per = int(sys.argv[5])
    deltas = [float(s) for s in sys.argv[6].split(',')]
    certify(L, R, N, n_per, deltas, log=lambda s: print(s, flush=True))
