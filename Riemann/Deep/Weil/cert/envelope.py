# pyright: reportMissingImports=false
# python-flint lives in a local venv (see README.md), not the system interpreter.

"""Envelope reduction (pointwise lower bound on the Weil symbol above frequency R):
 Q(f) >= eps*||f||^2 + (1/2pi) int_{-R}^{R} |fhat|^2 (S(r)-eps) dr + 2AB  =: <M_R f, f>,
 S(r) = Re psi(1/4+ir/2) - log pi - 2 sum_{n<e^L} Lambda(n) n^{-1/2} cos(r log n),
 eps = Phi(R) - P (Phi increasing; P = 2 sum Lambda(n) n^{-1/2}).
 M_R = eps*I + compact; Legendre Galerkin. This file: floating (non-rigorous) feasibility test."""
import sys, time, math
from flint import arb, arb_mat, acb_mat, acb, ctx
import weil

def phi(r):  # Re psi(1/4 + i r/2) - log pi
    return acb(arb(1)/4, r/2).digamma().real - arb.pi().log()

def primes_upto(L):
    out = []; n = 2
    while math.log(n) < float(L) - 1e-15:
        lam = weil.vm(n)
        if lam is not None: out.append((n, lam/arb(n).sqrt()))
        n += 1
    return out

BPREC = 1600
def sbessel_all(nmax, x):
    """j_0..j_nmax at x>0 by downward recurrence from arb Bessel values at the top."""
    old = ctx.prec; ctx.prec = BPREC
    x = arb(x)
    c = (arb.pi()/(2*x)).sqrt()
    top1 = c*x.bessel_j(arb(nmax+1)+arb(1)/2)
    top0 = c*x.bessel_j(arb(nmax)+arb(1)/2)
    out = [None]*(nmax+2); out[nmax+1] = top1; out[nmax] = top0
    for n in range(nmax, 0, -1):
        out[n-1] = (2*n+1)/x*out[n] - out[n+1]
    ctx.prec = old
    return [v + 0 for v in out]  # round to working precision (radius kept)

def build(L, R, N, per_unit, chunk=500, verbose=True):
    a = L/2
    pr = primes_upto(L)
    P = sum((2*c for _, c in pr), arb(0))
    eps = phi(R) - P
    # nodes on [0,R]: panels of length 1 (last partial)
    gl = weil.gauss(per_unit)
    nodes = []
    k = 0
    Rf = float(R)
    while k < Rf:
        lo = arb(k); hi = arb(min(k+1, Rf)) if k+1 < Rf else R
        hh = (hi - lo)/2; mid = (hi + lo)/2
        for (x, w) in gl: nodes.append((mid + hh*x, w*hh))
        k += 1
    if verbose: print(f"L={float(L)} R={float(R)} primes={[n for n,_ in pr]} P={P.str(6,radius=False)} eps={eps.str(6,radius=False)} N={N} nodes={len(nodes)}", flush=True)
    K = [arb_mat(N//2 + N%2, N//2 + N%2), arb_mat(N//2, N//2)]
    t0 = time.time()
    for c0 in range(0, len(nodes), chunk):
        ch = nodes[c0:c0+chunk]; m = len(ch)
        Be = arb_mat(K[0].nrows(), m); Bo = arb_mat(K[1].nrows(), m)
        We = arb_mat(K[0].nrows(), m); Wo = arb_mat(K[1].nrows(), m)
        for q, (r, w) in enumerate(ch):
            wt = phi(r) - eps
            for n, c in pr: wt -= 2*c*(r*arb(n).log()).cos()
            wt *= w
            js = sbessel_all(N, a*r)
            for j in range(N):
                s = (-1)**(j//2) * js[j]
                if j % 2 == 0: Be[j//2, q] = s; We[j//2, q] = s*wt
                else: Bo[j//2, q] = s; Wo[j//2, q] = s*wt
        K[0] += We*Be.transpose(); K[1] += Wo*Bo.transpose()
        if verbose: print(f"  chunk {c0//chunk} done {time.time()-t0:.0f}s", flush=True)
    sc = [math.sqrt(2*j+1) for j in range(N)]
    # poles
    A = [arb(0)]*N
    for (x, w) in weil.gauss(N + 60):
        u = weil.leg_on(N, x); e = (a*x/2).exp()
        for j in range(N): A[j] += u[j]*e*w*a.sqrt()
    Ms = []
    for p in (0, 1):
        idx = list(range(p, N, 2)); m = len(idx)
        M = arb_mat(m, m)
        for i, j in enumerate(idx):
            for k2, l in enumerate(idx):
                M[i, k2] = (2*a/arb.pi())*(arb(2*j+1)*(2*l+1)).sqrt()*K[p][i, k2] + (2 if p == 0 else -2)*A[j]*A[l]
            M[i, i] += eps
        Ms.append(M)
    return Ms, eps

if __name__ == "__main__":
    prec = int(sys.argv[1]); ctx.prec = prec
    L = weil.parseL(sys.argv[2]); R = arb(sys.argv[3]); N = int(sys.argv[4]); pu = int(sys.argv[5])
    Ms, eps = build(L, R, N, pu)
    for p, M in enumerate(Ms):
        print(("even" if p == 0 else "odd"), "lam_min(M_R) =", weil.lam_min(weil.sub(M, list(range(M.nrows())))).str(8, radius=False), flush=True)
