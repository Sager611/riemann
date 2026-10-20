# pyright: reportMissingImports=false
"""Rigorous (Arb ball arithmetic) refutation of Desogus, arXiv 2609.20367v3, Theorem 8.6 / ineq. (201).

Claim (201)+(MASTER-SCALAR+): if A_k > 0 (k >= 7) then for every f in the collar K_k,
   <H_k f, A_{k+1} H_k f> >= g_k |<e_k,f>|^2,  g_k >= 1.878  (Appendix A.5, eq. (210)).
Since H_k f minimises <v,A_{k+1}v> over v = (old part) + f, this implies for EVERY odd v supported in
(-a_{k+1}, a_{k+1}):   W(v) = <v, A_{k+1} v> >= g_k |<e_k, v>|^2.
Here W is the Weil functional (pole term -2<sinh(y/2),v>^2 on odd v), e_k = normalised sgn(y) e^{-|y|/2}
on the collar a_k < |y| < a_{k+1}, a_k = log(k)/2 (paper sec. 6.5, Lemma 6.4/6.8: q_Y(t)=t^{-1/2}).
We exhibit one explicit odd polynomial v (exact rational coefficients) and enclose W(v) and <e_k,v> rigorously.
W(v) = g(0)(-gamma - log pi - log(1-e^{-4a})) + int_0^{2a} [2e^{-2t} g(0) - 2 g(t) e^{-t/2}]/(1-e^{-2t}) dt
       - 2 (int v e^{y/2})^2 - 2 sum_{n<e^{2a}} Lambda(n) n^{-1/2} g(log n),   g(t) = int v(y) v(y-t) dy
(same convention as Riemann/Deep/Weil/cert/weil.py, checked against zeta zeros to 1e-28)."""
import sys, time
sys.path.insert(0, '/Users/adriansager/playground/riemann/.scratch/weil')
import flint
from flint import arb, acb, arb_mat, fmpq, fmpq_poly, acb_poly, ctx
import weil
from refute import coeffs

def galerkin_v(k, N, prec):
    ctx.prec = prec
    L = arb(k+1).log(); a = L/2; ak = arb(k).log()/2
    M = weil.weil_matrix(L, N)
    odd = list(range(1, N, 2)); Mo = weil.sub(M, odd)
    c = coeffs(N, a, ak, a, lambda y: (1 if y > 0 else -1)*(-abs(y)/2).exp())
    cv = arb_mat(len(odd), 1)
    for i, j in enumerate(odd): cv[i, 0] = c[j]
    x = Mo.solve(cv)
    est = (x.transpose()*Mo*x)[0, 0].mid()
    xs = [x[i, 0].mid() for i in range(len(odd))]
    mx = max(abs(float(t)) for t in xs)
    P = fmpq_poly([])
    for i, j in enumerate(odd):
        s = xs[i]*((arb(2*j+1)/2).sqrt())/a.sqrt()/mx       # v(y) = P(y/a)
        q = fmpq(int((s.mid()*arb(2)**400).floor().unique_fmpz()), 2**400)
        P += fmpq_poly.legendre_p(j)*q
    return P, est, mx

def G_poly(P):
    """G(s) = int_{s-1}^{1} P(x) P(x-s) dx, exact, via P(x-s) = sum_i P^(i)(x) (-s)^i / i!."""
    s = fmpq_poly([0, 1]); sm1 = fmpq_poly([-1, 1])
    G = fmpq_poly([]); D = P; fact = 1; i = 0
    while D != 0:
        R = (P*D).integral()
        G += (fmpq_poly([R(fmpq(1))]) - R(sm1)) * s**i * fmpq((-1)**i, fact)
        i += 1; fact *= i; D = D.derivative()
    return G

def to_acb(p):
    return acb_poly([acb(arb(fmpq(c))) for c in p.coeffs()])

def run(k, N, prec_g=384, prec=3000):
    t0 = time.time()
    P, est, mx = galerkin_v(k, N, prec_g)
    G = G_poly(P)
    ctx.prec = prec
    a = arb(k+1).log()/2; ak = arb(k).log()/2
    Gp = to_acb(G); Pp = to_acb(P); dG = to_acb(G.derivative())
    g = lambda t: a*Gp(t/a)
    g0 = a*arb(fmpq(G(fmpq(0))))           # = ||v||^2
    # archimedean integral, split [0,d] u [d,2a]
    d = arb(2)**-300
    f = lambda t, _: (2*(-2*t).exp()*g0 - 2*g(t)*(-t/2).exp()) / (1 - (-2*t).exp())
    I1 = acb.integral(f, d, 2*a)
    tb = arb(d/2, d/2)                        # ball [0,d]
    Nd = -4*(-2*tb).exp()*g0 - 2*dG(acb(tb)/a)*(-tb/2).exp() + g(acb(tb))*(-tb/2).exp()
    I0 = arb(0, (d*abs(Nd)).upper())          # |N(t)| <= t sup|N'|, 1-e^{-2t} >= t on [0,1/2]
    arch = g0*(-(arb.const_euler() + arb.pi().log() + (1 - (-4*a).exp()).log())) + I1.real + I0
    A = (acb.integral(lambda x, _: Pp(x)*(a*x/2).exp(), -1, 1)*a).real
    pole = -2*A*A
    primes = arb(0)
    for n in range(2, k+1):
        lam = weil.vm(n)
        if lam is not None:
            primes += lam/arb(n).sqrt()*2*g(acb(arb(n).log())).real
    W = arch + pole - primes
    ev = (acb.integral(lambda x, _: Pp(x)*(-a*x/2).exp(), ak/a, 1)*2*a).real
    ne2 = 2*((-ak).exp() - (-a).exp())      # ||sgn e^{-|y|/2}||^2 on collar
    ratio = W*ne2/(ev*ev)
    print(f"k={k} N={N}: deg v={P.degree()}, W(v) in {W.str(8)}, Galerkin est {float(est/mx**2):.8e}", flush=True)
    print(f"  <e_k,v> in {(ev/ne2.sqrt()).str(8)}; ||v||^2 in {g0.str(8)}", flush=True)
    print(f"  W(v)/|<e_k,v>|^2 in {ratio.str(8)}  (claim: >= g_k >= 1.878); upper < 1.878: {bool(ratio < arb('1.878'))}", flush=True)
    print(f"  t={time.time()-t0:.1f}s", flush=True)

if __name__ == '__main__':
    N = int(sys.argv[1])
    for k in sys.argv[2:]: run(int(k), N)
