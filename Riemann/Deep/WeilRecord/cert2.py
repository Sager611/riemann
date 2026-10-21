# pyright: reportMissingImports=false
"""Pointwise-envelope Weil certificate, v2 (fixes the window-1.875 assembly failure of ../weil/cert.py).

Same proof chain as ../weil/cert.py (steps (1)-(6) there); changes:
 * quadrature: panels [0,1/2],[1/2,1],[1,2],[2,4],...,then length lmax, with per-panel Bernstein ellipse rho
   and node count chosen from a RIGOROUS bound: the normalised integrand
   (2a/pi) sqrt((2j+1)(2k+1)) w(z) j_j(az) j_k(az) is bounded on the ellipse by (2a/pi) e^{2aB} max|w|,
   using |j_n(z)| <= e^{|Im z|}/sqrt(2n+1) (Poisson integral, Cauchy-Schwarz on int|P_n|), and max|w| from
   ball evaluation of the digamma extension on boxes covering the ellipse boundary (maximum modulus);
   the ellipse avoids the digamma poles z = +-i(2k+1/2) (checked).  ~2 nodes per unit length instead of 64.
 * assembly: Bessel values below 2^-400 are replaced by EXACT zeros and the dropped mass is added to the
   error budget (the old code used zero-centred balls with radii spanning ~10^4 decades; arb_mat block
   multiplication then returned infinite radii at window 1.875).  Every chunk is checked for finiteness.
 * node rounding to exact points: Cauchy estimate |F'| <= M/(A-h) on the real panel (A-h = distance from
   the panel to the ellipse boundary).
"""
import sys, time, math
import os
_H = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _H)
for _d in (os.path.join(_H, '..', 'Weil', 'cert'), os.path.join(_H, '..', 'weil')):
    if os.path.isdir(_d): sys.path.insert(1, _d)   # weil.py, cert.py, envelope.py
from flint import arb, arb_mat, acb, ctx
import weil, envelope, quad
from cert import phi, dfact, certified_shift

TINY_EXP = 400

def phi_ext(z):
    """holomorphic extension of Re psi(1/4+ir/2) - log pi"""
    return (acb(arb(1)/4) + acb(0, 1)*z/2).digamma()/2 + (acb(arb(1)/4) - acb(0, 1)*z/2).digamma()/2 - arb.pi().log()

def bessel_start(n, x, bprec, x_fn=None):
    """j_n(x) = sqrt(pi/(2x)) J_{n+1/2}(x) with relative accuracy ~2^-(bprec-32).  ROOT CAUSE of the old window-1.875
    failure: arb's J_nu(x) for large nu and x (series with cancellation ~e^x) silently returns balls of radius up to
    ~1e+97 at the old fixed precision; the downward recurrence then produced radii beyond double range ("inf").
    Here the working precision is raised until the start value is accurate."""
    p = bprec; old = ctx.prec
    while True:
        ctx.prec = p
        if x_fn is not None: x = x_fn()          # recompute x sharply at the raised precision
        v = (arb.pi()/(2*x)).sqrt()*x.bessel_j(arb(n) + arb(1)/2)
        ctx.prec = old
        if v != 0 and v.rel_accuracy_bits() >= bprec - 40: return v
        if p > 40000: raise ValueError(f"bessel start value not accurate: n={n}")
        p = int(p*1.5)

def sbessel(nmax, a, r, bprec, a_fn=None):
    """j_n(a r), n < nmax.  a_fn (optional): callable returning a at the current precision (a irrational)."""
    old = ctx.prec; ctx.prec = bprec
    x = a*r
    assert x.rad() == 0 or float(x.rad()) <= 2.0**(-bprec + 16)*float(abs(x).upper())
    x_fn = (lambda: a_fn()*r) if a_fn is not None else None
    out = [None]*(nmax+2)
    out[nmax+1] = bessel_start(nmax+1, x, bprec, x_fn)
    out[nmax] = bessel_start(nmax, x, bprec, x_fn)
    ix = 1/x
    for n in range(nmax, 0, -1):
        out[n-1] = (2*n+1)*ix*out[n] - out[n+1]
    ctx.prec = old
    tiny = arb(2)**(-TINY_EXP)
    radcap = arb(2)**(-old + 64)
    res = []
    for v in out[:nmax]:
        if abs(v).upper() < tiny: res.append(arb(0))
        else:
            v = v + 0
            if not (v.is_finite() and v.rad() < radcap): raise ValueError("inaccurate Bessel value (radius too large)")
            res.append(v)
    return res

def plan_panels(a, R, eps, pr, target, lmax, log):
    """Choose (rho, n) per panel so that the normalised entry error is <= target/len(panels)."""
    P = sum((2*c for _, c in pr), arb(0))
    pans = quad.panels(R, lmax)
    per = target/len(pans)
    plan = []; tot = arb(0); totP = arb(0)
    for (lo, hi) in pans:
        best = None
        for rho in ((1.5, 2) if float(hi) <= 2 else (2, 3, 4)):
            if not quad.digamma_poles_clear(lo, hi, rho): continue
            c, h, A, B = quad.ellipse_geom(lo, hi, rho)
            mphi = quad.max_abs_on_ellipse(lambda z: phi_ext(z) - eps, lo, hi, rho, nbox=32)
            if mphi is None: continue
            mw = mphi + sum((2*cc*(B*arb(n).log()).cosh() for n, cc in pr), arb(0))
            M = (2*a/arb.pi())*(2*a*B).exp()*mw
            n = 4
            while not (quad.gauss_err_factor(h, rho, n)*M < per): n += 4
            if best is None or n < best[1]: best = (rho, n, M, A - h, h)
        if best is None: raise ValueError("no admissible ellipse")
        rho, n, M, dmin, h = best
        plan.append((lo, hi, n, rho, M, dmin))
        tot += quad.gauss_err_factor(h, rho, n)*M
    return plan, tot

def certify(L, R, N, deltas, lmax=64, target=arb(10)**-60, chunk=500, log=print):
    L = arb(L); R = arb(R)
    a = L/2; x_max = a*R
    pr = envelope.primes_upto(L)
    P = arb(0)
    for _, c in pr: P += 2*c
    eps = phi(R) - P
    assert eps > 0, "R too small"
    bprec = ctx.prec + int(float(x_max)) + 128          # recurrence radius growth ~0.68 bit per unit x; start values adaptive
    log(f"L={L.str(8, radius=False)} a={float(a):.5f} R={float(R)} primes={[n for n,_ in pr]} P={P.str(8,radius=False)} eps={eps.str(8,radius=False)} N={N} prec={ctx.prec} bprec={bprec}")
    # (4) Legendre tail
    phi0 = phi(arb(0))
    supw = arb(max(abs(phi0).upper(), phi(R).upper())) + P + eps
    supw = arb(supw.upper())
    normA = (2*a.sinh()).sqrt()
    terms = []; df = dfact(N)
    for j in range(N, N + 400):
        if j > N: df *= (2*j+1)
        ej = ((2/arb.pi())*x_max**(2*j+1)).sqrt()/df
        Aj = (2*(2*j+1)*a).sqrt()*(a/2)**j*(a/2).exp()/df
        terms.append((supw*ej + 2*Aj*normA)**2)
    rat = terms[-1]/terms[-2]; assert rat < 1
    Ebound = 2*(sum(terms, arb(0)) + terms[-1]*rat/(1 - rat)).sqrt()
    log(f"  tail: ||E|| <= {Ebound.str(5, radius=False)}")
    # (5) quadrature plan + bound
    t0 = time.time()
    plan, qerr = plan_panels(a, R, eps, pr, target, lmax, log)
    nn = sum(p[2] for p in plan)
    half = (N + 1)//2
    Dq = arb(half)*qerr                                 # ||Delta||_F <= dim * max entry error
    log(f"  quadrature: {len(plan)} panels, {nn} nodes, entry err <= {qerr.str(3,radius=False)}, ||Delta|| <= {Dq.str(3,radius=False)} ({time.time()-t0:.0f}s)")
    # node rounding: on the real line |j_n|,|j_n'| <= 1/sqrt(2n+1) (Poisson integral), |psi'(1/4+ir/2)| <= psi'(1/4) < 18,
    # so the normalised integrand has |F'| <= (2a/pi) (9 + sum 2c log n + 2a supw).
    dF = (2*a/arb.pi())*(9 + sum((2*c*arb(n).log() for n, c in pr), arb(0)) + 2*a*supw)
    nodes = []; perr = arb(0)
    for (lo, hi, n, rho, M, dmin) in plan:
        hh = (hi - lo)/2; md = (hi + lo)/2; mr = 0.0
        for (x, w) in weil.gauss(n):
            rr = md + hh*x
            nodes.append((arb(rr.mid()), w*hh)); mr = max(mr, float(rr.rad()))
        perr += 2*hh*arb(mr)*dF
    flush_err = (2*a/arb.pi())*(2*N)*R*supw*3*arb(2)**(-TINY_EXP)
    Dbound = Dq + arb(half)*(perr + flush_err)
    log(f"  node rounding + flush: ||Delta|| <= {Dbound.str(3, radius=False)}")
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
        for p in (0, 1):
            if not all(K[p][i, i].is_finite() for i in range(K[p].nrows())): raise ValueError("non-finite K")
        log(f"    nodes {c0+m}/{len(nodes)} {time.time()-t0:.0f}s")
    Av = [(2*(2*j+1)*a).sqrt()*(arb.pi()/a).sqrt()*(a/2).bessel_i(arb(j)+arb(1)/2) for j in range(N)]
    assert all(v.is_finite() for v in Av)
    res = {}
    for p in (0, 1):
        idx = list(range(p, N, 2)); m = len(idx)
        G = arb_mat(m, m)
        for i, j in enumerate(idx):
            for k2, l in enumerate(idx):
                G[i, k2] = (2*a/arb.pi())*(arb(2*j+1)*(2*l+1)).sqrt()*K[p][i, k2] + (2 if p == 0 else -2)*Av[j]*Av[l]
            G[i, i] += eps
        rad = max(float(G[i, k2].rad()) for i in range(m) for k2 in range(m))
        if not math.isfinite(rad): raise ValueError("non-finite G radius")
        log(f"  {'even' if p==0 else 'odd'}: max entry radius {rad:.3e}, approx lam_min(G) = {weil.lam_min(G.mid()).str(6, radius=False)}")
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
                ok = d*(1 - 1e-3); log(f"  {'even' if p==0 else 'odd'}: lam_min >= {d}(1-1e-3) [residual+radius {b.str(3,radius=False)}, {time.time()-t1:.0f}s]"); break
            log(f"  {'even' if p==0 else 'odd'}: PD NOT certified at delta={d}")
        res[p] = ok
    r0, r1 = res[0], res[1]
    if r0 is None or r1 is None or not (eps > min(r0, r1) + Ebound):
        log(f"RESULT L={L.str(8,radius=False)}: not certified"); return False
    log(f"RESULT L={L.str(8,radius=False)}: CERTIFIED Q(f) >= {min(r0, r1)} ||f||^2 for all f supported in [-L/2, L/2]")
    return True

if __name__ == "__main__":
    ctx.prec = int(sys.argv[1])
    L = weil.parseL(sys.argv[2]); R = arb(sys.argv[3]); N = int(sys.argv[4])
    deltas = [float(s) for s in sys.argv[5].split(',')]
    lmax = int(sys.argv[6]) if len(sys.argv) > 6 else 64
    certify(L, R, N, deltas, lmax=lmax, log=lambda s: print(s, flush=True))
