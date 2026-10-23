# d4_consts.py -- own computation (ball arithmetic) of the cell-independent constants of the Barban-Vehov mean-square bound
# (Padhi Thm 6.5 / App. A; NOTE.md section 3.4).  Output out/d4_consts.json (arb strings).
#   A*(e1, e2, e12): upper bound of |A(w, s1, s2)| = prod_p |A_p| when the real parts of 1+w+s1, 1+w+s2, 1+w+s1+s2 are e1, e2, e12
#   (each possibly an interval): p <= 100 by branch-and-bound over the 3-torus of phases (independent phases), 100 < p <= 1e6 by
#   log|A_p| <= |A_p - 1| <= (r3 r1 + r3 r2 + r3^2 + r1 r2)/((1-r1)(1-r2)), p > 1e6 by pi(t) < 1.25506 t/log t (Rosser-Schoenfeld).
import math, json, sys, time, heapq
from flint import arb, acb, ctx
ctx.prec = 64
F = lambda v: arb(float(v))
eta, T2, W0, c, cp, s1, th1, b, th2, kap1 = [F(v) for v in (0.02, 400.0, 40.0, 0.02, 0.01, 0.01, 0.40, 0.38, 0.30, 0.38)]
PI = arb.pi(); one = arb(1)
def primes(n):
    sv = bytearray([1])*(n + 1); sv[0:2] = b'\x00\x00'
    for i in range(2, int(n**0.5) + 1):
        if sv[i]: sv[i*i::i] = bytearray(len(sv[i*i::i]))
    return [i for i in range(n + 1) if sv[i]]
PR = primes(10**6)
def Ap_abs_upper(r1, r2, r3, P1, P2, P3):
    x1 = r1*acb(0, P1).exp(); x2 = r2*acb(0, P2).exp(); z = r3*acb(0, P3).exp()
    num = abs((1 - x1 - x2 + z)*(1 - z)); den = abs((1 - x1)*(1 - x2))
    return num.upper()/den.lower() if den.lower() > 0 else float('inf')
def _AdA(r1, r2, r3, P1, P2, P3):
    """A_p and its partial derivatives in the three phases (x1 = r1 e^{iP1}, x2 = r2 e^{iP2}, z = r3 e^{iP3})."""
    I = acb(0, 1)
    x1 = r1*acb(0, P1).exp(); x2 = r2*acb(0, P2).exp(); z = r3*acb(0, P3).exp()
    N = (1 - x1 - x2 + z)*(1 - z); D = (1 - x1)*(1 - x2)
    A = N/D
    dN1 = -I*x1*(1 - z); dD1 = -I*x1*(1 - x2)
    dN2 = -I*x2*(1 - z); dD2 = -I*x2*(1 - x1)
    dN3 = I*z*(x1 + x2 - 2*z)
    dA1 = (dN1*D - N*dD1)/(D*D); dA2 = (dN2*D - N*dD2)/(D*D); dA3 = dN3/D
    return A, (dA1, dA2, dA3)
def sup_Ap(r1, r2, r3, tol=2e-4, maxit=60000, wmin=1e-4):
    """upper bound for sup over phases of |A_p|; r's are arb (balls allowed).  Box bound: g = |A|^2 <= g(centre) + sum_j r_j sup_box|dg/dphi_j|
    (mean value theorem; dg/dphi_j = 2 Re(conj(A) dA/dphi_j) evaluated in ball arithmetic on the box)."""
    rc = [arb(r.mid()) for r in (r1, r2, r3)]
    def ev(box):
        C = [arb((lo + hi)/2) for lo, hi in box]
        B = [arb((lo + hi)/2, (hi - lo)/2) for lo, hi in box]
        Ac, _ = _AdA(r1, r2, r3, *C)
        Ab, dAb = _AdA(r1, r2, r3, *B)
        g = (Ac.real**2 + Ac.imag**2).upper()
        acc = arb(g)
        for j in range(3):
            dg = 2*(Ab.real*dAb[j].real + Ab.imag*dAb[j].imag)
            acc += arb((abs(dg)*((box[j][1] - box[j][0])/2)).upper())
        v = float(acc.sqrt().upper())
        return v if v == v else float('inf')          # NaN (ball contains a pole of the enclosure) -> +inf: box gets split
    def evc(box):
        P = [arb((lo + hi)/2) for lo, hi in box]
        x1 = rc[0]*acb(0, P[0]).exp(); x2 = rc[1]*acb(0, P[1]).exp(); z = rc[2]*acb(0, P[2]).exp()
        return float((abs((1 - x1 - x2 + z)*(1 - z))/abs((1 - x1)*(1 - x2))).lower())
    pi = math.pi; heap = []; best = 0.0
    n1, n2, n3 = 4, 8, 8
    for i in range(n1):
        for j in range(n2):
            for k in range(n3):
                box = ((pi*i/n1, pi*(i + 1)/n1), (-pi + 2*pi*j/n2, -pi + 2*pi*(j + 1)/n2), (-pi + 2*pi*k/n3, -pi + 2*pi*(k + 1)/n3))
                heapq.heappush(heap, (-ev(box), box)); best = max(best, evc(box))
    it = 0
    while True:
        u, box = heapq.heappop(heap); u = -u; it += 1
        w = max(hi - lo for lo, hi in box)
        if u <= best*(1 + tol) or w < wmin or it > maxit:     # u = max of all remaining upper bounds: always a valid bound
            assert u >= best and u < float('inf'), (u, best)
            return u, best, it
        d = max(range(3), key=lambda q: box[q][1] - box[q][0]); lo, hi = box[d]; mid = (lo + hi)/2
        for half in ((lo, mid), (mid, hi)):
            nb = list(box); nb[d] = half; nb = tuple(nb)
            heapq.heappush(heap, (-ev(nb), nb)); best = max(best, evc(nb))
def interval(lo, hi): lo, hi = arb(lo), arb(hi); return (lo + hi)/2 + arb(0, ((hi - lo)/2).upper())
def Astar(e1, e2, e12, P=10**6, verbose=False, crude=False, tol=2e-4, wmin=1e-4):
    """e's: arb balls (intervals allowed).  Returns arb upper bound of sup |A|.  crude: |A_p| <= (1+r1+r2+r3)(1+r3)/((1-r1)(1-r2))
    with the largest moduli (valid for all phases and all real parts in the ranges)."""
    logsum = arb(0); info = []
    for p in PR:
        if p > 100: break
        pa = arb(p)
        if crude:
            r1, r2, r3 = [pa**(-arb(e.lower())) for e in (e1, e2, e12)]
            u = float(((1 + r1 + r2 + r3)*(1 + r3)/((1 - r1)*(1 - r2))).upper()); lo = 0; it = 0
        else:
            u, lo, it = sup_Ap(pa**(-e1), pa**(-e2), pa**(-e12), tol=tol, wmin=wmin)
        logsum += arb(u).log(); info.append((p, round(u, 6), round(lo, 6), it))
    e1l, e2l, e12l = [arb(e.lower()) for e in (e1, e2, e12)]       # moduli are decreasing in the exponents
    acc = arb(0)
    for p in PR:
        if p <= 100: continue
        pa = arb(p); x1, x2, z = pa**(-e1l), pa**(-e2l), pa**(-e12l)
        acc += (z*x1 + z*x2 + z*z + x1*x2)/((1 - x1)*(1 - x2))
    m = min(float((e12l + e1l).lower()), float((e12l + e2l).lower()), float((2*e12l).lower()), float((e1l + e2l).lower()))
    m = arb(m)
    Pa = arb(P)
    tail = 4*arb('1.25506')*m*Pa**(1 - m)/((m - 1)*Pa.log())/((1 - Pa**(-e1l))*(1 - Pa**(-e2l)))
    tot = (logsum + acc + tail).exp()
    if verbose: print('   p<=100 log', logsum.str(6), ' 100<p<=1e6', acc.str(6), ' tail', tail.str(6), ' worst primes', info[:4], flush=True)
    return arb(tot.upper())
def Astar_pieces(e1lo, e1hi, e2fun, e12fun, npieces, verbose=False):
    """sup over e1 in [e1lo, e1hi] (pieces), e2 = e2fun(e1-ball), e12 = e12fun(e1-ball)."""
    best = None
    for j in range(npieces):
        a = e1lo + (e1hi - e1lo)*j/npieces; bb = e1lo + (e1hi - e1lo)*(j + 1)/npieces
        E1 = interval(a, bb)
        v = Astar(E1, e2fun(E1), e12fun(E1))
        if verbose: print('  piece', j, v.str(6), flush=True)
        if best is None or v.upper() > best.upper(): best = v
    return best
def zeta(s): return arb(s).zeta()
def iota(s): s = arb(s); return s.zeta()/(2*s).zeta()
def gam_int(sig, weight, n=40000):
    """int_{|v| <= W0} |Gamma(sig + iv)| (1 + |v|)^weight dv: upper Riemann sum (|Gamma| decreasing in v (product formula))."""
    h = W0/n; tot = arb(0)
    for i in range(n):
        g = acb(sig, i*h).gamma().abs_upper()
        tot += arb(g)*((1 + (i + 1)*h) if weight else 1)
    return 2*tot*h
def Jint(sig, lo, hi, weight):
    """2 int_lo^hi e^{eta^2 (sig^2 - t^2)/2} (1+t)^weight/(sig^2 + t^2) dt (upper sum; tail beyond 3000 explicit)."""
    top = min(hi, 3000.0); tot = arb(0); t0 = lo; E2 = eta*eta/2
    while t0 < top:
        hh = 1e-4 if t0 < 0.1 else (1e-3 if t0 < 2 else 0.01)
        t1 = min(t0 + hh, top); a0, a1 = arb(t0), arb(t1)
        tot += (E2*(sig*sig - a0*a0)).exp()*((1 + a1) if weight else 1)/(sig*sig + a0*a0)*(a1 - a0)
        t0 = t1
    if hi > 3000:
        a = arb(3000)
        tot += (E2*sig*sig).exp()*(arb(2)/a if weight else 1/(a*a))*(-(E2*a*a)).exp()/(eta*eta*a)
    return 2*tot
if __name__ == "__main__":
    import os
    t0 = time.time(); OUTF = 'out/d4_consts.json'
    out = json.load(open(OUTF)) if os.path.exists(OUTF) else {}          # merge: never drop constants of earlier runs
    def put(k, v): out[k] = v.str(30, radius=True); print(k, v.str(10), '%.0fs' % (time.time() - t0), flush=True); json.dump(out, open(OUTF, 'w'), indent=1)
    W = set((sys.argv[1] if len(sys.argv) > 1 else 'all').split(','))
    if W & {'all', 'A'}:
        sig2 = 1 + c - th1 - cp
        put('A1', Astar(1 + c - th1, 1 + c - cp, sig2, verbose=True))
        put('Ax', Astar(1 + cp, 1 + c - cp, one, verbose=True))
        put('AV', Astar(1 - b, 1 + c + b, one, verbose=True))
        put('AU', Astar(1 + b, 1 + c - b, one, verbose=True))
    if W & {'all', 'A', 'AK'}:
        # AK: sup over e2 in [1 + kap1 - th2, 1 + kap1 + c].  For fixed x1, z the factor A_p is analytic in x2 on |x2| < 1 (only pole
        # x2 = 1), so by maximum modulus sup_{|x2| <= r} |A_p| = sup_{|x2| = r}: the bound is attained at the largest modulus
        # r2 = p^{-e2}, i.e. at the smallest e2.  One evaluation at e2 = 1 + kap1 - th2 covers the whole range.
        put('AK', Astar(1 - kap1, 1 + kap1 - th2, one, verbose=True))
    if W & {'all', 'A', 'boxes'}:
        # boxes along contours that only enter terms carrying e^{-eta^2 T2^2/2} or the tiny R2: crude bound
        put('A1box', Astar(interval(1 + c - th1, 1 + c + s1), 1 + c - cp, interval(1 + c - th1 - cp, 1 + c + s1 - cp), crude=True, verbose=True))
        put('AVbox', Astar(interval(1 - b, 1 + cp), interval(1 + c - cp, 1 + c + b), one, crude=True))
        put('AUbox', Astar(interval(1 + cp, 1 + b), interval(1 + c - b, 1 + c - cp), one, crude=True))
        put('AK2box', Astar(interval(1 - kap1, 1 + cp), interval(1 - th2 - cp, 1 + c + kap1), one, crude=True))
    if W & {'all', 'rest'}:
        put('gam_c_w1', gam_int(c, True)); put('gam_c', gam_int(c, False))
        put('J_mth1_0_T2_w1', Jint(-th1, 0.0, 400.0, True)); put('J_mcp_0_inf_w1', Jint(-cp, 0.0, math.inf, True))
        put('J_mcp_0_inf', Jint(-cp, 0.0, math.inf, False)); put('J_s1_T2_inf_w1', Jint(s1, 400.0, math.inf, True))
        for nm, s in [('z1c', 1 + c), ('z1c_s1_cp', 1 + c + s1 - cp), ('z2_2c', 2 + 2*c)]: put(nm, zeta(s))
        for nm, s in [('iz1c_cp', 1 + c - cp), ('iz1c_s1', 1 + c + s1), ('iz1cp', 1 + cp), ('iz1c_b', 1 + c + b), ('iz1b', 1 + b), ('iz1k_t2', 1 + kap1 - th2), ('iz13', arb('1.3'))]:
            put(nm, iota(s))
