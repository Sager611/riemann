# certlow.py (morirate): Chekmachev's certificate (Theorem 5.10 + Definition 6.1; ChekAudit's re-implementation cert.py) for the
# DELETED form Q' = Q - 2<f, cos(gamma0 x)>^2 (even sector): the deleted zero is one more negative rank-one column
# ('eval', gamma0, sqrt 2) in B, so Q' - eps >= D + U M U* - B'B'* >= 0 iff lmax(B'* P^-1 B') <= 1 (Woodbury), as for Q.
# gamma0 is an Arb ball certified to contain the first zero (sign change of Hardy's Z). Bisection on eps.
# Usage: certlow.py c bits eps_lo eps_hi [steps]
import sys, time, math
from flint import arb, acb, arb_mat, acb_mat, ctx

def vm(n):
    for p in range(2, n+1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return arb(p).log() if m == 1 else None
    return None

class Sym:
    def __init__(s, c):
        s.c = c; s.L = arb(c).log(); s.pi = arb.pi(); s.lnpi = s.pi.log()
        s.atoms = []
        for n in range(2, c):
            lam = vm(n)
            if lam is not None: s.atoms.append((arb(n).log(), lam/arb(n).sqrt()))
        s.sig2 = 2*sum((w for _, w in s.atoms), arb(0))
        # tail coefficients c^{-a_n}, a_n = 2n+1/2, enough terms + geometric remainder
        s.tail = []
        eps = arb(2)**(-(ctx.prec+30)); n = 0
        while True:
            an = arb(2*n) + arb(1)/2; t = arb(c)**(-an)
            s.tail.append((an, t))
            if (t < arb(2)**(-(ctx.prec+30))) and n > 2: break
            n += 1
        an = arb(2*len(s.tail)) + arb(1)/2
        q = arb(c)**(-2)
        s.trem = arb(c)**(-an)*an/(1-q)*2   # bound for sum_{n>=N} c^{-a_n} a_n (an grows, ratio <= q*(1+2/an) -> use factor 2 safety)
    def r(s, k): return s.pi*k/s.L
    def hplus(s, r):
        return acb(arb(1)/4, r/2).digamma().real - s.lnpi
    def T(s, k):
        r = s.r(k); r2 = r*r
        tot = arb(0)
        for an, t in s.tail: tot += t*an/(an*an + r2)
        rem = s.trem/(r2 + arb(1)/16)
        tot = tot + arb(0, rem.upper())
        return -2*(1 if k % 2 == 0 else -1)*tot
    def mu(s, k):
        r = s.r(k)
        ps = arb(0)
        for ln, w in s.atoms: ps += w*(r*ln).cos()
        return s.hplus(r) - s.T(k) - 2*ps
    def monotone(s, k):
        r = s.r(k)
        TT = arb(0)
        for an, t in s.tail: TT += t*(1/r + an/(r*r))
        TT += s.trem*2/r
        return s.hplus(r) - s.sig2 - 2*TT

def scan_end(S, start, thr):
    # smallest k >= start with monotone(k+1) >= thr (lower endpoint)
    ok = lambda k: S.monotone(k+1).lower() >= thr
    if ok(start): return start
    lo, hi = start, start*2
    while not ok(hi): lo, hi = hi, hi*2
    while hi - lo > 1:
        m = (lo+hi)//2
        if ok(m): hi = m
        else: lo = m
    return hi

def certify(c, sector, epsstr, gball, delete=True):
    eps = arb(epsstr) if epsstr != '0' else arb(0)
    t0 = time.time()
    S = Sym(c); L = S.L; pi = S.pi
    K2 = max(201, math.ceil(5*c*c/3));  K2 += (K2 % 2 == 0)
    J = K2
    tau, taue, t1 = 1, 0, arb(1)
    ks = scan_end(S, K2, tau); kse = scan_end(S, 2*J, taue)
    # rigorous scan (balls at low precision): odd k in (K2, ks], even 2j in (2J, kse]
    selo, sele = [], []
    oldp = ctx.prec; ctx.prec = 80
    S80 = Sym(c)
    for k in range(K2+2, ks+1, 2):
        if S80.mu(k).lower() < tau: selo.append(k)
    j0 = J+1
    for j in range(j0, kse//2 + 1):
        if S80.mu(2*j).lower() < taue: sele.append(j)
    ctx.prec = oldp
    tscan = time.time()-t0
    Kset = set(range(1, K2+1, 2)) | set(selo)
    Eset = set(range(0, J+1)) | set(sele)
    for k in selo:
        for j in range((k-1)//2 - 32, (k-1)//2 + 33):
            if j >= 0: Eset.add(j)
    for j in sele:
        for jj in range(j-3, j+4): Eset.add(jj)
    # symbol balls
    mus = {}
    def mu(k):
        if k not in mus: mus[k] = S.mu(k)
        return mus[k]
    # window rule (2): even 2j > K2 with mu < 0 -> odd |k-2j|<=2 to K, even |j'-j|<=3 to E (iterate)
    changed = True
    done = set()
    while changed:
        changed = False
        for j in sorted(Eset):
            if j in done or 2*j <= K2: continue
            done.add(j)
            if mu(2*j).mid() < 0:
                for k in (2*j-1, 2*j+1):
                    if k not in Kset: Kset.add(k); changed = True
                for jj in range(j-3, j+4):
                    if jj not in Eset: Eset.add(jj); changed = True
    if sector == 'odd': Eset.discard(0)
    Kl = sorted(Kset); El = sorted(Eset)
    w = {k: ((mu(k) - t1)/L).lower() for k in Kl}
    deff = arb(1)/2 - eps
    dl = {j: ((mu(2*j) + t1)/2 - eps).lower() for j in El}
    Kp = [k for k in Kl if w[k] > 0]; Km = [k for k in Kl if not (w[k] > 0)]
    G = [j for j in El if dl[j] > 0]; Jm = [j for j in El if not (dl[j] > 0)]
    # basis normalisation
    def Nj(j): return (1/L).sqrt() if (j == 0 and sector == 'even') else (2/L).sqrt()
    sL4, cL4, sL2 = (L/4).sinh(), (L/4).cosh(), (L/2).sinh()
    def sgk(k): return 1 if ((k-1)//2) % 2 == 0 else -1
    # inner products
    def ip_phi_node(j, k):   # <phi_j, e_k> or <phi^-_j, s_k>, k odd
        sj = 1 if j % 2 == 0 else -1
        if sector == 'even':
            return Nj(j)*(L/pi)*sj*2*k*sgk(k)/(k*k - 4*j*j)
        else:
            return Nj(j)*sj*4*j*(L/pi)*sgk(k)/(k*k - 4*j*j)
    def ip_phi_v(j):
        sj = 1 if j % 2 == 0 else -1
        om = 2*pi*j/L
        if sector == 'even': return Nj(j)*sj*sL4/(om*om + arb(1)/4)
        else: return Nj(j)*(-2)*om*sj*sL4/(om*om + arb(1)/4)
    def ip_node_v(k):
        r = S.r(k)
        if sector == 'even': return 2*r*sgk(k)*cL4/(r*r + arb(1)/4)
        else: return sgk(k)*cL4/(r*r + arb(1)/4)
    vv = sL2 + L/2 if sector == 'even' else sL2 - L/2
    # columns: B = [sqrt(deff-dl_j) phi_j (Jm)] + [sqrt(-w_k) node_k (Km)] (+ sqrt2 v in odd)
    #          U = [node_k (Kp)] (+ v in even), M = diag(w_k, 2)
    Bcols = [('phi', j, (deff - dl[j]).sqrt()) for j in Jm] + [('node', k, (-w[k]).sqrt()) for k in Km]
    Ucols = [('node', k, arb(1)) for k in Kp]
    Mdiag = [w[k] for k in Kp]
    if sector == 'even': Ucols.append(('v', None, arb(1))); Mdiag.append(arb(2))
    else: Bcols.append(('v', None, arb(2).sqrt()))
    if delete:
        assert sector == 'even'
        Bcols.append(('eval', gball, arb(2).sqrt()))
    h = L/2
    def ip_eval_phi(g, j):
        om = 2*pi*j/L
        return Nj(j)*(((om - g)*h).sin()/(om - g) + ((om + g)*h).sin()/(om + g))
    def ip_eval_node(g, k):
        r = S.r(k)
        return ((r - g)*h).sin()/(r - g) + ((r + g)*h).sin()/(r + g)
    def ip_eval_v(g):
        return 2*((h/2).sinh()*(g*h).cos()/2 + g*(h/2).cosh()*(g*h).sin())/(g*g + arb(1)/4)
    def ip_eval_eval(g):
        return h + (2*g*h).sin()/(2*g)
    def ip(a, b):   # <a,b> for a,b column descriptors (unscaled)
        ta, ia, _ = a; tb, ib, _ = b
        if ta == 'eval' or tb == 'eval':
            if ta == 'eval' and tb == 'eval': return ip_eval_eval(ia)
            if tb == 'eval': ta, ia, tb, ib = tb, ib, ta, ia
            if tb == 'phi': return ip_eval_phi(ia, ib)
            if tb == 'node': return ip_eval_node(ia, ib)
            if tb == 'v': return ip_eval_v(ia)
        if ta == 'phi' and tb == 'phi': return arb(1) if ia == ib else arb(0)
        if ta == 'phi' and tb == 'node': return ip_phi_node(ia, ib)
        if ta == 'node' and tb == 'phi': return ip_phi_node(ib, ia)
        if ta == 'node' and tb == 'node': return L/2 if ia == ib else arb(0)
        if ta == 'phi' and tb == 'v': return ip_phi_v(ia)
        if ta == 'v' and tb == 'phi': return ip_phi_v(ib)
        if ta == 'node' and tb == 'v': return ip_node_v(ia)
        if ta == 'v' and tb == 'node': return ip_node_v(ib)
        return vv
    chi = [1/dl[j] - 1/deff for j in G]
    allc = Bcols + Ucols
    nA = len(allc)
    # matrix Phi (|G| x nA) of <phi_j, col>, then Gram = (1/deff) <a,b> + Phi^T diag(chi) Phi
    Phi = arb_mat(len(G), nA); PhiC = arb_mat(len(G), nA)
    for gi, j in enumerate(G):
        for ci, col in enumerate(allc):
            t = ip(('phi', j, None), col)*col[2]
            Phi[gi, ci] = t; PhiC[gi, ci] = t*chi[gi]
    Gram = Phi.transpose()*PhiC
    for a in range(nA):
        for b in range(nA):
            Gram[a, b] = Gram[a, b] + ip(allc[a], allc[b])*allc[a][2]*allc[b][2]/deff
    nB = len(Bcols); nU = len(Ucols)
    X = arb_mat(nB, nB); Y = arb_mat(nU, nB); Sm = arb_mat(nU, nU)
    for a in range(nB):
        for b in range(nB): X[a, b] = Gram[a, b]
    for a in range(nU):
        for b in range(nB): Y[a, b] = Gram[nB+a, b]
        for b in range(nU): Sm[a, b] = Gram[nB+a, nB+b] + ((1/Mdiag[a]) if a == b else 0)
    tasm = time.time()-t0
    Z = Sm.solve(Y)
    W = X - Y.transpose()*Z
    IW = arb_mat(nB, nB)
    for a in range(nB):
        for b in range(nB): IW[a, b] = (1 if a == b else 0) - W[a, b]
    R = arb_mat(nB, nB)
    for a in range(nB):
        for b in range(nB): R[a, b] = IW[a, b].mid()
    ev = sorted([e.real.mid() for e in acb_mat(R).eig(algorithm="approx")], key=lambda z: float(z))
    lmin = ev[0]
    # rigorous: Lemma 5.5 with nu = 0.999*lmin (if lmin>0): Cholesky of midpoint (IW - s0 I), residual in balls
    rig = None
    if lmin > 0:
        nu = lmin*arb('0.999'); s0 = lmin*arb('0.9995')
        n = nB
        A = [[R[i, j] - (s0 if i == j else 0) for j in range(n)] for i in range(n)]
        Lc = [[arb(0)]*n for _ in range(n)]
        okc = True
        for i in range(n):
            for j in range(i+1):
                s = A[i][j]
                for kk in range(j): s -= Lc[i][kk]*Lc[j][kk]
                if i == j:
                    if not (s.mid() > 0): okc = False; break
                    Lc[i][j] = s.mid().sqrt().mid()
                else:
                    Lc[i][j] = (s/Lc[j][j]).mid()
            if not okc: break
        if okc:
            Lm = arb_mat(n, n)
            for i in range(n):
                for j in range(i+1): Lm[i, j] = Lc[i][j]
            Fm = IW - Lm*Lm.transpose()
            emax = arb(0)
            for i in range(n):
                rs = arb(0)
                for j in range(n):
                    v = Fm[i, j] - (s0 if i == j else 0)
                    rs += abs(v).upper()
                if rs.upper() > emax.upper(): emax = rs.upper()
            rig = (nu, bool((emax <= s0 - nu)))
    return lmin, (rig[1] if rig else False), nB, time.time()-t0

def Zfun(t):
    s_ = acb(arb(1)/2, t)
    th = acb(arb(1)/4, t/2).lgamma().imag - t/2*arb.pi().log()
    z = s_.zeta()*acb(th.cos(), th.sin())
    return z.real

def zero_ball(k, prec):
    import mpmath as mp
    mp.mp.prec = prec + 40
    g = mp.zetazero(k).imag
    gs = mp.nstr(g, int(prec*0.30) + 5)
    r = arb(10)**(-int(prec*0.30) + 20)
    zl = Zfun(arb(gs) - r); zr = Zfun(arb(gs) + r)
    ok = bool((zl < 0 and zr > 0) or (zl > 0 and zr < 0))
    return arb(gs).union(arb(gs) - r).union(arb(gs) + r), ok, float(g)

if __name__ == '__main__':
    import json
    c = int(sys.argv[1]); ctx.prec = int(sys.argv[2]); lo = float(sys.argv[3]); hi = float(sys.argv[4])
    steps = int(sys.argv[5]) if len(sys.argv) > 5 else 14
    gball, ok, gf = zero_ball(1, ctx.prec)
    print(json.dumps({'c': c, 'prec': ctx.prec, 'gamma1': gf, 'sign_change_certified': ok}), flush=True)
    m0, r0, nB, tt = certify(c, 'even', '0', gball)
    print(json.dumps({'eps': 0, 'margin': m0.str(8, radius=False), 'rigorous': r0, 'order': nB, 't': round(tt, 1)}), flush=True)
    best = None
    for it in range(steps):
        mid = math.sqrt(lo*hi)
        m, rg, nB, tt = certify(c, 'even', repr(mid), gball)
        good = bool(m > 0) and rg
        print(json.dumps({'eps': mid, 'margin': m.str(6, radius=False), 'rigorous': rg, 'ok': good, 't': round(tt, 1)}), flush=True)
        if good: lo = mid; best = mid
        else: hi = mid
    print(json.dumps({'certified_eps_deleted': best, 'upper_fail': hi}), flush=True)
