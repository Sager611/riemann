# cert.py (weildefl): deflated mixed-precision certificate  Q(f) >= eps ||f||^2  for all f in L^2(-L/2, L/2) of one parity, L = ln c.
# Chain (DERIVATION.md):
#  (A) lattice minorant (Chekmachev, Zenodo 23848597, Thm 5.10 + Def. 6.1; ../WeilZeta/DERIVATION.md sec. 5, rigorous scan RIGOR.md):
#      Q(f) - eps ||f||^2 >= (P - N)(f), P = P0 + U D U^T > 0, N = B B^T, and P - N >= 0  <=>  T := I - X + Y^T S^{-1} Y >= 0
#      (X = B^T P0^{-1} B, Y = U^T P0^{-1} B, S = D^{-1} + U^T P0^{-1} U; closed forms, assembled here exactly as in ../WeilZeta/cert_zeta.py).
#  (B) deflation lemma (DERIVATION.md sec. 3): for ANY exact m x k matrix Yt and delta > 0, T >= 0 if
#      (i)  z^T T z >= delta |z|^2 for every z orthogonal to Yt, and
#      (ii) Yt^T T Yt - (P_perp T Yt)^T (P_perp T Yt)/delta >= 0   (k x k), P_perp = I - Yt (Yt^T Yt)^{-1} Yt^T.
#  (i) is proved at LOW precision p_lo: T >= T_lo(Ws) = I - X + Y^T Ws + Ws^T Y - Ws^T S Ws for any Ws (Chekmachev Lemma 5.3), and
#      R_aug = T_lo + kappa Yt Yt^T - delta I >= 0 by a floating Cholesky + Arb residual (Chekmachev Lemma 5.5); z ⊥ Yt gives (i).
#  (ii) needs HIGH precision p_hi but only k columns: Yt^T T Yt >= Lk = Yt^T Yt - Yt^T X Yt + G^T W + W^T G - W^T S W (G = Y Yt, any W),
#      |P_perp T Yt a| <= |P_perp E1 a| + rho |a|, E1 = -X Yt + Y^T W, rho = ||Y||_F max(D) ||G - S W||_F  (S >= D^{-1}).
# Yt = the near-null invariant subspace of T, from a float64 eigensolve of T, refined at p_hi (correction equation solved in float64).
# usage: cert.py c sector eps p_hi [p_lo] [K2]     env: WZ_RSCAN (rigorous scan binary), WZ_SCAN (scan cache dir), WD_KTAU (deflation cut)
import sys, os, time, math, json, resource
import numpy as np
from flint import arb, acb, arb_mat, acb_mat, ctx, fmpq
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'WeilZeta'))
import cert_zeta as cz     # SymZeta, scan_end, rscan, fscan, parse_c (unchanged; audited in ../WeilZeta, ../ChekAudit)

class SymDH(cz.SymZeta):
    """ADVERSARIAL TEST ONLY: Davenport-Heilbronn function (conductor 5, odd gamma factor Gamma((s+1)/2), real coefficients a_n, no Euler
    product, zeros off the line). Lattice symbol of the same shape: mu(k) = log(5/pi) + Re psi(3/4 + i r_k/2) + 2(-1)^k sum c^{-b_n} b_n/(b_n^2 + r^2)
    - 2 sum_{m<c} (b_m/sqrt m) cos(r_k ln m), b_n = 2n + 3/2, b_m = coefficients of -f'/f (signed, every m); no pole terms."""
    def __init__(s, cq):
        s.key, s.q, s.a = 'dh', 5, 1
        s.cq = cq; s.c = arb(cq); s.L = s.c.log(); s.pi = arb.pi()
        s.const = (arb(5)/s.pi).log()
        r5 = arb(5).sqrt(); kap = ((10 - 2*r5).sqrt() - 2)/(r5 - 1)
        nmax = int(math.ceil(int(cq.p)/int(cq.q))) + 1
        av = [arb(0)]*(nmax + 1)
        for n in range(1, nmax + 1): av[n] = {1: arb(1), 2: kap, 3: -kap, 4: arb(-1), 0: arb(0)}[n % 5]
        bv = [arb(0)]*(nmax + 1)
        for n in range(2, nmax + 1):
            t = av[n]*arb(n).log()
            for d in range(2, n + 1):
                if n % d == 0: t -= av[d]*bv[n//d]
            bv[n] = t
        s.atoms = []
        n = 2
        while fmpq(n) < cq:
            if not (bv[n] == 0): s.atoms.append((n, arb(n).log(), bv[n]/arb(n).sqrt()))
            n += 1
        s.sig2 = 2*sum((abs(w) for _, _, w in s.atoms), arb(0))
        s.tail = []; nn = 0
        while True:
            b = arb(2*nn) + arb(3)/2; t = s.c**(-b)
            s.tail.append((b, t)); nn += 1
            if t < arb(2)**(-(ctx.prec + 40)) and nn > 2: break
        bN = arb(2*nn) + arb(3)/2
        s.trem = s.c**(-bN)/(bN*(1 - s.c**(-2)))
    def h(s, r): return s.const + acb(arb(3)/4, r/2).digamma().real
FORM = os.environ.get('WD_FORM', 'zeta')
if FORM == 'dh': cz.SymZeta = SymDH      # cz.rscan and the 96-bit re-check build their symbol through cz.SymZeta

def rss_mb(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20   # macOS: bytes

class MatGen:
    """Matrix-free access to X, Y, S (lean mode, WD_LEAN=1): entries from the closed forms on demand, products in row blocks of bs rows,
    so only O(n) high-precision data (Cauchy sums F, F2, column data) are kept. Same numbers as the stored matrices."""
    def __init__(s, gram, nB, nU, Md, bs=48):
        s.gram, s.nB, s.nU, s.Md, s.bs = gram, nB, nU, Md, bs
    def entry(s, kind, i, j):
        if kind == 'X': return s.gram(i, j)
        if kind == 'Y': return s.gram(s.nB + i, j)
        v = s.gram(s.nB + i, s.nB + j)
        return v + 1/s.Md[i] if i == j else v
    def block(s, kind, i0, i1):
        ncol = s.nB if kind in ('X', 'Y') else s.nU
        M = arb_mat(i1 - i0, ncol)
        for i in range(i0, i1):
            for j in range(ncol): M[i - i0, j] = s.entry(kind, i, j)
        return M
    def mul(s, kind, V):
        kc = V.ncols()
        if kind in ('X', 'S'):                 # symmetric: each entry evaluated once per product (upper triangle)
            n = s.nB if kind == 'X' else s.nU
            out = arb_mat(n, kc)
            for i0 in range(0, n, s.bs):
                i1 = min(n, i0 + s.bs); b = i1 - i0
                Bd = arb_mat(b, b)
                for i in range(b):
                    for j in range(i, b):
                        x = s.entry(kind, i0 + i, i0 + j); Bd[i, j] = x; Bd[j, i] = x
                Vd = arb_mat(b, kc)
                for i in range(b):
                    for j in range(kc): Vd[i, j] = V[i0 + i, j]
                P = Bd*Vd
                if i1 < n:
                    Bo = arb_mat(b, n - i1)
                    for i in range(b):
                        for j in range(i1, n): Bo[i, j - i1] = s.entry(kind, i0 + i, j)
                    Vo = arb_mat(n - i1, kc)
                    for i in range(i1, n):
                        for j in range(kc): Vo[i - i1, j] = V[i, j]
                    P = P + Bo*Vo
                    Q2 = Bo.transpose()*Vd
                    for i in range(n - i1):
                        for j in range(kc): out[i1 + i, j] += Q2[i, j]
                for i in range(b):
                    for j in range(kc): out[i0 + i, j] += P[i, j]
            return out
        if kind == 'YT':                       # Y^T V, V: nU x kc
            out = arb_mat(s.nB, kc)
            for i0 in range(0, s.nU, s.bs):
                i1 = min(s.nU, i0 + s.bs); Vb = arb_mat(i1 - i0, kc)
                for i in range(i0, i1):
                    for j in range(kc): Vb[i - i0, j] = V[i, j]
                out = out + s.block('Y', i0, i1).transpose()*Vb
            return out
        nrow = s.nB if kind == 'X' else s.nU
        out = arb_mat(nrow, kc)
        for i0 in range(0, nrow, s.bs):
            i1 = min(nrow, i0 + s.bs); P = s.block(kind, i0, i1)*V
            for i in range(i1 - i0):
                for j in range(kc): out[i0 + i, j] = P[i, j]
        return out
    def doubles(s):
        """float64 copies and error bounds of X, Y, S (one pass, symmetry used for X and S) and an upper bound of ||Y||_F."""
        nB, nU = s.nB, s.nU
        Xd = np.empty((nB, nB)); Yd = np.empty((nU, nB)); Sd = np.empty((nU, nU)); e = [arb(0)]*3; fy = arb(0)
        def put(D, i, j, x, t):
            d = float(x.mid()); D[i, j] = d; r = abs(x - arb(d)).upper()
            if r > e[t]: e[t] = r
            return d
        for i in range(nB):
            for j in range(i, nB):
                d = put(Xd, i, j, s.entry('X', i, j), 0); Xd[j, i] = d
        for i in range(nU):
            for j in range(nB):
                x = s.entry('Y', i, j); put(Yd, i, j, x, 1); fy += x*x
            for j in range(i, nU):
                d = put(Sd, i, j, s.entry('S', i, j), 2); Sd[j, i] = d
        ee = [math.nextafter(float(arb(t).upper()), math.inf) for t in e]
        return (Xd, ee[0]), (Yd, ee[1]), (Sd, ee[2]), arb(fy.upper()).sqrt()

def a2s(x, nd): return x.str(nd, radius=True, more=True)
def build(cq, sector, eps, K2arg=None, perturb=None, lean=False):
    """Node sets, weights and the Woodbury matrices X, Y, S at ctx.prec, exactly as ../WeilZeta/cert_zeta.py main().
    perturb = (n, factor): multiply the atom weight w_n by factor (adversarial test only)."""
    t0 = time.time()
    S = cz.SymZeta(cq); L = S.L; pi = S.pi
    if perturb is not None:
        n0, fac = perturb
        S.atoms = [(n, ln, (w*(1 + arb(fac)) if n == n0 else w)) for (n, ln, w) in S.atoms]
        S.sig2 = 2*sum((abs(w) for _, _, w in S.atoms), arb(0))
    cf = int(cq.p)/int(cq.q)
    K2 = max(201, math.ceil(5*cf*cf/3)); K2 += (K2 % 2 == 0)
    if K2arg is not None: K2 = K2arg + (K2arg % 2 == 0)
    J = K2
    T1s = os.environ.get('WZ_T1', '1'); t1 = arb(T1s); t1f = float(T1s)
    ks = cz.scan_end(S, K2, t1f); kse = cz.scan_end(S, 2*J, 0)
    rexe = os.environ.get('WZ_RSCAN')
    tag = ('' if perturb is None else f"_p{perturb[0]}x{perturb[1]}") + ('' if FORM == 'zeta' else '_' + FORM)
    cache = os.path.join(os.environ.get('WZ_SCAN', '.'), f"scan_zeta_{str(cq).replace('/', 'o')}_{K2}{tag}" + ("_rig" if rexe else "_flt") + ".json")
    if os.path.exists(cache):
        sc = json.load(open(cache)); assert sc['ks'] == ks and sc['kse'] == kse
        selo, sele, sinfo = sc['selo'], sc['sele'], sc.get('scan', {})
    else:
        sinfo = dict(method='float')
        if rexe: fo, fe, sinfo = cz.rscan(S, K2, J, ks, kse, rexe, t1f)
        else: fo, fe = cz.fscan(S, K2, J, ks, kse, t1f)
        oldp = ctx.prec; ctx.prec = 96
        S96 = cz.SymZeta(cq)
        if perturb is not None: S96.atoms = [(n, ln, (w*(1 + arb(perturb[1])) if n == perturb[0] else w)) for (n, ln, w) in S96.atoms]
        selo = [k for k in fo if S96.mu(k).lower() < t1f]
        sele = [j for j in fe if S96.mu(2*j).lower() < 0]
        ctx.prec = oldp
        json.dump(dict(ks=ks, kse=kse, selo=selo, sele=sele, flagged=[len(fo), len(fe)], scan=sinfo), open(cache, 'w'))
    tscan = time.time() - t0
    NBO, NBEK, NBEJ = 32, 1, 3
    Kset = set(range(1, K2+1, 2)) | set(selo); Eset = set(range(0, J+1)) | set(sele)
    for k in selo:
        for j in range((k-1)//2 - NBO, (k-1)//2 + NBO + 1):
            if j >= 0: Eset.add(j)
    for j in sele:
        for jj in range(j-NBEJ, j+NBEJ+1):
            if jj >= 0: Eset.add(jj)
    mus = {}
    def mu(k):
        if k not in mus: mus[k] = S.mu(k)
        return mus[k]
    changed = True; done = set()
    while changed:
        changed = False
        for j in sorted(Eset):
            if j in done or 2*j <= K2: continue
            done.add(j)
            if mu(2*j).mid() < 0:
                for k in range(2*j-2*NBEK+1, 2*j+2*NBEK, 2):
                    if k > 0 and k not in Kset: Kset.add(k); changed = True
                for jj in range(j-NBEJ, j+NBEJ+1):
                    if jj >= 0 and jj not in Eset: Eset.add(jj); changed = True
    if sector == 'odd': Eset.discard(0)
    Kl = sorted(Kset); El = sorted(Eset)
    w = {k: ((mu(k) - t1)/L).lower() for k in Kl}
    deff = t1/2 - eps
    assert deff > 0
    dl = {j: ((mu(2*j) + t1)/2 - eps).lower() for j in El}
    Kp = [k for k in Kl if w[k] > 0]; Km = [k for k in Kl if not (w[k] > 0)]
    G = [j for j in El if dl[j] > 0]; Jm = [j for j in El if not (dl[j] > 0)]
    tsym = time.time() - t0
    def Nj(j): return (1/L).sqrt() if (j == 0 and sector == 'even') else (2/L).sqrt()
    sgn = lambda j: 1 if j % 2 == 0 else -1
    sgk = lambda k: 1 if ((k-1)//2) % 2 == 0 else -1
    sL4, cL4, sL2 = (L/4).sinh(), (L/4).cosh(), (L/2).sinh()
    yv = -(L*L)/(4*pi*pi)
    if sector == 'even':
        A = lambda j: Nj(j)*(L/pi)*sgn(j); bk = lambda k: arb(2*k*sgk(k)); bv = -(L/pi)*sL4
        ip_node_v = lambda k: 2*S.r(k)*sgk(k)*cL4/(S.r(k)**2 + arb(1)/4); vv = sL2 + L/2
    else:
        A = lambda j: Nj(j)*sgn(j)*4*j*(L/pi); bk = lambda k: arb(sgk(k)); bv = sL4
        ip_node_v = lambda k: sgk(k)*cL4/(S.r(k)**2 + arb(1)/4); vv = sL2 - L/2
    cols = [('phi', j, (deff - dl[j]).sqrt(), None, None) for j in Jm]
    cols += [('node', k, (-w[k]).sqrt(), arb(k*k), bk(k)) for k in Km]
    if sector == 'odd' and FORM == 'zeta': cols.append(('v', None, arb(2).sqrt(), yv, bv))
    nB = len(cols)
    ucols = [('node', k, arb(1), arb(k*k), bk(k)) for k in Kp]; Md = [w[k] for k in Kp]
    if sector == 'even' and FORM == 'zeta': ucols.append(('v', None, arb(1), yv, bv)); Md.append(arb(2))
    allc = cols + ucols
    ckpt = os.environ.get('WD_CKPT'); nd = int(ctx.prec*0.302) + 15
    if ckpt and os.path.exists(ckpt):          # Cauchy sums from an earlier run of the same (c, sector, eps, prec): balls re-read from decimal
        ck = json.load(open(ckpt))             # strings with radius, which enclose the saved balls
        assert ck['key'] == [str(cq), sector, str(eps), ctx.prec, len(allc), len(G)], 'checkpoint mismatch'
        F = {int(i): arb(v) for i, v in ck['F'].items()}; F2 = {int(i): arb(v) for i, v in ck['F2'].items()}
    else:
        GA2 = [((1/dl[j] - 1/deff)*A(j)**2, arb(4*j*j)) for j in G]
        F = {}; F2 = {}
        for ci, col in enumerate(allc):
            if col[0] == 'phi': continue
            y = col[3]; f1 = arb(0); f2 = arb(0)
            for (ca, x) in GA2:
                d = 1/(y - x); f1 += ca*d; f2 += ca*d*d
            F[ci] = f1; F2[ci] = f2
        if ckpt:
            json.dump(dict(key=[str(cq), sector, str(eps), ctx.prec, len(allc), len(G)], F={str(i): a2s(v, nd) for i, v in F.items()},
                           F2={str(i): a2s(v, nd) for i, v in F2.items()}), open(ckpt, 'w'))
    tF = time.time() - t0
    def ipf(a, b):
        ta, ia = a[0], a[1]; tb, ib = b[0], b[1]
        if ta == 'phi' and tb == 'phi': return arb(1) if ia == ib else arb(0)
        if ta == 'phi': return A(ia)*b[4]/(b[3] - arb(4*ia*ia))
        if tb == 'phi': return A(ib)*a[4]/(a[3] - arb(4*ib*ib))
        if ta == 'node' and tb == 'node': return L/2 if ia == ib else arb(0)
        if ta == 'node' and tb == 'v': return ip_node_v(ia)
        if ta == 'v' and tb == 'node': return ip_node_v(ib)
        return vv
    def gram(p, qq):
        a = allc[p]; b = allc[qq]
        v = ipf(a, b)/deff
        if a[0] != 'phi' and b[0] != 'phi':
            if p == qq: v += a[4]*a[4]*F2[p]
            else: v += a[4]*b[4]*(F[p] - F[qq])/(b[3] - a[3])
        return v*a[2]*b[2]
    nU = len(ucols)
    if lean:
        maxMd = max((m_.upper() for m_ in Md), default=arb(0))
        meta = dict(K2=K2, kscan=ks, kscan_e=kse, sel_odd=len(selo), sel_even=len(sele), scan=sinfo.get('method', '?'), scanE=sinfo.get('E'),
                    Kp=len(Kp), Km=len(Km), Jm=len(Jm), G=len(G), order=nB, nU=nU, atoms=len(S.atoms),
                    t_scan=round(tscan, 1), t_sym=round(tsym - tscan, 1), t_cauchy=round(tF - tsym, 1), t_asm=0.0,
                    _cols=[(col[0], (-1 if col[1] is None else int(col[1])), float(col[2].mid())) for col in cols])
        return MatGen(gram, nB, nU, Md), None, None, arb(maxMd), meta
    X = arb_mat(nB, nB); Y = arb_mat(nU, nB); Sm = arb_mat(nU, nU)
    for p in range(nB):
        for qq in range(p, nB):
            v = gram(p, qq); X[p, qq] = v; X[qq, p] = v
    for a in range(nU):
        for b in range(nB): Y[a, b] = gram(nB + a, b)
        for b in range(a, nU):
            v = gram(nB + a, nB + b)
            if a == b: v = v + 1/Md[a]
            Sm[a, b] = v; Sm[b, a] = v
    tasm = time.time() - t0
    maxMd = max((m_.upper() for m_ in Md), default=arb(0))
    colinfo = [(col[0], (-1 if col[1] is None else int(col[1])), float(col[2].mid())) for col in cols]
    meta = dict(K2=K2, kscan=ks, kscan_e=kse, sel_odd=len(selo), sel_even=len(sele), scan=sinfo.get('method', '?'), scanE=sinfo.get('E'),
                Kp=len(Kp), Km=len(Km), Jm=len(Jm), G=len(G), order=nB, nU=nU, atoms=len(S.atoms),
                t_scan=round(tscan, 1), t_sym=round(tsym - tscan, 1), t_cauchy=round(tF - tsym, 1), t_asm=round(tasm - tF, 1))
    meta['_cols'] = colinfo
    return X, Y, Sm, arb(maxMd), meta

def tonp(M):
    return np.array([[float(M[i, j].mid()) for j in range(M.ncols())] for i in range(M.nrows())])
def tonp_err(M):
    """float64 midpoints D and a rigorous bound e >= max_ij |M_ij - D_ij| (ball radius plus rounding to binary64)."""
    r, c = M.nrows(), M.ncols(); D = np.empty((r, c)); e = arb(0)
    for i in range(r):
        for j in range(c):
            x = M[i, j]; d = float(x.mid()); D[i, j] = d
            t = abs(x - arb(d)).upper()
            if t > e: e = t
    return D, math.nextafter(float(arb(e).upper()), math.inf)
U53 = 2.0**-53
def gam(n): return n*U53/(1 - n*U53)
def fro(A): return math.nextafter(float(np.linalg.norm(A))*(1 + 1e-12), math.inf)   # Frobenius norm, rounded up (rel. error << 1e-12)
def compl_double(Xd, eX, Yd, eY, Sd, eS, Wsd, Ytd, eYt, dlt, kappa=1.0):
    """(i) in binary64 with a rigorous error bound (DERIVATION.md sec. 4b): returns (ok, lam, s0, nu, Etot, ||E||_F)."""
    m, nU = Xd.shape[0], Sd.shape[0]; k = Ytd.shape[1]
    A1 = Yd.T @ Wsd; B1 = Sd @ Wsd; A2 = Wsd.T @ B1; C = Ytd @ Ytd.T
    R = np.eye(m) - Xd + A1 + A1.T - A2 + kappa*C - dlt*np.eye(m)
    Rs = (R + R.T)/2
    fW, fY, fS, fX, fYt = fro(Wsd), fro(Yd), fro(Sd), fro(Xd), fro(Ytd)
    E_in = m*eX + 2*math.sqrt(nU*m)*eY*fW + fW*fW*nU*eS + abs(kappa)*2*(fYt + math.sqrt(m*k)*eYt)*math.sqrt(m*k)*eYt
    E_mm = 2*gam(nU)*fY*fW + gam(nU)*fS*fW*fW + gam(nU)*fW*fro(B1) + abs(kappa)*gam(k)*fYt*fYt
    E_sum = 8*U53*(math.sqrt(m) + fX + 2*fro(A1) + fro(A2) + abs(kappa)*fro(C) + dlt*math.sqrt(m)) + U53*fro(R)
    Etot = 1.01*(E_in + E_mm + E_sum)
    lam = float(np.linalg.eigvalsh(Rs)[0])
    if not lam > 0: return False, lam, None, None, Etot, None
    for fac in (0.5, 0.25, 0.1):
        s0 = lam*fac; nu = s0/2
        try: Lc = np.linalg.cholesky(Rs - s0*np.eye(m))
        except np.linalg.LinAlgError: continue
        LL = Lc @ Lc.T; Ed = Rs - s0*np.eye(m) - LL
        errE = gam(m)*fro(Lc)**2 + 3*U53*(fro(Rs) + s0*math.sqrt(m) + fro(LL))
        bound = 1.01*(fro(Ed) + errE) + Etot
        if bound <= s0 - nu: return True, lam, s0, nu, Etot, bound
    return False, lam, None, None, Etot, None
def toarb(A):
    A = np.asarray(A, dtype=float); r, c = A.shape; M = arb_mat(r, c)
    for i in range(r):
        for j in range(c): M[i, j] = arb(float(A[i, j]))
    return M
def frob(M):   # rigorous UPPER bound of the Frobenius norm
    s = arb(0)
    for i in range(M.nrows()):
        for j in range(M.ncols()): s += M[i, j]*M[i, j]      # (x**2 is NaN in python-flint for balls containing 0)
    return arb(s.upper()).sqrt()
def main():
    cq = cz.parse_c(sys.argv[1]); sector = sys.argv[2]; epsstr = sys.argv[3]; p_hi = int(sys.argv[4])
    p_lo = int(sys.argv[5]) if len(sys.argv) > 5 and sys.argv[5] != '-' else 96
    K2arg = int(sys.argv[6]) if len(sys.argv) > 6 and sys.argv[6] != '-' else None
    pert = os.environ.get('WD_PERTURB')          # adversarial: "n:factor"
    perturb = None if not pert else (int(pert.split(':')[0]), pert.split(':')[1])   # "n:d" -> w_n *= (1 + d)
    ctx.prec = p_hi
    eps = arb(0) if epsstr == '0' else (arb(2)**(-int(epsstr[3:])) if epsstr.startswith('2^-') else arb(epsstr))
    t0 = time.time()
    lean = os.environ.get('WD_LEAN') == '1'
    X, Y, Sm, maxMd, meta = build(cq, sector, eps, K2arg, perturb, lean=lean)
    colinfo = meta.pop('_cols')
    tb = time.time()
    if lean:
        gen = X; m, nU = gen.nB, gen.nU
        (Xd, eX), (Yd, eY), (Sd, eS), frobY = gen.doubles()
        mulX = lambda V: gen.mul('X', V); mulY = lambda V: gen.mul('Y', V); mulYT = lambda V: gen.mul('YT', V); mulS = lambda V: gen.mul('S', V)
        os.environ['WD_COMPL'] = 'double'
    else:
        m, nU = X.nrows(), Sm.nrows()
        (Xd, eX), (Yd, eY), (Sd, eS) = tonp_err(X), tonp_err(Y), tonp_err(Sm)
        mulX = lambda V: X*V; mulY = lambda V: Y*V; mulYT = lambda V: Y.transpose()*V; mulS = lambda V: Sm*V
        frobY = None
    meta['t_doubles'] = round(time.time() - tb, 1)
    # ---- float64 picture of T
    Wsd = np.linalg.solve(Sd, Yd)
    Td = np.eye(m) - Xd + Yd.T @ Wsd; Td = (Td + Td.T)/2
    ev, V = np.linalg.eigh(Td)
    tau = float(os.environ.get('WD_KTAU', '1e-2'))
    k = int(np.sum(ev < tau))
    if k == 0: k = 1
    if k >= m: k = m - 1
    delta0 = float(ev[k]); dlt = arb(delta0/2)       # delta = half the first non-deflated float eigenvalue
    kappa = arb(1)
    Yt = V[:, :k].copy()
    Kaug = Td + V[:, :k] @ V[:, :k].T
    cK = np.linalg.cholesky(Kaug)
    def ksolve(B): return np.linalg.solve(cK.T, np.linalg.solve(cK, B))
    td = time.time()
    # ---- refinement of the invariant subspace (k columns only) at p_mid = p_hi/2 + 64 bits: it only has to reach
    #      |P_perp T Yt| ~ 2^-(p_hi/2 + 16), since the Schur term is quadratic in it (DERIVATION.md sec. 4)
    condS = float(np.linalg.cond(Sd))
    def to_mid(M): return M.mid()
    def solveS(Gm, W0, tol):
        Wm = W0 if W0 is not None else toarb(np.linalg.solve(Sd, tonp(Gm)))
        R = Gm - mulS(Wm); nit = 0
        while nit < 30:
            nr = float(frob(R).mid())
            if nr == 0 or nr < tol: break
            Wm = (Wm + toarb(np.linalg.solve(Sd, tonp(R)))).mid(); nit += 1
            R = Gm - mulS(Wm)
        return Wm, R, nit
    hist = []; nS = 0
    p_mid = min(p_hi, p_hi//2 + 64); ctx.prec = p_mid
    YtA = toarb(Yt); Wm = None
    target = 2.0**(-(p_hi//2 + 16))
    npr_prev = 1e-13
    for outer in range(16):
        G_ = mulY(YtA)
        Wm, Rk, nit = solveS(G_, Wm, max(target*1e-6, npr_prev*1e-17)); nS += nit   # S-solve only as accurate as this step needs
        TY = YtA - mulX(YtA) + mulYT(Wm)
        Gk = YtA.transpose()*YtA; Gki = Gk.inv()
        H = Gki*(YtA.transpose()*TY)
        Rr = TY - YtA*H
        PR = Rr - YtA*(Gki*(YtA.transpose()*Rr))
        npr = float(frob(PR).mid()); hist.append('%.2e' % npr); npr_prev = npr
        if npr < target or outer == 15: break
        # correction equation T_perp E - E H = PR, solved in float64 (E orthogonal to Yt), then Yt <- Yt - E
        PRd = tonp(PR); Hd = tonp(H); Ytd = tonp(YtA)
        E = ksolve(PRd)
        for _ in range(10): E = ksolve(PRd + E @ Hd)
        E = E - Ytd @ np.linalg.solve(Ytd.T @ Ytd, Ytd.T @ E)
        YtA = (YtA - toarb(E)).mid()
    ctx.prec = p_hi
    tr = time.time()
    # ---- (ii): rigorous k x k block at p_hi
    G_ = mulY(YtA)
    Wm, Rk, nit = solveS(G_, Wm, target*1e-6); nS += nit
    XY = mulX(YtA); SW = mulS(Wm)
    Gk = YtA.transpose()*YtA; Gki = Gk.inv()
    Lk = Gk - YtA.transpose()*XY + G_.transpose()*Wm + Wm.transpose()*G_ - Wm.transpose()*SW
    E1 = mulYT(Wm) - XY
    PE1 = E1 - YtA*(Gki*(YtA.transpose()*E1))
    rho = (frob(Y) if frobY is None else frobY)*maxMd*frob(Rk)
    Pk = 2*(PE1.transpose()*PE1)
    for i in range(k): Pk[i, i] += 2*rho*rho
    Sk = Lk - Pk*(1/dlt)
    # positivity of the k x k block (Chekmachev Lemma 5.5 style, robust for graded blocks): s0 = margin/2, an mpmath Cholesky L of
    # mid(Sk) - s0 I, and the Arb residual Ek = Sk - L L^T - s0 I; if max row sum |Ek| <= s0 - nu2 then lambda_min(Sk) >= nu2 > 0.
    import mpmath as mp
    dps = int(p_hi*0.3) + 10; mp.mp.dps = dps
    def tomp(M):
        A = mp.matrix([[mp.mpf(M[i, j].mid().str(dps, radius=False)) for j in range(M.ncols())] for i in range(M.nrows())])
        return (A + A.T)/2
    def mineig(A):
        x = mp.matrix([1]*A.rows)
        for _ in range(60):
            x = mp.lu_solve(A, x); x = x/mp.norm(x)
        return (x.T*A*x)[0]/(x.T*x)[0]
    Skm = tomp(Sk); sk_min = mineig(Skm); lk_min = mineig(tomp(Lk))
    ok2 = False; nu2 = None; ek = None
    if sk_min > 0:
        for fac in ('0.5', '0.25', '0.1'):
            s0k = sk_min*mp.mpf(fac)
            try: Lmp = mp.cholesky(Skm - s0k*mp.eye(k))
            except Exception: continue
            LA = arb_mat([[arb(mp.nstr(Lmp[i, j], dps, strip_zeros=False)) if j <= i else arb(0) for j in range(k)] for i in range(k)])
            s0a = arb(mp.nstr(s0k, dps)); nu2a = s0a/2
            Ek = Sk - LA*LA.transpose()
            ek = max(sum((abs(Ek[i, j] - (s0a if i == j else 0)) for j in range(k)), arb(0)).upper() for i in range(k))
            if ek <= s0a - nu2a: ok2 = True; nu2 = nu2a; break
    piv = [nu2] if nu2 is not None else []
    rho_s = rho.str(3, radius=False); pe1_s = frob(PE1).str(3, radius=False)
    tk = time.time()
    # ---- (i): complement coercivity. Default: binary64 with a rigorous error bound (assumes IEEE binary64, round to nearest, classical
    #      dgemm: |fl(AB) - AB| <= gamma_n |A||B|). WD_COMPL=arb: the same matrix in Arb at p_lo bits (no assumption).
    mode = os.environ.get('WD_COMPL', 'double')
    Ytd, eYt = tonp_err(YtA)
    if os.environ.get('WD_DUMP'):     # optional: deflated subspace and the negative columns, for ladder.py (no role in the proof)
        np.savez(os.environ['WD_DUMP'], Yt=Ytd, ctype=np.array([c_[0] for c_ in colinfo]), cidx=np.array([c_[1] for c_ in colinfo]),
                 cscale=np.array([c_[2] for c_ in colinfo]), ev=ev, L=float(arb(cq).log().mid()))
    if mode == 'double':
        del X, Y, Sm, XY, SW, E1, PE1, G_, Wm, Rk
        mulX = mulY = mulYT = mulS = None
        import gc; gc.collect()
        ok1, lam, s0, nu, Etot, emax_f = compl_double(Xd, eX, Yd, eY, Sd, eS, Wsd, Ytd, eYt, float(dlt.mid()))
        emax = None if emax_f is None else arb(emax_f)
        compl_info = dict(Etot='%.2e' % Etot, eX='%.1e' % eX, eY='%.1e' % eY, eS='%.1e' % eS, eYt='%.1e' % eYt, normW='%.2e' % fro(Wsd))
    else:
        ctx.prec = p_lo
        Ws = toarb(Wsd)
        Tlo = arb_mat(m, m)
        P1 = Y.transpose()*Ws; P2 = Ws.transpose()*(Sm*Ws)
        YY = YtA*YtA.transpose()
        for i in range(m):
            for j in range(m):
                Tlo[i, j] = (1 if i == j else 0) - X[i, j] + P1[i, j] + P1[j, i] - P2[i, j] + kappa*YY[i, j] - (dlt if i == j else 0)
        del P1, P2, YY
        Ra = (tonp(Tlo) + tonp(Tlo).T)/2
        lam = float(np.linalg.eigvalsh(Ra)[0])
        ok1 = False; s0 = nu = emax = None
        if lam > 0:
            for fac in (0.5, 0.25, 0.1):
                s0 = lam*fac; nu = s0/2
                try: Lc = np.linalg.cholesky(Ra - s0*np.eye(m))
                except np.linalg.LinAlgError: continue
                LA = toarb(Lc); Rres = Tlo - LA*LA.transpose()
                emax = max(sum((abs(Rres[i, j] - (s0 if i == j else 0)) for j in range(m)), arb(0)).upper() for i in range(m))
                if emax <= arb(s0 - nu): ok1 = True; break
        compl_info = {}
    t1 = time.time()
    res = dict(form=FORM, c=str(cq), sector=sector, eps=epsstr, perturb=pert, p_hi=p_hi, p_lo=p_lo, **meta, k=k, delta0=delta0, delta=float(dlt.mid()),
               float_ev=['%.2e' % e for e in ev[:k + 3]], refine=hist, p_mid=p_mid, condS='%.1e' % condS, S_iters=nS,
               ii_ok=bool(ok2), ii_margin=mp.nstr(sk_min, 12), Lk_min=mp.nstr(lk_min, 12), ii_nu=(nu2.str(6, radius=False) if nu2 is not None else None), ii_resid=(ek.str(3, radius=False) if ek is not None else None),
               rho=rho_s, PE1=pe1_s,
               i_mode=mode, i_ok=bool(ok1), i_lam=lam, i_nu=nu, i_resid=(emax.str(3, radius=False) if emax is not None else None), **compl_info,
               rigorous=bool(ok1 and ok2),
               lean=lean, t_float=round(td - tb, 1), t_refine=round(tr - td, 1), t_block_hi=round(tk - tr, 1), t_compl_lo=round(t1 - tk, 1), t_total=round(t1 - t0, 1), rss_mb=round(rss_mb(), 1))
    print(json.dumps(res), flush=True)

if __name__ == '__main__':
    main()
