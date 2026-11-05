# Scalable version of cert.py: float64 scan (uint64 phases, margin 1e-7) + Arb balls at flagged/retained nodes,
# Gram sums over G by the Cauchy formulas (Prop. 5.2, re-derived), exact Woodbury W = X - Y^T S^{-1} Y in Arb,
# positivity of Id - W by Lemma 5.5 (approx. Cholesky + Arb residual).  Usage: cert2.py c sector eps bits
import sys, time, math, random
import numpy as np, mpmath as mp
from flint import arb, acb, arb_mat, ctx
from cert import Sym, scan_end

def fscan(c, L_mp, K2, J, ks, kse):
    mp.mp.dps = 40
    pp = []
    for n in range(2, c):
        f = [p for p in range(2, n+1) if n % p == 0 and all(p % q for q in range(2, p))]
        if len(f) == 1: pp.append((n, mp.log(f[0])/mp.sqrt(n)))
    B = np.array([int(mp.nint(mp.log(n)/(2*L_mp)*mp.mpf(2)**64)) for n, _ in pp], dtype=np.uint64)
    W = np.array([float(w) for _, w in pp]); Lf = float(L_mp); lnpi = float(mp.log(mp.pi))
    an = np.array([2*n + 0.5 for n in range(6)]); can = np.array([float(mp.mpf(c)**(-(2*n + mp.mpf(1)/2))) for n in range(6)])
    def mu(k):
        k = k.astype(np.uint64); r = np.pi*k.astype(np.float64)/Lf; w = 0.25 + 0.5j*r
        psi = np.log(w) - 1/(2*w) - 1/(12*w**2) + 1/(120*w**4) - 1/(252*w**6) + 1/(240*w**8)
        sgn = np.where(k % 2 == 0, 1.0, -1.0)
        T = -2*sgn*sum(can[i]*an[i]/(an[i]**2 + r*r) for i in range(6))
        s = np.zeros_like(r)
        for i in range(len(B)): s += W[i]*np.cos(2*np.pi*((k*B[i]).astype(np.float64)/2.0**64))
        return psi.real - lnpi - T - 2*s
    selo, sele = [], []; CH = 2_000_000
    k0 = K2 + 2
    while k0 <= ks:
        k = np.arange(k0, min(k0 + 2*CH, ks + 1), 2, dtype=np.int64); m = mu(k); selo.extend(k[m < 1 + 1e-7].tolist()); k0 += 2*CH
    j0 = J + 1
    while 2*j0 <= kse:
        j = np.arange(j0, min(j0 + CH, kse//2 + 1), dtype=np.int64); m = mu(2*j); sele.extend(j[m < 1e-7].tolist()); j0 += CH
    return selo, sele

def main():
    c = int(sys.argv[1]); sector = sys.argv[2]; epsstr = sys.argv[3]; ctx.prec = int(sys.argv[4])
    eps = arb(epsstr) if epsstr != '0' else arb(0)
    t0 = time.time()
    S = Sym(c); L = S.L; pi = S.pi
    K2 = max(201, math.ceil(5*c*c/3)); K2 += (K2 % 2 == 0); J = K2
    ks = scan_end(S, K2, 1); kse = scan_end(S, 2*J, 0)
    mp.mp.dps = 50
    fo, fe = fscan(c, mp.log(c), K2, J, ks, kse)
    oldp = ctx.prec; ctx.prec = 96; S96 = Sym(c)
    selo = [k for k in fo if S96.mu(k).lower() < 1]
    sele = [j for j in fe if S96.mu(2*j).lower() < 0]
    ctx.prec = oldp
    tscan = time.time() - t0
    Kset = set(range(1, K2+1, 2)) | set(selo); Eset = set(range(0, J+1)) | set(sele)
    for k in selo:
        for j in range((k-1)//2 - 32, (k-1)//2 + 33):
            if j >= 0: Eset.add(j)
    for j in sele:
        for jj in range(j-3, j+4): Eset.add(jj)
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
                for k in (2*j-1, 2*j+1):
                    if k not in Kset: Kset.add(k); changed = True
                for jj in range(j-3, j+4):
                    if jj not in Eset: Eset.add(jj); changed = True
    if sector == 'odd': Eset.discard(0)
    Kl = sorted(Kset); El = sorted(Eset)
    w = {k: ((mu(k) - 1)/L).lower() for k in Kl}
    deff = arb(1)/2 - eps
    dl = {j: ((mu(2*j) + 1)/2 - eps).lower() for j in El}
    Kp = [k for k in Kl if w[k] > 0]; Km = [k for k in Kl if not (w[k] > 0)]
    G = [j for j in El if dl[j] > 0]; Jm = [j for j in El if not (dl[j] > 0)]
    tsym = time.time() - t0
    sL4, cL4, sL2 = (L/4).sinh(), (L/4).cosh(), (L/2).sinh()
    def Nj(j): return (1/L).sqrt() if (j == 0 and sector == 'even') else (2/L).sqrt()
    sgn = lambda j: 1 if j % 2 == 0 else -1
    sgk = lambda k: 1 if ((k-1)//2) % 2 == 0 else -1
    yv = -(L*L)/(4*pi*pi)
    if sector == 'even':
        A = lambda j: Nj(j)*(L/pi)*sgn(j); bk = lambda k: arb(2*k*sgk(k)); bv = -(L/pi)*sL4
        ip_node_v = lambda k: 2*S.r(k)*sgk(k)*cL4/(S.r(k)**2 + arb(1)/4); vv = sL2 + L/2
    else:
        A = lambda j: Nj(j)*sgn(j)*4*j*(L/pi); bk = lambda k: arb(sgk(k)); bv = sL4
        ip_node_v = lambda k: sgk(k)*cL4/(S.r(k)**2 + arb(1)/4); vv = sL2 - L/2
    # column descriptors: (type, index, scale, y, b)
    cols = [('phi', j, (deff - dl[j]).sqrt(), None, None) for j in Jm]
    cols += [('node', k, (-w[k]).sqrt(), arb(k*k), bk(k)) for k in Km]
    if sector == 'odd': cols.append(('v', None, arb(2).sqrt(), yv, bv))
    nB = len(cols)
    ucols = [('node', k, arb(1), arb(k*k), bk(k)) for k in Kp]; Md = [w[k] for k in Kp]
    if sector == 'even': ucols.append(('v', None, arb(1), yv, bv)); Md.append(arb(2))
    allc = cols + ucols; nA = len(allc)
    # Cauchy sums F(y), F2(y) over G
    GA2 = [((1/dl[j] - 1/deff)*A(j)**2, arb(4*j*j)) for j in G]
    F = {}; F2 = {}
    for ci, col in enumerate(allc):
        if col[0] == 'phi': continue
        y = col[3]; f1 = arb(0); f2 = arb(0)
        for (ca, x) in GA2:
            d = 1/(y - x); f1 += ca*d; f2 += ca*d*d
        F[ci] = f1; F2[ci] = f2
    tF = time.time() - t0
    def ipf(a, b):   # <a,b> unscaled
        ta, ia = a[0], a[1]; tb, ib = b[0], b[1]
        if ta == 'phi' and tb == 'phi': return arb(1) if ia == ib else arb(0)
        if ta == 'phi' and tb in ('node', 'v'): return A(ia)*b[4]/(b[3] - arb(4*ia*ia))
        if tb == 'phi' and ta in ('node', 'v'): return A(ib)*a[4]/(a[3] - arb(4*ib*ib))
        if ta == 'node' and tb == 'node': return L/2 if ia == ib else arb(0)
        if ta == 'node' and tb == 'v': return ip_node_v(ia)
        if ta == 'v' and tb == 'node': return ip_node_v(ib)
        return vv
    Gm = arb_mat(nA, nA)
    for p in range(nA):
        a = allc[p]
        for q in range(p, nA):
            b = allc[q]
            v = ipf(a, b)/deff
            if a[0] != 'phi' and b[0] != 'phi':
                if p == q: v += a[4]*a[4]*F2[p]
                else: v += a[4]*b[4]*(F[p] - F[q])/(b[3] - a[3])
            v = v*a[2]*b[2]
            Gm[p, q] = v; Gm[q, p] = v
    nU = len(ucols)
    X = arb_mat(nB, nB); Y = arb_mat(nU, nB); Sm = arb_mat(nU, nU)
    for a in range(nB):
        for b in range(nB): X[a, b] = Gm[a, b]
    for a in range(nU):
        for b in range(nB): Y[a, b] = Gm[nB+a, b]
        for b in range(nU): Sm[a, b] = Gm[nB+a, nB+b] + ((1/Md[a]) if a == b else 0)
    tasm = time.time() - t0
    Z = Sm.solve(Y)
    Wm = X - Y.transpose()*Z
    IW = arb_mat(nB, nB)
    for a in range(nB):
        for b in range(nB): IW[a, b] = (1 if a == b else 0) - Wm[a, b]
    tW = time.time() - t0
    R = arb_mat(nB, nB)
    for a in range(nB):
        for b in range(nB): R[a, b] = IW[a, b].mid()
    random.seed(1)
    x = arb_mat(nB, 1, [arb(random.uniform(-1, 1)) for _ in range(nB)])
    for it in range(3):
        yv_ = R.solve(x, algorithm="approx")
        nrm = sum((yv_[i, 0]**2 for i in range(nB)), arb(0)).sqrt()
        x = arb_mat(nB, 1, [(yv_[i, 0]/nrm).mid() for i in range(nB)])
    Rx = R*x; lmin = sum((x[i, 0]*Rx[i, 0] for i in range(nB)), arb(0)).mid()
    # Lemma 5.5 with nu = 0.99*lmin
    nu = lmin*arb('0.99'); s0 = lmin*arb('0.995'); n = nB
    Lc = [[None]*n for _ in range(n)]; okc = True
    for i in range(n):
        for j in range(i+1):
            s = R[i, j] - (s0 if i == j else 0)
            for kk in range(j): s -= Lc[i][kk]*Lc[j][kk]
            if i == j:
                if not (s > 0): okc = False; break
                Lc[i][j] = s.sqrt().mid()
            else: Lc[i][j] = (s/Lc[j][j]).mid()
        if not okc: break
    rig = None
    if okc:
        Lm = arb_mat(n, n)
        for i in range(n):
            for j in range(i+1): Lm[i, j] = Lc[i][j]
        Fm = IW - Lm*Lm.transpose(); emax = arb(0)
        for i in range(n):
            rs = arb(0)
            for j in range(n): rs += abs(Fm[i, j] - (s0 if i == j else 0)).upper()
            if rs.upper() > emax.upper(): emax = rs.upper()
        rig = bool(emax <= s0 - nu)
    print(f"c={c} sector={sector} eps={epsstr} prec={ctx.prec} K2=J={K2} kscan={ks} kscan_e={kse} sel_odd={len(selo)} sel_even={len(sele)} |K+|={len(Kp)} |K-|={len(Km)} |J-|={len(Jm)} |G|={len(G)} order={nB}")
    print(f"   margin 1-lmax(W) ~ {lmin.str(10, radius=False)}; Lemma 5.5 with nu = 0.99*margin: {rig}")
    print(f"   times: scan {tscan:.0f}s sym {tsym:.0f}s F {tF:.0f}s asm {tasm:.0f}s W {tW:.0f}s total {time.time()-t0:.0f}s", flush=True)

if __name__ == '__main__':
    main()
