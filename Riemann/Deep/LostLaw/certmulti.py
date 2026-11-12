# certmulti.py (lostlaw): Chekmachev's certificate (Thm 5.10 + Def. 6.1; ChekAudit cert.py; MoriRate certlow.py) for model forms
#   Q_S = Q - sum_{k in D} 2<f, cos(g_k x)>^2 + sum_{a in A} 2<f, cos(a x)>^2      (even sector)
# Deleted zero = extra negative rank-one column ('eval', g_k, sqrt 2) in B; added point = extra positive column ('eval', a, 1) in U
# with M-weight 2. Q_S - eps >= D + U M U* - B'B'* >= 0 iff lmax(B'* P^-1 B') <= 1 (Woodbury); positivity by Lemma 5.5.
# Changes vs certlow (same inequalities): scan once per c (cached); dl_j = d0_j - eps with d0_j = ((mu(2j)+1)/2).lower();
# Gram(eps) = G0 + eps*G1 + E, G0 = P0 + 2 IP, G1 = P1 + 4 IP, P_i = Phi^T diag(chi_i) Phi accumulated in chunks over G
# (chi_0 = 1/d0 - 2, chi_1 = 1/d0^2 - 4), E a ball of radius eps^2 (kappa n_a n_b + 16 |IP_ab|) (exact second-order remainder
# bound, n_a^2 = sum_G Phi_ja^2, kappa = 1/(dmin^2 (dmin - eps)) + 16); eps by linear extrapolation of the margin, then a
# rigorous Lemma 5.5 check at the returned eps. Usage: certmulti.py c prec sets.json gal.log [names,...]
import sys, os, json, time, math, resource
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from flint import arb, acb, arb_mat, acb_mat, ctx
from certlow import Sym, scan_end, Zfun

def rss(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
def log(d): print(json.dumps(d), flush=True)

def zero_ball(gs, rexp):
    r = arb(10)**(-rexp); g = arb(gs)
    zl = Zfun(g - r); zr = Zfun(g + r)
    return g.union(g - r).union(g + r), bool((zl < 0 and zr > 0) or (zl > 0 and zr < 0))

def main():
    c = int(sys.argv[1]); ctx.prec = int(sys.argv[2])
    sets = json.loads(open(sys.argv[3]).read())
    gal = {}
    for l in open(sys.argv[4]):
        try: d = json.loads(l)
        except Exception: continue
        if 'set' in d: gal[d['set']] = float(d['lam1'])
        if 'E1' in d: gal['full'] = float(d['E1'])
    only = sys.argv[5].split(',') if len(sys.argv) > 5 else None
    sets = [{'name': 'full', 'del': []}] + sets
    if only: sets = [s for s in sets if s['name'] in only]
    zs = [l.strip() for l in open(os.path.join(HERE, 'zeros_dps179.txt')) if l.strip()]
    t0 = time.time()
    S = Sym(c); L = S.L; pi = S.pi
    K2 = int(float(os.environ.get('K2MUL', '1'))*max(201, math.ceil(5*c*c/3))); K2 += (K2 % 2 == 0); J = K2
    cache = os.path.join(HERE, 'logs', 'scan_c%d.json' % c) if K2 == max(201, math.ceil(5*c*c/3)) + (max(201, math.ceil(5*c*c/3)) % 2 == 0) else os.path.join(HERE, 'logs', 'scan_c%d_K%d.json' % (c, K2))
    if os.path.exists(cache):
        sc = json.load(open(cache)); selo, sele = sc['selo'], sc['sele']
    else:
        ks = scan_end(S, K2, 1); kse = scan_end(S, 2*J, 0)
        oldp = ctx.prec; ctx.prec = 80; S80 = Sym(c)
        selo = [k for k in range(K2+2, ks+1, 2) if S80.mu(k).lower() < 1]
        sele = [j for j in range(J+1, kse//2 + 1) if S80.mu(2*j).lower() < 0]
        ctx.prec = oldp
        json.dump({'selo': selo, 'sele': sele, 'ks': ks, 'kse': kse}, open(cache, 'w'))
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
    Kl = sorted(Kset); El = sorted(Eset)
    w = {k: ((mu(k) - 1)/L).lower() for k in Kl}
    d0 = {j: ((mu(2*j) + 1)/2).lower() for j in El}
    Kp = [k for k in Kl if w[k] > 0]; Km = [k for k in Kl if not (w[k] > 0)]
    G = [j for j in El if d0[j] > 0]; Jm = [j for j in El if not (d0[j] > 0)]
    dmin = min(d0[j] for j in G)
    log({'c': c, 'prec': ctx.prec, 'K2': K2, 'sel_odd': len(selo), 'sel_even': len(sele), 'Kp': len(Kp), 'Km': len(Km), 'G': len(G),
         'Jm': len(Jm), 'dmin': arb(dmin).str(5, radius=False), 't_scan': round(tscan, 1), 't_sym': round(time.time() - t0, 1)})
    # columns
    def Nj(j): return (1/L).sqrt() if j == 0 else (2/L).sqrt()
    sL4, cL4, sL2 = (L/4).sinh(), (L/4).cosh(), (L/2).sinh()
    def sgk(k): return 1 if ((k-1)//2) % 2 == 0 else -1
    h = L/2
    def sincd(d):   # sin(d h)/d, safe near 0
        if abs(d.mid()) < arb(10)**(-30):
            z = d*h; return h*(1 + arb(0, (z*z/6).upper()))
        return (d*h).sin()/d
    def ip_phi_node(j, k):
        sj = 1 if j % 2 == 0 else -1
        return Nj(j)*(L/pi)*sj*2*k*sgk(k)/(k*k - 4*j*j)
    def ip_phi_v(j):
        sj = 1 if j % 2 == 0 else -1; om = 2*pi*j/L
        return Nj(j)*sj*sL4/(om*om + arb(1)/4)
    def ip_node_v(k):
        r = S.r(k); return 2*r*sgk(k)*cL4/(r*r + arb(1)/4)
    vv = sL2 + L/2
    def ip_eval_phi(g, j):
        om = 2*pi*j/L; return Nj(j)*(sincd(om - g) + sincd(om + g))
    def ip_eval_node(g, k):
        r = S.r(k); return sincd(r - g) + sincd(r + g)
    def ip_eval_v(g):
        return 2*((h/2).sinh()*(g*h).cos()/2 + g*(h/2).cosh()*(g*h).sin())/(g*g + arb(1)/4)
    def ip_eval_eval(g1, g2):
        return sincd(g1 - g2) + sincd(g1 + g2)
    def ip(a, b):
        ta, ia, _ = a; tb, ib, _ = b
        if ta == 'eval' or tb == 'eval':
            if ta == 'eval' and tb == 'eval': return ip_eval_eval(ia, ib)
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
    Bbase = [('phi', j, (arb(1)/2 - d0[j]).sqrt()) for j in Jm] + [('node', k, (-w[k]).sqrt()) for k in Km]
    Ubase = [('node', k, arb(1)) for k in Kp] + [('v', None, arb(1))]
    Mbase = [w[k] for k in Kp] + [arb(2)]
    rexp = min(170, int(ctx.prec*0.30) - 10)
    dk = sorted(set(k for s in sets for k in s.get('del', [])))
    ak = sorted(set((a[0], a[1]) for s in sets for a in s.get('add', [])))
    zb = {}
    for k in dk:
        gb, ok = zero_ball(zs[k-1], rexp); zb[k] = gb
        if not ok: log({'warn': 'zero ball not certified', 'k': k})
    Bx = [('eval', zb[k], arb(2).sqrt()) for k in dk]
    Ux = [('eval', arb(zs[a[0]-1]) + arb(a[1]), arb(1)) for a in ak]
    allc = Bbase + Bx + Ubase + Ux
    nA = len(allc); nBb = len(Bbase); nUb = len(Ubase)
    iBx = {k: nBb + i for i, k in enumerate(dk)}
    iU0 = nBb + len(Bx); iUx = {a: iU0 + nUb + i for i, a in enumerate(ak)}
    # IP (scaled), G0 = P0 + 2 IP, G1 = P1 + 4 IP
    ta = time.time()
    IP = arb_mat(nA, nA)
    for a in range(nA):
        for b in range(a, nA):
            v = ip(allc[a], allc[b])*allc[a][2]*allc[b][2]; IP[a, b] = v; IP[b, a] = v
    # memory-lean assembly: P0 starts at 2 IP (full precision), P1 = 4 IP + Phi^T diag(chi1) Phi kept at 128 bits (balls);
    # chunk products added in place (no nA x nA sums kept twice)
    ipn = [abs(IP[a, a]).sqrt() for a in range(nA)]
    G0 = 2*IP
    oldp = ctx.prec; ctx.prec = 128
    G1 = arb_mat(nA, nA)
    for a in range(nA):
        for b in range(nA): G1[a, b] = +(4*IP[a, b])
    ctx.prec = oldp
    del IP
    n2 = [arb(0)]*nA
    CH = 192
    for s0 in range(0, len(G), CH):
        ch = G[s0:s0+CH]; m = len(ch)
        Pc = arb_mat(m, nA); C0 = arb_mat(m, nA); x1s = []
        for gi, j in enumerate(ch):
            x0 = 1/arb(d0[j]) - 2; x1s.append(1/arb(d0[j])**2 - 4)
            for ci, cl in enumerate(allc):
                t = ip(('phi', j, None), cl)*cl[2]
                Pc[gi, ci] = t; C0[gi, ci] = t*x0
                n2[ci] += t*t
        PT = Pc.transpose()
        prod = PT*C0
        for a in range(nA):
            for b in range(nA): G0[a, b] = G0[a, b] + prod[a, b]
        del prod, C0
        ctx.prec = 128
        PT1 = arb_mat(nA, m); C1 = arb_mat(m, nA)
        for gi in range(m):
            for ci in range(nA):
                v = +Pc[gi, ci]; PT1[ci, gi] = v; C1[gi, ci] = v*x1s[gi]
        prod = PT1*C1
        for a in range(nA):
            for b in range(nA): G1[a, b] = G1[a, b] + prod[a, b]
        ctx.prec = oldp
        del prod, PT1, C1, Pc, PT
    nr = [x.sqrt() for x in n2]
    log({'nA': nA, 't_assembly': round(time.time() - ta, 1), 'rss': rss()})

    def margin(sidx_B, sidx_U, Mdiag, eps, rigorous):
        epsb = arb(eps) if eps != 0 else arb(0)
        kap = 1/(arb(dmin)**2*(arb(dmin) - epsb)) + 16
        nB = len(sidx_B); nU = len(sidx_U)
        e2 = epsb*epsb
        def gv(ia, ib):
            v = G0[ia, ib] + epsb*G1[ia, ib]
            if eps != 0: v = v + arb(0, (e2*(kap*nr[ia]*nr[ib] + 16*ipn[ia]*ipn[ib])).upper())
            return v
        X = arb_mat(nB, nB); Y = arb_mat(nU, nB); Sm = arb_mat(nU, nU)
        for a in range(nB):
            for b in range(nB): X[a, b] = gv(sidx_B[a], sidx_B[b])
        for a in range(nU):
            for b in range(nB): Y[a, b] = gv(sidx_U[a], sidx_B[b])
            for b in range(nU): Sm[a, b] = gv(sidx_U[a], sidx_U[b]) + ((1/Mdiag[a]) if a == b else 0)
        if not rigorous:
            Z = Sm.solve(Y, algorithm="approx"); Wm = X - Y.transpose()*Z
        elif os.environ.get('SOLVEMODE') == 'resid':
            # a-posteriori: S*^-1 Y* = Zt + S*^-1 (Y* - S* Zt); exact S* >= diag(1/Mdiag) (Gram_UU is PSD), so ||S*^-1|| <= max Mdiag;
            # |(Y*^T S*^-1 R*)_ab| <= ||Y_a|| max(M) ||R_b||, added as a ball to W = X - Y^T Zt
            Za = Sm.solve(Y, algorithm="approx")
            Zt = arb_mat(nU, nB)
            for i in range(nU):
                for j in range(nB): Zt[i, j] = Za[i, j].mid()
            del Za
            Res = Y - Sm*Zt
            Wm = X - Y.transpose()*Zt
            mM = max(Mdiag[a].upper() for a in range(nU))
            ny = [sum((abs(Y[i, a]).upper()**2 for i in range(nU)), arb(0)).sqrt().upper() for a in range(nB)]
            nrr = [sum((abs(Res[i, b]).upper()**2 for i in range(nU)), arb(0)).sqrt().upper() for b in range(nB)]
            if os.environ.get('CERTDEBUG'): log({'dbg': 'resid', 'mM': float(mM), 'ny_max': float(max(ny)), 'nrr_max': float(max(nrr)), 'Zt_max': float(max(abs(Zt[i, j].mid()) for i in range(nU) for j in range(nB)))})
            del Res, Zt
            for a in range(nB):
                for b in range(nB): Wm[a, b] = Wm[a, b] + arb(0, (ny[a]*mM*nrr[b]).upper())
        else:
            Z = Sm.solve(Y); Wm = X - Y.transpose()*Z
        del X, Y, Sm
        IW = arb_mat(nB, nB)
        for a in range(nB):
            for b in range(nB): IW[a, b] = (1 if a == b else 0) - Wm[a, b]
        R = arb_mat(nB, nB)
        for a in range(nB):
            for b in range(nB): R[a, b] = IW[a, b].mid()
        ev = sorted([e.real.mid() for e in acb_mat(R).eig(algorithm="approx")], key=lambda z: float(z))
        lmin = ev[0]
        if not rigorous or not (lmin > 0): return lmin, False
        nu = lmin*arb('0.999'); s0 = lmin*arb('0.9995')
        A = [[R[i, j] - (s0 if i == j else 0) for j in range(nB)] for i in range(nB)]
        Lc = [[arb(0)]*nB for _ in range(nB)]
        for i in range(nB):
            for j in range(i+1):
                s = A[i][j]
                for kk in range(j): s -= Lc[i][kk]*Lc[j][kk]
                if i == j:
                    if not (s.mid() > 0): return lmin, False
                    Lc[i][j] = s.mid().sqrt().mid()
                else:
                    Lc[i][j] = (s/Lc[j][j]).mid()
        Lm = arb_mat(nB, nB)
        for i in range(nB):
            for j in range(i+1): Lm[i, j] = Lc[i][j]
        Fm = IW - Lm*Lm.transpose()
        emax = arb(0)
        for i in range(nB):
            rs = arb(0)
            for j in range(nB):
                v = Fm[i, j] - (s0 if i == j else 0); rs += abs(v).upper()
            if rs.upper() > emax.upper(): emax = rs.upper()
        if os.environ.get('CERTDEBUG'): log({'dbg': 'lemma55', 'emax': float(emax), 's0-nu': float(s0 - nu), 'lmin': float(lmin)})
        return lmin, bool(emax <= s0 - nu)

    for s in sets:
        ts = time.time()
        D = s.get('del', []); Ad = [(a[0], a[1]) for a in s.get('add', [])]
        sB = list(range(nBb)) + [iBx[k] for k in D]
        sU = list(range(iU0, iU0 + nUb)) + [iUx[a] for a in Ad]
        Md = Mbase + [arb(2)]*len(Ad)
        lg = gal.get(s['name'])
        m0, _ = margin(sB, sU, Md, 0, False)
        e1 = 0.5*lg if lg else 1e-300
        m1, _ = margin(sB, sU, Md, e1, False)
        slope = (m0 - m1)/arb(e1)
        est = float((m0/slope).mid())
        cert = None; tries = []
        for f in (0.985, 0.95, 0.85, 0.6):
            e = f*est
            m, ok = margin(sB, sU, Md, e, True)
            tries.append((f, float(m.mid()) if hasattr(m, 'mid') else float(m), ok))
            if ok and m > 0: cert = e; break
        log({'set': s['name'], 'del': D, 'add': s.get('add', []), 'margin0': float(m0), 'eps_extrap': est, 'eps_cert': cert,
             'galerkin': lg, 'ratio_cert_gal': (cert/lg if (cert and lg) else None), 'tries': tries, 't': round(time.time() - ts, 1), 'rss': rss()})
    log({'wall': round(time.time() - t0, 1), 'rss': rss()})

if __name__ == '__main__':
    main()
