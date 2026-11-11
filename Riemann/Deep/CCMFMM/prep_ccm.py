# prep_slt.py (explore28): node sets, weights and Cauchy sources of Chekmachev's minorant (exactly as WeilZeta/cert_zeta.py),
# written as particle data for slt.c (structured LDL^T of the Cauchy-like matrix J'' = G + diag(delta)).
#   index a: point y_a, charges (u_a, w_a), diagonal; G_ab = (u_a w_b - w_a u_b)/(y_a - y_b) (a != b) = <c_a, P0^{-1} c_b>.
#   node k: y = k^2, u = beta_k, w = beta_k F(k^2), diag = (L/2)/deff + beta_k^2 F2(k^2), delta = 1/w_k
#   pole v: y = y_v, u = beta_v, w = beta_v F(y_v) + kappa/deff, diag = |v|^2/deff + beta_v^2 F2(y_v), delta = +1/2 (even) or -1/2 (odd)
#   phi j (dl_j <= 0): y = 4 j^2, u = 0, w = A_j/deff, diag = 1/deff, delta = 1/(dl_j - deff)
#   F(y) = sum_{j in G} gam_j/(y - 4j^2), F2(y) = sum gam_j/(y - 4j^2)^2, gam_j = (1/dl_j - 1/deff) A_j^2  (computed in slt.c)
# Usage: prep_slt.py c sector eps bits out [K2] [adv]    (adv = 1: adversarial moved-zero-pair symbol, see adversarial.py)
# CCM (ccmfmm): sector 'defl' = the even sector plus ../WeilSimple's deflation column M g g^T in the positive part, g = sum_l g_l phi_l
#   read from CCM_G (gstate.py JSON; exact dyadics), M = CCM_M (default 1). Writes out.brd for slt.c (env SLT_BRD):
#   "n J''_gg delta_g", then r_a = <g, P0^{-1} c_a> for every particle a (order of out.idx):
#   node k: beta_k Phi(k^2), pole v: beta_v Phi(y_v), bad phi_j: g_j/deff (j <= NG) else 0, Phi(y) = sum_l g_l A_l/(p_l (y - 4 l^2)),
#   p_l = dl_l (l in G) or deff; J''_gg = sum_l g_l^2/p_l + 1/M, delta_g = 1/M. Lines marked CCM.
import sys, os, time, math, json
_h = os.path.dirname(os.path.abspath(__file__))
for _p in (os.path.join(_h, '..', 'WeilZeta'), os.path.join(_h, '..', '..', 'Riemann', 'Deep', 'WeilZeta')):
    if os.path.exists(os.path.join(_p, 'cert_zeta.py')):
        sys.path.insert(0, _p); break
from flint import arb, acb, ctx, fmpq, fmpz as fmpz_
from cert_zeta import SymZeta, scan_end, rscan, fscan, parse_c

def main():
    cq = parse_c(sys.argv[1]); sector = sys.argv[2]; epsstr = sys.argv[3]; ctx.prec = int(sys.argv[4]); out = sys.argv[5]
    defl = (sector == 'defl')   # CCM
    if defl: sector = 'even'
    adv = len(sys.argv) > 7 and sys.argv[7] == '1'
    if epsstr == '0': eps = arb(0)
    elif epsstr.startswith('2^-'): eps = arb(2)**(-int(epsstr[3:]))
    else: eps = arb(epsstr)
    t0 = time.time()
    S = SymZeta(cq); L = S.L; pi = S.pi
    cf = int(cq.p)/int(cq.q)
    K2 = max(201, math.ceil(5*cf*cf/3)); K2 += (K2 % 2 == 0)
    if len(sys.argv) > 6 and sys.argv[6] != '-': K2 = int(sys.argv[6]); K2 += (K2 % 2 == 0)
    J = K2; t1 = arb(1); t1f = 1.0
    ks = scan_end(S, K2, t1f); kse = scan_end(S, 2*J, 0)
    rexe = os.environ.get('WZ_RSCAN')
    cache = os.path.join(os.environ.get('WZ_SCAN', '.'), f"scan_zeta_{str(cq).replace('/', 'o')}_{K2}" + ("_rig" if rexe else "_flt") + ".json")
    if os.path.exists(cache):
        sc = json.load(open(cache)); assert sc['ks'] == ks and sc['kse'] == kse
        selo, sele = sc['selo'], sc['sele']
    else:
        if rexe: fo, fe, sinfo = rscan(S, K2, J, ks, kse, rexe, t1f)
        else: fo, fe = fscan(S, K2, J, ks, kse, t1f); sinfo = dict(method='float')
        oldp = ctx.prec; ctx.prec = 96; S96 = SymZeta(cq)
        selo = [k for k in fo if S96.mu(k).lower() < t1f]
        sele = [j for j in fe if S96.mu(2*j).lower() < 0]
        ctx.prec = oldp
        json.dump(dict(ks=ks, kse=kse, selo=selo, sele=sele, flagged=[len(fo), len(fe)], scan=sinfo), open(cache, 'w'))
    tscan = time.time() - t0
    dmu = None
    if adv:   # explore27's moved pair: gamma_1, gamma_2 -> (gamma_1 + gamma_2)/2 +- i delta, delta = 0.3 (balls via acb)
        import mpmath
        mpmath.mp.prec = ctx.prec + 20
        g1 = arb(str(mpmath.zetazero(1).imag)); g2 = arb(str(mpmath.zetazero(2).imag)); gm = (g1 + g2)/2; dd = arb('0.3')
        def zw(gam, r):   # int_{-L}^{L} e^{i gam x} e^{-i r x} dx = 2 sin((gam - r) L)/(gam - r)
            z = gam - r
            return 2*(z*L).sin()/z
        def dmu(r):
            tot = acb(0)
            for s in (1, -1):
                tot -= acb(zw(s*g1, r)) + acb(zw(s*g2, r))
                tot += zw(acb(s*gm, dd), acb(r)) + zw(acb(s*gm, -dd), acb(r))
            return tot.real
    NBO, NBEK, NBEJ = 32, 1, 3
    Kset = set(range(1, K2+1, 2)) | set(selo); Eset = set(range(0, J+1)) | set(sele)
    for k in selo:
        for j in range((k-1)//2 - NBO, (k-1)//2 + NBO + 1):
            if j >= 0: Eset.add(j)
    for j in sele:
        for jj in range(j-NBEJ, j+NBEJ+1):
            if jj >= 0: Eset.add(jj)
    lean = os.environ.get('LEAN') == '1'
    mus = {}
    def mu(k):
        if lean: return S.mu(k) + (dmu(S.r(k)) if dmu is not None else 0)
        if k not in mus:
            m = S.mu(k)
            if dmu is not None: m = m + dmu(S.r(k))
            mus[k] = m
        return mus[k]
    # neighbour rule (Chekmachev Def. 6.1); it only adds indices (always valid), so in lean mode its sign tests use 128-bit balls
    oldp = ctx.prec
    if lean: ctx.prec = 128; S128 = SymZeta(cq)
    def musign(k):
        if not lean: return mu(k).mid() < 0
        m = S128.mu(k)
        if m.contains(0): ctx.prec = oldp; r = S.mu(k).mid() < 0; ctx.prec = 128; return r
        return m < 0
    changed = True; done = set()
    while changed:
        changed = False
        for j in sorted(Eset):
            if j in done or 2*j <= K2: continue
            done.add(j)
            if musign(2*j):
                for k in range(2*j-2*NBEK+1, 2*j+2*NBEK, 2):
                    if k > 0 and k not in Kset: Kset.add(k); changed = True
                for jj in range(j-NBEJ, j+NBEJ+1):
                    if jj >= 0 and jj not in Eset: Eset.add(jj); changed = True
    ctx.prec = oldp
    del done
    if sector == 'odd': Eset.discard(0)
    Kl = sorted(Kset); El = sorted(Eset); del Kset, Eset
    deff = t1/2 - eps
    assert deff > 0
    w = {k: ((mu(k) - t1)/L).lower() for k in Kl}
    Kp = [k for k in Kl if w[k] > 0]; Km = [k for k in Kl if not (w[k] > 0)]
    assert all(w[k] != 0 for k in Km)
    def Nj(j): return (1/L).sqrt() if (j == 0 and sector == 'even') else (2/L).sqrt()
    sgn = lambda j: 1 if j % 2 == 0 else -1
    sgk = lambda k: 1 if ((k-1)//2) % 2 == 0 else -1
    sL4, cL4, sL2 = (L/4).sinh(), (L/4).cosh(), (L/2).sinh()
    yv = -(L*L)/(4*pi*pi)
    if sector == 'even':
        A = lambda j: Nj(j)*(L/pi)*sgn(j); bk = lambda k: 2*k*sgk(k); bv = -(L/pi)*sL4
        kap = (L/pi)*cL4; vv = sL2 + L/2; dpole = arb(1)/2       # pole column in U, weight +2
    else:
        A = lambda j: Nj(j)*sgn(j)*4*j*(L/pi); bk = lambda k: sgk(k); bv = sL4
        kap = (L/pi)**2*cL4; vv = sL2 - L/2; dpole = -arb(1)/2   # pole column in B, weight -2
    nd = int(ctx.prec*0.30103) + 12
    S_ = lambda x: x.str(nd, radius=True, more=True).replace(' ', '')
    nsrc = 0; Jm = []
    with open(out + '.src', 'w') as fs:
        for j in El:
            dlj = ((mu(2*j) + t1)/2 - eps).lower()
            if dlj > 0:
                fs.write(f"{4*j*j} {S_((1/dlj - 1/deff)*A(j)**2)}\n"); nsrc += 1
            else: Jm.append((j, dlj))
    # indices: U block (K+ nodes, even pole), then B block (bad phi_j, K- nodes, odd pole)
    rows = []
    for k in Kp: rows.append(f"N 0 {k*k} {bk(k)} {S_((L/2)/deff)} {S_(1/w[k])}")
    if sector == 'even': rows.append(f"V 0 {S_(yv)} {S_(bv)} {S_(kap/deff)} {S_(vv/deff)} {S_(dpole)}")
    nU = len(rows)
    for j, dlj in Jm: rows.append(f"P 1 {4*j*j} {S_(A(j)/deff)} {S_(1/deff)} {S_(1/(dlj - deff))}")
    for k in Km: rows.append(f"N 1 {k*k} {bk(k)} {S_((L/2)/deff)} {S_(1/w[k])}")
    if sector == 'odd': rows.append(f"V 1 {S_(yv)} {S_(bv)} {S_(kap/deff)} {S_(vv/deff)} {S_(dpole)}")
    with open(out + '.idx', 'w') as fi:
        fi.write(f"{len(rows)} {nU} {nsrc}\n" + "\n".join(rows) + "\n")
    gin = None   # CCM: border of g
    if defl:
        gj = json.load(open(os.environ['CCM_G'])); assert gj['c'] == str(cq)
        gc = [arb(fmpz_(m)) * arb(2)**int(e) if int(e) >= 0 else arb(fmpz_(m)) / arb(2)**(-int(e)) for m, e in gj['g']]
        NG = len(gc) - 1; Mw = arb(os.environ.get('CCM_M', '1'))
        assert all(l in El for l in range(NG + 1))
        dls = {l: ((mu(2*l) + t1)/2 - eps).lower() for l in range(NG + 1)}
        pl = [(dls[l] if dls[l] > 0 else deff) for l in range(NG + 1)]
        cl = [gc[l]*A(l)/pl[l] for l in range(NG + 1)]
        def Phi(y): return sum((cl[l]/(y - 4*l*l) for l in range(NG + 1)), arb(0))
        rb = []
        for k in Kp: rb.append(bk(k)*Phi(arb(k*k)))
        if sector == 'even': rb.append(bv*Phi(yv))
        for j, dlj in Jm: rb.append(gc[j]/deff if j <= NG else arb(0))
        for k in Km: rb.append(bk(k)*Phi(arb(k*k)))
        assert len(rb) == len(rows)
        Jgg = sum((gc[l]*gc[l]/pl[l] for l in range(NG + 1)), arb(0)) + 1/Mw
        with open(out + '.brd', 'w') as fb:
            fb.write(f"{len(rows)} {S_(Jgg)} {S_(1/Mw)}\n" + "\n".join(S_(x) for x in rb) + "\n")
        gin = dict(NG=NG, M=str(Mw.mid().str(5)), U=gj['U'], U_up=gj['U_up'], Jgg=Jgg.str(10), g0=gc[0].str(8), gnorm2=sum((x*x for x in gc), arb(0)).str(12))
    info = dict(c=str(cq), sector=('defl' if defl else sector), g=gin, eps=epsstr, prec=ctx.prec, K2=K2, kscan=ks, kscan_e=kse, sel_odd=len(selo), sel_even=len(sele),
                Kp=len(Kp), Km=len(Km), Jm=len(Jm), G=nsrc, order=len(rows) - nU, nU=nU, n=len(rows), adv=adv, tscan=round(tscan), tprep=round(time.time() - t0))
    print(json.dumps(info), flush=True)

if __name__ == '__main__':
    main()
