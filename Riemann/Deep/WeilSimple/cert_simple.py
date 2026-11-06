# cert_simple.py (weilsimple): cert_zeta.py of ../WeilZeta (copied unchanged except the parts marked WS) plus a deflated
# even certificate: sector 'defl' proves Q(f) - eps ||f||^2 + M <f, g>^2 >= 0 for all even f in L^2(-L/2, L/2), with g the
# approximate ground state of the even block of QW_lambda^NG (weil_matrix.py; env WS_NG, WS_M), and prints
# U = Q(g)/||g||^2 (Arb). Then: even f orthogonal to g have Q(f) >= eps ||f||^2 (Courant-Fischer: lambda_2^even >= eps).
# Original header of cert_zeta.py:
# cert_zeta.py: Chekmachev's periodized lattice minorant (Zenodo 23848597, Thm 5.10 + Def. 6.1) for the Weil form of ZETA.
# Port of Riemann/Deep/WeilChi/cert_chi.py (rigorous scan rscan.c, blocked assembly) back to zeta: trivial character, constant
# -log pi, archimedean b_n = 2n + 1/2, prime powers Lambda(n)/sqrt n (n < c), and the rank-one pole column of
# Riemann/Deep/ChekAudit/evidence/cert2.py (+2<f,cosh(x/2)>^2 for even f, -2<f,sinh(x/2)>^2 for odd f).
# Woodbury step: Chekmachev's Lemma 5.3 (W <= W_up = X - Y^T Z - Z^T Y + Z^T S Z for ANY Z, since S > 0) with Z an approximate
# solve, so no rigorous (memory-hungry) arb_mat solve is needed. Derivation: DERIVATION.md.
#   Q(f) = pole(f) + (1/2L) sum_k mu(k) |F(r_k)|^2,  r_k = pi k/L,  L = ln c,  f in L^2(-L/2, L/2),
#   mu(k) = -log pi + Re psi(1/4 + i r_k/2) + 2 (-1)^k sum_n c^{-b_n} b_n/(b_n^2 + r_k^2) - 2 sum_{m<c} Lambda(m) m^{-1/2} cos(r_k ln m).
# Usage: cert_zeta.py c sector eps bits [K2]      c = integer or p/q; eps = 0 or 2^-k or decimal; sector even|odd
# Env: WZ_RSCAN=<rscan binary> (rigorous scan; else numpy float scan), WZ_SCAN=<cache dir>, WZ_T1, WZ_NBO, WZ_NBEK, WZ_NBEJ, WZ_CHOL_BS
import sys, time, math, random, json, os
import numpy as np, mpmath as mp
from flint import arb, acb, arb_mat, acb_mat, ctx, fmpq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from weil_matrix import WeilMat   # WS

def ground_vec(T, n, its=8):   # WS: inverse iteration on mid(T); returns exact dyadic coefficients (midpoints)
    rnd = random.Random(7); v = arb_mat(n, 1)
    for i in range(n): v[i, 0] = arb(rnd.uniform(-1, 1))
    for it in range(its):
        y = T.solve(v, algorithm='approx')
        nrm = sum((y[i, 0]*y[i, 0] for i in range(n)), arb(0)).sqrt()
        for i in range(n): v[i, 0] = (y[i, 0]/nrm).mid()
    return [v[i, 0] for i in range(n)]

def rayleigh(T, g):   # WS: rigorous enclosure of g^T T g / g^T g
    n = len(g); num = arb(0)
    for i in range(n):
        s = arb(0)
        for j in range(n): s += T[i, j]*g[j]
        num += g[i]*s
    return num, sum((x*x for x in g), arb(0))

def ppow(n):   # (p, k) if n = p^k else None
    for p in range(2, n+1):
        if n % p == 0:
            m = n; e = 0
            while m % p == 0: m //= p; e += 1
            return (p, e) if m == 1 else None
    return None

class SymZeta:
    def __init__(s, cq):
        s.key, s.q, s.a = 'zeta', 1, 0
        s.cq = cq; s.c = arb(cq); s.L = s.c.log(); s.pi = arb.pi()
        s.const = -s.pi.log()
        s.atoms = []           # (n, ln n, Lambda(n)/sqrt n), n < c
        n = 2
        while fmpq(n) < cq:
            pe = ppow(n)
            if pe is not None: s.atoms.append((n, arb(n).log(), arb(pe[0]).log()/arb(n).sqrt()))
            n += 1
        s.sig2 = 2*sum((abs(w) for _, _, w in s.atoms), arb(0))
        s.tail = []; nn = 0
        while True:
            b = arb(2*nn) + arb(1)/2; t = s.c**(-b)
            s.tail.append((b, t)); nn += 1
            if t < arb(2)**(-(ctx.prec + 40)) and nn > 2: break
        bN = arb(2*nn) + arb(1)/2
        s.trem = s.c**(-bN)/(bN*(1 - s.c**(-2)))     # sum_{n>=N} c^{-b_n} b_n/(b_n^2+r^2) <= this
    def r(s, k): return s.pi*k/s.L
    def h(s, r): return s.const + acb(arb(1)/4, r/2).digamma().real
    def tailsum(s, r):
        r2 = r*r; tot = arb(0)
        for b, t in s.tail: tot += t*b/(b*b + r2)
        return tot + arb(0, s.trem.upper())
    def mu(s, k):
        r = s.r(k); ps = arb(0)
        for _, ln, w in s.atoms: ps += w*(r*ln).cos()
        return s.h(r) + 2*(1 if k % 2 == 0 else -1)*s.tailsum(r) - 2*ps
    def monotone(s, k):   # lower bound of mu(k') for all k' >= k >= 1 (h increasing, b/(b^2+r^2) decreasing in r)
        r = s.r(k)
        return s.h(r) - s.sig2 - 2*(s.tailsum(r) + s.trem)

def scan_end(S, start, thr):
    ok = lambda k: S.monotone(k+1).lower() >= thr
    if ok(start): return start
    lo, hi = start, start*2
    while not ok(hi): lo, hi = hi, hi*2
    while hi - lo > 1:
        m = (lo + hi)//2
        if ok(m): hi = m
        else: lo = m
    return hi

def fscan(S, K2, J, ks, kse, t1f=1.0):
    # float64 pre-scan (NOT rigorous: libm cos), exact uint64 phases; flags at margin 1e-6, confirmed later in Arb
    mp.mp.dps = 40
    cmp_ = mp.mpf(int(S.cq.p))/mp.mpf(int(S.cq.q)); Lm = mp.log(cmp_)
    B = np.array([int(mp.nint(mp.log(n)/(2*Lm)*mp.mpf(2)**64)) % 2**64 for n, _, _ in S.atoms], dtype=np.uint64)
    W = np.array([float(w.mid()) for _, _, w in S.atoms]); Lf = float(Lm); cst = float(S.const.mid())
    bn = np.array([2*n + 0.5 for n in range(8)]); cbn = np.array([float(cmp_**(-(2*n + mp.mpf(1)/2))) for n in range(8)])
    def mu(k):
        k = k.astype(np.uint64); r = np.pi*k.astype(np.float64)/Lf; w = 0.25 + 0.5j*r
        psi = np.log(w) - 1/(2*w) - 1/(12*w**2) + 1/(120*w**4) - 1/(252*w**6) + 1/(240*w**8) - 1/(132*w**10)
        sgn = np.where(k % 2 == 0, 1.0, -1.0)
        T = 2*sgn*sum(cbn[i]*bn[i]/(bn[i]**2 + r*r) for i in range(8))
        s = np.zeros_like(r)
        for i in range(len(B)): s += W[i]*np.cos(2*np.pi*((k*B[i]).astype(np.float64)/2.0**64))
        return psi.real + cst + T - 2*s
    selo, sele = [], []; CHK = 2_000_000
    k0 = K2 + 2
    while k0 <= ks:
        k = np.arange(k0, min(k0 + 2*CHK, ks + 1), 2, dtype=np.int64); m = mu(k); selo.extend(k[m < t1f + 1e-6].tolist()); k0 += 2*CHK
    j0 = J + 1
    while 2*j0 <= kse:
        j = np.arange(j0, min(j0 + CHK, kse//2 + 1), dtype=np.int64); m = mu(2*j); sele.extend(j[m < 1e-6].tolist()); j0 += CHK
    return selo, sele

def rscan(S, K2, J, ks, kse, exe, t1f=1.0):
    # RIGOROUS pre-scan (rscan.c; proof of the error bound in RIGOR.md): every node k in (K2, ks] (odd) or (2J, kse] (even) with
    # mu(k) < threshold is printed by rscan.c, because mu(k) >= mu_c(k) - E and E < margin.  Returns flagged lists and the bound data.
    import subprocess, mpmath as mpm
    cq = S.cq
    oldp = ctx.prec; ctx.prec = 192
    S2 = SymZeta(cq); pi = arb.pi(); L = arb(cq).log(); two64 = arb(2)**64
    u = 2.0**-53
    def dn(x): return math.nextafter(float(x.lower()), -math.inf)
    def up(x): return math.nextafter(float(x.upper()), math.inf)
    tab = []
    for a in range(65536):
        ang = 2*pi*a/65536
        tab.append(f"{float(ang.cos().mid()).hex()} {float(ang.sin().mid()).hex()}")
    K = float((2*pi/two64).mid())
    mpm.mp.dps = 60; Lm = mpm.log(mpm.mpf(int(cq.p))/mpm.mpf(int(cq.q)))
    at = []; emax = arb(0); dW = arb(0); sigW = 0.0
    for n, _, w in S2.atoms:
        Bn = int(mpm.floor(mpm.log(n)/(2*Lm)*mpm.mpf(2)**64 + mpm.mpf(1)/2))
        en = abs(arb(n).log()/(2*L)*two64 - Bn)
        if en.upper() > emax.upper(): emax = arb(en.upper())
        Wn = float(w.mid()); dW += abs(arb(Wn) - w).upper(); sigW += abs(Wn)
        at.append(f"{Bn % 2**64} {Wn.hex()}")
        assert 0 <= Bn < 2**64
    assert emax < arb('0.5000001')
    def blocks(k0, k1):
        out = []; Hmax = 0.0; Tmax = 0.0
        if k0 > k1: return out, Hmax, Tmax
        k1 -= (k1 - k0) % 2
        k = k0
        while k <= k1:
            m = max(0, int(k*0.0005)); kb = min(k1, k + 2*m)
            r = S2.r(k); H = dn(S2.h(r)); T = up(2*S2.tailsum(r) + 2*S2.trem)
            out.append(f"{k} {kb} {H.hex()} {T.hex()}"); Hmax = max(Hmax, abs(H)); Tmax = max(Tmax, abs(T))
            k = kb + 2
        return out, Hmax, Tmax
    bo, Ho, To = blocks(K2 + 2, ks); be, He, Te = blocks(2*(J + 1), kse)
    na = len(S2.atoms); kmax = max(ks, kse)
    gam = na*u/(1 - na*u)
    # a-priori bound (RIGOR.md): |cos_c - cos(theta_used)| <= 1e-15; phase error <= 2 pi k e_n 2^-64; sum and final rounding
    E = 2*(sigW*(1 + u)*(gam*(1 + u) + u + 1e-15) + float(dW.upper())*1.0001 + sigW*1.0001*2*math.pi*kmax*float(emax.upper())*2.0**-64)
    E += 3*u*(max(Ho, He) + max(To, Te) + 2*sigW*1.0001)
    E *= 1.01
    margin = 1e-6
    assert E < margin/10, E
    ctx.prec = oldp
    hdr = "\n".join(tab) + f"\n{K.hex()}\n{na}\n" + "\n".join(at) + "\n"
    def run(bl, thr):
        if not bl: return []
        inp = hdr + f"{float(thr).hex()} {margin.hex()}\n{len(bl)}\n" + "\n".join(bl) + "\n"
        out = subprocess.run([exe], input=inp, capture_output=True, text=True, check=True)
        return [int(v) for v in out.stdout.split()]
    selo = run(bo, t1f)
    sele = [k//2 for k in run(be, 0.0)]
    info = dict(method='rigorous', E=E, margin=margin, emax=float(emax.upper()), blocks=[len(bo), len(be)], kmax=kmax)
    return selo, sele, info

def parse_c(s):
    if '/' in s: p, q = s.split('/'); return fmpq(int(p), int(q))
    return fmpq(int(s))

def chol_python(R, n, s0):
    # approximate Cholesky of mid(R) - s0 Id (exact arb entries), plain triple loop; None if a pivot is not positive
    Lc = [[None]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1):
            s = R[i, j] - (s0 if i == j else 0)
            Li = Lc[i]; Lj = Lc[j]
            for kk in range(j): s -= Li[kk]*Lj[kk]
            if i == j:
                if not (s > 0): return None
                Li[j] = s.sqrt().mid()
            else: Li[j] = (s/Lj[j]).mid()
    Lm = arb_mat(n, n)
    for i in range(n):
        row = Lc[i]
        for j in range(i+1): Lm[i, j] = row[j]
        Lc[i] = None
    return Lm

def chol_block(R, n, s0, bs):
    # right-looking block Cholesky of mid(R) - s0 Id, floating point at the working precision (no error bounds needed:
    # the residual is enclosed in Arb afterwards). Panel solves and the trailing update run in C (arb_mat ops).
    Sr = [[(R[i, j] - (s0 if i == j else 0)).mid() for j in range(i+1)] for i in range(n)]   # lower triangle, updated in place
    Lm = arb_mat(n, n)
    for i0 in range(0, n, bs):
        i1 = min(n, i0 + bs); b = i1 - i0; m = n - i1
        L11 = arb_mat(b, b)
        for i in range(b):
            for j in range(i+1):
                s = Sr[i0+i][i0+j]
                for kk in range(j): s -= L11[i, kk]*L11[j, kk]
                if i == j:
                    if not (s > 0): return None
                    L11[i, j] = s.sqrt().mid()
                else: L11[i, j] = (s/L11[j, j]).mid()
        for i in range(b):
            for j in range(i+1): Lm[i0+i, i0+j] = L11[i, j]
        if m == 0: break
        PT = arb_mat(b, m)                      # panel transposed: PT[c, r] = S[i1 + r, i0 + c]
        for r in range(m):
            row = Sr[i1+r]
            for c in range(b): PT[c, r] = row[i0+c]
        X = L11.solve(PT, algorithm="approx").mid()    # L11 X = PT  =>  L21 = X^T
        for r in range(m):
            for c in range(b): Lm[i1+r, i0+c] = X[c, r]
        U = X.transpose()*X                     # L21 L21^T, m x m (its midpoints are used below; no extra copy)
        for r in range(m):
            row = Sr[i1+r]
            for c in range(r+1): row[i1+c] = (row[i1+c] - U[r, c]).mid()
        del PT, X, U
    return Lm

def main():
    cq = parse_c(sys.argv[1]); sector = sys.argv[2]; epsstr = sys.argv[3]; ctx.prec = int(sys.argv[4])
    defl = (sector == 'defl')   # WS
    if defl: sector = 'even'
    if epsstr == '0': eps = arb(0)
    elif epsstr.startswith('2^-'): eps = arb(2)**(-int(epsstr[3:]))
    else: eps = arb(epsstr)
    t0 = time.time()
    S = SymZeta(cq); L = S.L; pi = S.pi
    cf = int(cq.p)/int(cq.q)
    K2 = max(201, math.ceil(5*cf*cf/3)); K2 += (K2 % 2 == 0)
    if len(sys.argv) > 5 and sys.argv[5] != '-': K2 = int(sys.argv[5]); K2 += (K2 % 2 == 0)
    J = K2
    T1s = os.environ.get('WZ_T1', '1'); t1 = arb(T1s); t1f = float(T1s)
    ks = scan_end(S, K2, t1f); kse = scan_end(S, 2*J, 0)
    rexe = os.environ.get('WZ_RSCAN')
    cache = os.path.join(os.environ.get('WZ_SCAN', '.'), f"scan_zeta_{str(cq).replace('/', 'o')}_{K2}" + ("" if T1s == '1' else f"_t{T1s}") + ("_rig" if rexe else "_flt") + ".json")
    if os.path.exists(cache):
        sc = json.load(open(cache)); assert sc['ks'] == ks and sc['kse'] == kse
        selo, sele, sinfo = sc['selo'], sc['sele'], sc.get('scan', {})
    else:
        sinfo = dict(method='float')
        if rexe: fo, fe, sinfo = rscan(S, K2, J, ks, kse, rexe, t1f)
        else: fo, fe = fscan(S, K2, J, ks, kse, t1f)
        oldp = ctx.prec; ctx.prec = 96; S96 = SymZeta(cq)
        selo = [k for k in fo if S96.mu(k).lower() < t1f]
        sele = [j for j in fe if S96.mu(2*j).lower() < 0]
        ctx.prec = oldp
        json.dump(dict(ks=ks, kse=kse, selo=selo, sele=sele, flagged=[len(fo), len(fe)], scan=sinfo), open(cache, 'w'))
    tscan = time.time() - t0
    NBO = int(os.environ.get('WZ_NBO', '32')); NBEK = int(os.environ.get('WZ_NBEK', '1')); NBEJ = int(os.environ.get('WZ_NBEJ', '3'))
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
    gc = None; Uenc = None; NG = None   # WS: ground-state vector g and U = Q(g)/||g||^2
    if defl:
        NG = int(os.environ.get('WS_NG', '80')); Mw = arb(os.environ.get('WS_M', '1'))
        WM = WeilMat(cq); Tg, _ = WM.block(NG, 'even')
        gc = ground_vec(Tg, NG + 1)
        num, gg = rayleigh(Tg, gc); Uenc = num/gg
        del Tg
        Gset = set(G)
        assert all(l in dl for l in range(NG + 1))
        pl = [(dl[l] if l in Gset else deff) for l in range(NG + 1)]
    tU = time.time() - t0
    def Nj(j): return (1/L).sqrt() if (j == 0 and sector == 'even') else (2/L).sqrt()
    sgn = lambda j: 1 if j % 2 == 0 else -1
    sgk = lambda k: 1 if ((k-1)//2) % 2 == 0 else -1
    sL4, cL4, sL2 = (L/4).sinh(), (L/4).cosh(), (L/2).sinh()
    yv = -(L*L)/(4*pi*pi)
    if sector == 'even':     # phi_j = N_j cos(2 pi j x/L), e_k = cos(pi k x/L), v = cosh(x/2) (pole term +2<f,v>^2)
        A = lambda j: Nj(j)*(L/pi)*sgn(j); bk = lambda k: arb(2*k*sgk(k)); bv = -(L/pi)*sL4
        ip_node_v = lambda k: 2*S.r(k)*sgk(k)*cL4/(S.r(k)**2 + arb(1)/4); vv = sL2 + L/2
    else:                    # phi_j = (2/L)^{1/2} sin(2 pi j x/L), e_k = sin(pi k x/L), v = sinh(x/2) (pole term -2<f,v>^2)
        A = lambda j: Nj(j)*sgn(j)*4*j*(L/pi); bk = lambda k: arb(sgk(k)); bv = sL4
        ip_node_v = lambda k: sgk(k)*cL4/(S.r(k)**2 + arb(1)/4); vv = sL2 - L/2
    # columns: (type, index, scale, y, b);  <phi_j, col> = A(j) b/(y - 4 j^2)
    cols = [('phi', j, (deff - dl[j]).sqrt(), None, None) for j in Jm]
    cols += [('node', k, (-w[k]).sqrt(), arb(k*k), bk(k)) for k in Km]
    if sector == 'odd': cols.append(('v', None, arb(2).sqrt(), yv, bv))
    nB = len(cols)
    ucols = [('node', k, arb(1), arb(k*k), bk(k)) for k in Kp]; Md = [w[k] for k in Kp]
    if sector == 'even': ucols.append(('v', None, arb(1), yv, bv)); Md.append(arb(2))
    if defl: ucols.append(('g', None, arb(1), None, None)); Md.append(Mw)   # WS: + M g g^T in the positive part
    allc = cols + ucols; nA = len(allc)
    GA2 = [((1/dl[j] - 1/deff)*A(j)**2, arb(4*j*j)) for j in G]
    F = {}; F2 = {}
    for ci, col in enumerate(allc):
        if col[0] in ('phi', 'g'): continue   # WS
        y = col[3]; f1 = arb(0); f2 = arb(0)
        for (ca, x) in GA2:
            d = 1/(y - x); f1 += ca*d; f2 += ca*d*d
        F[ci] = f1; F2[ci] = f2
    tF = time.time() - t0
    def ipf(a, b):   # <a, b>, unscaled
        ta, ia = a[0], a[1]; tb, ib = b[0], b[1]
        if ta == 'g' or tb == 'g':   # WS: g = sum_l gc[l] phi_l
            if ta == 'g' and tb == 'g': return sum((x*x for x in gc), arb(0))
            o = b if ta == 'g' else a
            return sum((gc[l]*ipf(('phi', l), o) for l in range(len(gc))), arb(0))
        if ta == 'phi' and tb == 'phi': return arb(1) if ia == ib else arb(0)
        if ta == 'phi': return A(ia)*b[4]/(b[3] - arb(4*ia*ia))
        if tb == 'phi': return A(ib)*a[4]/(a[3] - arb(4*ib*ib))
        if ta == 'node' and tb == 'node': return L/2 if ia == ib else arb(0)
        if ta == 'node' and tb == 'v': return ip_node_v(ia)
        if ta == 'v' and tb == 'node': return ip_node_v(ib)
        return vv
    def gram(p, qq):   # <col_p, P0^{-1} col_q> (scaled), P0 = deff Id + sum_G (dl_j - deff) phi_j phi_j^T
        a = allc[p]; b = allc[qq]
        if a[0] == 'g' or b[0] == 'g':   # WS: <g, P0^{-1} b> = sum_l gc[l]/p_l <phi_l, b>, P0 phi_l = p_l phi_l
            if a[0] == 'g' and b[0] == 'g': v = sum((gc[l]*gc[l]/pl[l] for l in range(len(gc))), arb(0))
            else:
                o = b if a[0] == 'g' else a
                v = sum((gc[l]/pl[l]*ipf(('phi', l), o) for l in range(len(gc))), arb(0))
            return v*a[2]*b[2]
        v = ipf(a, b)/deff
        if a[0] != 'phi' and b[0] != 'phi':
            if p == qq: v += a[4]*a[4]*F2[p]
            else: v += a[4]*b[4]*(F[p] - F[qq])/(b[3] - a[3])
        return v*a[2]*b[2]
    nU = len(ucols)
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
    del F, F2, GA2
    tasm = time.time() - t0
    # W = X - Y^T S^{-1} Y <= W_up = X - Y^T Z - Z^T Y + Z^T S Z = X - (Z^T Y)^T + Z^T (S Z - Y)  for any Z (Lemma 5.3; S > 0)
    # Column blocks of B keep the memory at about X + Y + S + Z (no rigorous solve, no full products).
    if nU > 0 and nB > 0:
        WS = int(os.environ.get('WZ_WS', '256')); WB = int(os.environ.get('WZ_WB', '128'))
        Z = arb_mat(nU, nB)
        for b0 in range(0, nB, WS):     # approximate solves S Z_b = Y_b (floating point at the working precision)
            b1 = min(nB, b0 + WS); bs = b1 - b0
            Yb = arb_mat(nU, bs)
            for a in range(nU):
                for b in range(bs): Yb[a, b] = Y[a, b0 + b]
            Zb = Sm.solve(Yb, algorithm="approx")
            for a in range(nU):
                for b in range(bs): Z[a, b0 + b] = Zb[a, b].mid()
            del Yb, Zb
        for b0 in range(0, nB, WB):     # IW[:, b] = e_b - X[:, b] + (Z^T Y)[b, :]^T - (Rz^T Z)[b, :]^T,  Rz = S Z - Y
            b1 = min(nB, b0 + WB); bs = b1 - b0
            Zb = arb_mat(nU, bs); Yb = arb_mat(nU, bs)
            for a in range(nU):
                for b in range(bs): Zb[a, b] = Z[a, b0 + b]; Yb[a, b] = Y[a, b0 + b]
            Rzb = Sm*Zb - Yb
            P1 = Zb.transpose()*Y
            P2 = Rzb.transpose()*Z
            for bb in range(bs):
                for a in range(nB): X[a, b0 + bb] = (1 if a == b0 + bb else 0) - X[a, b0 + bb] + P1[bb, a] - P2[bb, a]
            del Zb, Yb, Rzb, P1, P2
        del Z, Y, Sm
    else:
        for a in range(nB):
            for b in range(nB): X[a, b] = (1 if a == b else 0) - X[a, b]
        del Y, Sm
    IW = X
    tW = time.time() - t0
    maxrad = max((IW[a, b].rad() for a in range(nB) for b in range(nB)), default=arb(0))
    R = IW          # approximate solves and the Cholesky use midpoints only; no copy (memory)
    random.seed(1)
    rig = None; lmin = None; ritz = []; nu = None; emax = None
    BS = int(os.environ.get('WZ_CHOL_BS', '0'))
    if nB > 0:
        # block inverse iteration (p vectors) + Rayleigh-Ritz: estimates of the smallest eigenvalues of mid(IW)
        p = min(6, nB)
        Q = [[arb(random.uniform(-1, 1)) for _ in range(nB)] for _ in range(p)]
        for it in range(6):
            V = arb_mat(nB, p)
            for j in range(p):
                for i in range(nB): V[i, j] = Q[j][i]
            Yv = R.solve(V, algorithm="approx"); Q = []
            for j in range(p):
                v = [Yv[i, j] for i in range(nB)]
                for _ in range(2):
                    for qv in Q:
                        d = sum((v[i]*qv[i] for i in range(nB)), arb(0))
                        v = [v[i] - d*qv[i] for i in range(nB)]
                nrm = sum((vi*vi for vi in v), arb(0)).sqrt()
                Q.append([(vi/nrm).mid() for vi in v])
        V = arb_mat(nB, p)
        for j in range(p):
            for i in range(nB): V[i, j] = Q[j][i]
        H = V.transpose()*(R*V)
        Hm = arb_mat(p, p)
        for i in range(p):
            for j in range(p): Hm[i, j] = ((H[i, j] + H[j, i])/2).mid()
        ritz = sorted([e.real.mid() for e in acb_mat(Hm).eig(algorithm="approx")], key=lambda z: float(z))
        lmin = ritz[0]
        n = nB; Lm = None
        tchol0 = time.time()
        if lmin > 0:
            for fac in ('0.5', '0.125', '0.01'):
                s0 = lmin*arb(fac); nu = s0/2
                Lm = chol_block(R, n, s0, BS) if (BS > 0 and n > BS) else chol_python(R, n, s0)
                if Lm is not None: break
        tchol = time.time() - tchol0
        if Lm is not None:
            del R
            LmT = Lm.transpose(); emax = arb(0); RB = 128
            for i0 in range(0, n, RB):        # residual IW - Lm Lm^T - s0 Id in row blocks (memory)
                i1 = min(n, i0 + RB)
                Lb = arb_mat(i1 - i0, n)
                for i in range(i0, i1):
                    for j in range(i+1): Lb[i - i0, j] = Lm[i, j]
                Pb = Lb*LmT
                for i in range(i0, i1):
                    rs = arb(0)
                    for j in range(n): rs += abs(IW[i, j] - Pb[i - i0, j] - (s0 if i == j else 0)).upper()
                    if rs.upper() > emax.upper(): emax = rs.upper()
                del Lb, Pb
            del Lm, LmT
            rig = bool(emax <= s0 - nu)      # then lambda_min(IW) >= s0 - emax >= nu > 0
        else: rig = False
    else:
        rig = True   # no negative part: P - N = P > 0 trivially (deff > 0, dl_j > 0 on G, w_k > 0 on K+, pole weight 2 > 0)
        tchol = 0
    res = dict(char='zeta', c=str(cq), sector=('defl' if defl else sector), NG=NG, M=os.environ.get('WS_M', '1') if defl else None,
               U=(Uenc.str(12, radius=True) if Uenc is not None else None), U_up=(arb(Uenc.upper()).str(12, radius=False) if Uenc is not None else None), eps=epsstr, prec=ctx.prec, K2=K2, t1=T1s, nb=[NBO, NBEK, NBEJ], kscan=ks, kscan_e=kse, sel_odd=len(selo), sel_even=len(sele),
               scan=sinfo.get('method', '?'), scanE=sinfo.get('E'), Kp=len(Kp), Km=len(Km), Jm=len(Jm), G=len(G), order=nB, nU=nU, atoms=len(S.atoms), sig2=float(S.sig2.mid()),
               margin=(lmin.str(12, radius=False) if lmin is not None else None), rigorous=rig, maxrad=maxrad.str(3, radius=False), ritz=[z.str(6, radius=False) for z in ritz[:3]], nu=(nu.str(6, radius=False) if nu is not None else None),
               resid=(emax.str(3, radius=False) if emax is not None else None),
               mu1=mu(1).str(8, radius=False), mu0=mu(0).str(8, radius=False), times=[round(tscan), round(tsym), round(tF), round(tasm), round(tW), round(tchol), round(time.time()-t0)])
    print(json.dumps(res), flush=True)

if __name__ == '__main__':
    main()
