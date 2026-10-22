# lib5.py -- rigorous paired-probe source floors (dbn5).
# Arb (python-flint) for every transcendental value; numpy for integer counting and for
# interval arithmetic with IEEE-754 outward rounding (nextafter after every operation).
import numpy as np, math, json, os, time
from fractions import Fraction as Fr
from flint import arb, fmpz, fmpq, ctx  # pyright: ignore[reportMissingImports]

ctx.prec = 96
def _repo():
    p = os.path.dirname(os.path.abspath(__file__))
    while p != os.path.dirname(p):
        if os.path.exists(os.path.join(p, 'lakefile.lean')) or os.path.exists(os.path.join(p, '.git')): return p
        p = os.path.dirname(p)
    raise SystemExit('repository root not found')
HERE = os.environ.get('DBN5_WORK', os.path.join(_repo(), '.scratch', 'dbn5'))
try:
    os.makedirs(HERE, exist_ok=True)
except OSError as _e:
    raise SystemExit('cannot create work dir %s: %s' % (HERE, _e))
N0 = 690988
X = 6000000185827
NG = 1920                       # global log grid: step 1/NG
ETA = fmpq(1, NG)
SUB = 4                         # sub-bin = 4 grid cells
BIN = 96                        # lambda bin = 96 grid cells  (w = 1/20)
PSI = 40                        # coarse cell for Psi_g
SPB = BIN // SUB                # sub-bins per bin
W = fmpq(BIN, NG)
DIVS = [1, 2, 3, 5, 6, 10, 15, 30]
MU = {1: 1, 2: -1, 3: -1, 5: -1, 6: 1, 10: 1, 15: 1, 30: -1}
PF = {1: [], 2: [2], 3: [3], 5: [5], 6: [2, 3], 10: [2, 5], 15: [3, 5], 30: [2, 3, 5]}
INF = np.inf
U53 = 2.0**-53

def up(x): return np.nextafter(x, INF)
def dn(x): return np.nextafter(x, -INF)
def fu(a): return float(np.nextafter(np.nextafter(float(a.upper()), INF), INF))
def fd(a): return float(np.nextafter(np.nextafter(float(a.lower()), -INF), -INF))
def arr_up(vals): return up(up(np.array([float(v.upper()) for v in vals], dtype=np.float64)))
def arr_dn(vals): return dn(dn(np.array([float(v.lower()) for v in vals], dtype=np.float64)))
def sum_up(x):
    x = np.asarray(x, dtype=np.float64).ravel()
    if x.size == 0: return 0.0
    assert np.all(x >= 0)
    s = float(np.sum(x))
    return float(up(up(s*(1.0 + 4.0*(x.size + 2)*U53))))
def aq(q):
    if isinstance(q, Fr): return arb(fmpq(q.numerator, q.denominator))
    return arb(q)

V0 = arb(N0).log()
LOGP = {p: arb(p).log() for p in (2, 3, 5)}
LOGD = {d: arb(d).log() for d in DIVS}
QD = {d: sum((LOGP[p]**2 for p in PF[d]), arb(0)) for d in DIVS}

# ---------------- global integer thresholds G[i] = ceil(exp(i/NG)) ----------------
def thresholds(imax):
    fn = os.path.join(HERE, 'G_thresholds.npy')
    if os.path.exists(fn):
        Gs = np.load(fn)
        if len(Gs) >= imax + 1: return Gs[:imax+1]
    out = np.zeros(imax+1, dtype=np.int64)
    old = ctx.prec; ctx.prec = 160
    e1 = arb(ETA).exp(); v = arb(1)
    for i in range(imax+1):
        if i % 4096 == 0: v = (arb(i)*arb(ETA)).exp()     # refresh to keep balls tight
        c = v.ceil().unique_fmpz()
        if c is None: raise RuntimeError('ambiguous threshold %d' % i)
        out[i] = int(c)
        v = v*e1
    ctx.prec = old
    assert np.all(np.diff(out) >= 0) and out[0] == 1
    np.save(fn, out)
    return out
SMAX_FIN = 1.6                  # largest finite subrange end
IMAX = int(math.ceil((math.log(30*N0 + 1) + SMAX_FIN + 0.02)*NG)) + 2*BIN + 8
G = thresholds(IMAX)
IV0 = int(math.ceil(math.log(N0)*NG)) + SUB + 4       # table length for u in [0, v0]
def gidx(k):
    return np.searchsorted(G, k, side='right') - 1     # log k in [i/NG, (i+1)/NG)
def coprime_upto(x):
    x = np.asarray(x, dtype=np.int64)
    return sum(MU[d]*(x//d) for d in DIVS)
SMOOTH = sorted(set(2**a*3**b*5**c for a in range(21) for b in range(13) for c in range(9) if 2**a*3**b*5**c <= N0))
def rad(g): return math.prod(p for p in (2, 3, 5) if g % p == 0)

# ---------------- global label index data (independent of the cell) ----------------
JMAX = IMAX // SUB - 2
def build_index():
    fn = os.path.join(HERE, 'index_data.npz')
    if os.path.exists(fn):
        z = np.load(fn); return {k: z[k] for k in z.files}
    jj = np.arange(JMAX+1, dtype=np.int64)
    mlo = G[SUB*jj]; mhi = G[SUB*(jj+1)] - 1
    H = {k: [] for k in ('g', 'cnt', 'ilo', 'ihi', 'bin')}
    for g in SMOOTH:
        kmax = N0 // g
        klo = (mlo + g - 1)//g; khi = np.minimum(mhi//g, kmax)
        if g == 1: klo = np.maximum(klo, 2)
        ok = klo <= khi
        cnt = np.zeros_like(klo)
        cnt[ok] = coprime_upto(khi[ok]) - coprime_upto(klo[ok]-1)
        ok &= cnt > 0
        H['g'].append(np.full(ok.sum(), g)); H['cnt'].append(cnt[ok])
        H['ilo'].append(gidx(klo[ok])); H['ihi'].append(gidx(khi[ok]) + 1); H['bin'].append(jj[ok]//SPB)
    T = {k: [] for k in ('d', 'cnt', 'ilo', 'ihi', 'bin')}
    for d in DIVS[1:]:
        nlo = np.maximum((mlo + d - 1)//d, N0//d + 1); nhi = np.minimum(mhi//d, N0)
        ok = nlo <= nhi
        T['d'].append(np.full(ok.sum(), d)); T['cnt'].append((nhi - nlo + 1)[ok])
        T['ilo'].append(gidx(nlo[ok])); T['ihi'].append(gidx(nhi[ok]) + 1); T['bin'].append(jj[ok]//SPB)
    R = {k: [] for k in ('d', 'cnt', 'j', 'kr')}
    jr = np.arange(0, IV0//SUB, dtype=np.int64)
    nlo = G[SUB*jr]; nhi = np.minimum(G[SUB*(jr+1)] - 1, N0)
    ok = nlo <= nhi
    for d in DIVS:
        kr = []
        for j in jr[ok]:
            if d == 1:
                kr.append(-((SUB*int(j) + SUB)*20 + NG - 1)//NG)       # floor(-(4j+4)/96)
            else:
                v = (LOGD[d] - arb(fmpq(SUB*int(j) + SUB, NG)))*20
                f = v.floor().unique_fmpz()
                if f is None: raise RuntimeError('ambiguous kr')
                kr.append(int(f))
        R['d'].append(np.full(ok.sum(), d)); R['cnt'].append((nhi - nlo + 1)[ok]); R['j'].append(jr[ok]); R['kr'].append(np.array(kr, dtype=np.int64))
    out = {}
    for name, D in (('H', H), ('T', T), ('R', R)):
        for k, v in D.items(): out[name + '_' + k] = np.concatenate(v).astype(np.int64)
    np.savez(fn, **out)
    return out
IDX = build_index()
def _kr_table():
    fn = os.path.join(HERE, 'kr_table.npy')
    jmax = IMAX//SUB - 2
    if os.path.exists(fn):
        t = np.load(fn)
        if t.shape[1] >= jmax + 1: return t
    out = np.zeros((len(DIVS), jmax + 1), dtype=np.int64)
    for di, d in enumerate(DIVS):
        for j in range(jmax + 1):
            if d == 1:
                out[di, j] = -((SUB*j + SUB)*20 + NG - 1)//NG
            else:
                v = (LOGD[d] - arb(fmpq(SUB*j + SUB, NG)))*20
                f = v.floor().unique_fmpz()
                if f is None: raise RuntimeError('ambiguous kr')
                out[di, j] = int(f)
    np.save(fn, out)
    return out
KRT = _kr_table()
KF = (IMAX//SUB - 2)//SPB + 2
KRMIN = int(KRT.min()) - 1; KRMAX = int(KRT.max()) + 1

# ---------------- per (cell, height, frequency subrange): rigorous bin masses ----------------
def _recur_table(a2, sig, n):
    """upper bounds of exp(a2*u^2 - sig*u) at u = i/NG, i = 0..n-1 (arb recurrence, refreshed)."""
    eta = arb(ETA)
    al = (a2*eta*eta - sig*eta).exp(); be = (2*a2*eta*eta).exp()
    out = []; E = arb(1); B = arb(1)
    for i in range(n):
        if i % 2048 == 0 and i > 0:
            u = arb(i)*eta; E = (a2*u*u - sig*u).exp(); B = (2*a2*eta*eta*i).exp()
        out.append(E)
        E = E*al*B; B = B*be
    return arr_up(out)

class Rect:
    def __init__(self, h, cell, s1, s2):
        self.h = Fr(h); self.cell = cell; self.s1 = Fr(s1); self.s2 = None if s2 is None else Fr(s2)
        self.tL = Fr(cell, 4000); self.tR = Fr(cell + 1, 4000)
        h_ = aq(self.h); tL = aq(self.tL); tR = aq(self.tR); eps = arb(fmpq(1, 10**12))
        self.L1 = 2*(V0 + aq(self.s1)) - eps
        self.L2 = None if s2 is None else 2*(V0 + aq(self.s2)) + eps
        self.linked_refl = self.L2 is not None
        L1 = self.L1
        sigU = (1 + h_)/2 + tL*L1/4; sigV = (1 - h_)/2 + tL*L1/4
        if self.L2 is not None:
            nm = (self.L2/2).exp().floor().unique_fmpz()
            self.Nmax = (int(nm) if nm is not None else int((self.L2/2).exp().upper().floor().unique_fmpz()) + 1) + 1
        else:
            self.Nmax = N0
        self.IV = int(gidx(self.Nmax)) + 2*SUB + 4
        self.EU = _recur_table(tL/4, sigU, self.IV + 2)
        self.EV = _recur_table(tL/4, sigV, self.IV + 2)
        # family constants
        self.Ag = {g: fu((-sigU*arb(g).log()).exp()) for g in SMOOTH}
        Td = {d: (tL*QD[d]/4 - sigU*LOGD[d]).exp() for d in DIVS}
        self.Td = {d: fu(Td[d]) for d in DIVS}
        self.Rd = {d: fu(Td[d]*(-h_*L1/2).exp()) for d in DIVS}
        # Euler factor data (rho_max carries Raymond's e^{1e-5} allowance)
        self.rmax = {}; self.rmin = {}
        for p in (2, 3, 5):
            lp = LOGP[p]
            self.rmax[p] = (tL*(lp*lp - L1*lp)/4 - (1 + h_)/2*lp + arb(fmpq(1, 10**5))).exp()
            self.rmin[p] = arb(0) if self.L2 is None else (tR*(lp*lp - self.L2*lp)/4 - (1 + h_)/2*lp).exp()
        En = arb(1)
        for p in (2, 3, 5): En = En*(1 + self.rmax[p])
        Ey = arb(0)
        for p in (2, 3, 5):
            pr = arb(1)
            for q in (2, 3, 5):
                if q != p: pr = pr*(1 + self.rmax[q])
            Ey += LOGP[p]*self.rmax[p]*pr/2
        self.En = En; self.Ey = Ey
        if self.L2 is None:
            Aa, Da = self.activation()
            self.Aact = En*Aa; self.Bact = 2*(En*Da + Ey*Aa)
        else:
            self.Aact = arb(0); self.Bact = arb(0)       # activated labels are explicit linked labels (masses())

    def activation(self):
        """Raymond, Analytic Foundations (8.1)-(8.3): envelopes for every activated label (N0 < n <= N)."""
        h_ = aq(self.h); tL = aq(self.tL); v0 = V0; dl = arb(fmpq(1, 10**12))
        a0 = tL*v0/2; b = (1 - h_)/2; bp = (1 + h_)/2; c0 = a0 + h_; k = a0 - b
        C0 = (-tL*v0*v0/4 + b*v0).exp()*arb(fmpq(1, 10**5)).exp()
        def J0(be, z):
            if be == 0: return z
            return ((be*z).exp() - 1)/be
        def J1(be, z):
            if be == 0: return z*z/2
            return z*(be*z).exp()/be - ((be*z).exp() - 1)/(be*be)
        bq = Fr(1) - self.h; bq = bq/2; bpq = (Fr(1) + self.h)/2
        bA = arb(0) if bq == 0 else b
        l = aq(self.s1)
        if self.s2 is not None:
            u = aq(self.s2) + dl
            A_ = (-a0*l).exp()*J0(bA if bq != 0 else 0, u) + (-c0*l).exp()*J0(bp, u)
            D_ = (-a0*l).exp()*(v0*J0(bA if bq != 0 else 0, u) + J1(bA if bq != 0 else 0, u))/2 \
                 + (-c0*l).exp()*((v0/2 + u)*J0(bp, u) - J1(bp, u)/2)
        else:
            assert k > 0
            z = l + dl
            A_ = (-a0*l).exp()*J0(bA if bq != 0 else 0, z) + (-k*l).exp()/k + (-c0*l).exp()*J0(bp, z) + (-k*l).exp()/k
            D_ = (-a0*l).exp()*(v0*J0(bA if bq != 0 else 0, z) + J1(bA if bq != 0 else 0, z))/2 \
                 + (-k*l).exp()*((v0 + l)/(2*k) + 1/(2*k*k)) \
                 + (-c0*l).exp()*((v0/2 + l + dl)*J0(bp, z) - J1(bp, z)/2) \
                 + (-k*l).exp()*((v0 + l + 2*dl)/(2*k) + 1/(2*k*k))
        return C0*A_, C0*D_

    def masses(self, Psi):
        """Psi: per-cell table object (PsiTable). Returns dict with bin data (all floats are rigorous bounds)."""
        I = IDX
        EU = self.EU; EV = self.EV
        # head families
        g = I['H_g']; cnt = I['H_cnt'].astype(np.float64); ilo = I['H_ilo']; ihi = I['H_ihi']
        Ag = np.array([self.Ag[x] for x in SMOOTH])
        gi = np.searchsorted(np.array(SMOOTH), g)
        Em = np.maximum(EU[ilo], EU[ihi])
        amp = Ag[gi]*Em
        big = g > 1
        ca = ilo[big]//PSI; cb = -((-ihi[big])//PSI)
        psib = up(Psi.Pp[gi[big], cb] - Psi.Pn[gi[big], ca])
        amp[big] = amp[big]*psib*Psi.Dt[cb]
        amp = up(amp*(1 + 8*U53))
        mH = up(amp*cnt*(1 + 4*U53))
        # tail families
        d = I['T_d']; Td = np.array([self.Td[x] for x in DIVS])
        di = np.searchsorted(np.array(DIVS), d)
        ampT = up(Td[di]*np.maximum(EU[I['T_ilo']], EU[I['T_ihi']])*(1 + 4*U53))
        mT = up(ampT*I['T_cnt'].astype(np.float64)*(1 + 4*U53))
        # reflected
        dr = I['R_d']; Rd = np.array([self.Rd[x] for x in DIVS]); dri = np.searchsorted(np.array(DIVS), dr)
        j = I['R_j']
        ampR = up(Rd[dri]*np.maximum(EV[SUB*j], EV[SUB*j + SUB])*(1 + 4*U53))
        mR = up(ampR*I['R_cnt'].astype(np.float64)*(1 + 4*U53))
        # activated labels N0 < n <= Nmax (finite subranges): first branch (d,n) for all d|30, reflected (d,n)
        actF_bin = []; actF_m = []; actF_a = []; actR_bin = []; actR_m = []; actR_a = []
        self.bd_amp = 0.0; self.bd_lam = 0.0
        if self.L2 is not None and self.Nmax > N0:
            Nm = self.Nmax
            jlo = int(gidx(N0 + 1))//SUB; jhi = int(gidx(30*Nm))//SUB + 1
            jj = np.arange(jlo, jhi + 1, dtype=np.int64)
            mlo = G[SUB*jj]; mhi = G[SUB*(jj + 1)] - 1
            for d in DIVS:
                nlo = np.maximum((mlo + d - 1)//d, N0 + 1); nhi = np.minimum(mhi//d, Nm)
                ok = nlo <= nhi
                if not ok.any(): continue
                a_ = up(self.Td[d]*np.maximum(EU[gidx(nlo[ok])], EU[gidx(nhi[ok]) + 1])*(1 + 4*U53))
                actF_bin.append(jj[ok]//SPB); actF_a.append(a_); actF_m.append(up(a_*(nhi - nlo + 1)[ok].astype(np.float64)*(1 + 4*U53)))
            j2lo = int(gidx(N0 + 1))//SUB; j2hi = int(gidx(Nm))//SUB
            jr = np.arange(j2lo, j2hi + 1, dtype=np.int64)
            nlo = np.maximum(G[SUB*jr], N0 + 1); nhi = np.minimum(G[SUB*(jr + 1)] - 1, Nm); ok = nlo <= nhi
            for di, d in enumerate(DIVS):
                a_ = up(self.Rd[d]*np.maximum(EV[SUB*jr[ok]], EV[SUB*jr[ok] + SUB])*(1 + 4*U53))
                actR_bin.append(KRT[di, jr[ok]]); actR_a.append(a_); actR_m.append(up(a_*(nhi - nlo + 1)[ok].astype(np.float64)*(1 + 4*U53)))
            # boundary labels (present at one probe only): one n* per family, amplitude <= max single amplitude
            iA = int(gidx(N0 + 1)); iB = int(gidx(Nm)) + 1
            eU = float(np.max(EU[iA:iB + 1])); eV = float(np.max(EV[iA:iB + 1]))
            self.bd_amp = float(up(sum(up(self.Td[d]*eU) for d in DIVS) + sum(up(self.Rd[d]*eV) for d in DIVS)))
        mAF = np.concatenate(actF_m) if actF_m else np.zeros(0); bAF = np.concatenate(actF_bin) if actF_bin else np.zeros(0, dtype=np.int64)
        aAF = np.concatenate(actF_a) if actF_a else np.zeros(0)
        mAR = np.concatenate(actR_m) if actR_m else np.zeros(0); bAR = np.concatenate(actR_bin) if actR_bin else np.zeros(0, dtype=np.int64)
        aAR = np.concatenate(actR_a) if actR_a else np.zeros(0)
        ntot = mH.size + mT.size + mR.size + mAF.size + mAR.size
        fac = 1.0 + 4.0*(ntot + 2)*U53
        MF = up(up((np.bincount(I['H_bin'], mH, minlength=KF) + np.bincount(I['T_bin'], mT, minlength=KF) + np.bincount(bAF, mAF, minlength=KF))*fac))
        aF = np.zeros(KF); np.maximum.at(aF, I['H_bin'], amp); np.maximum.at(aF, I['T_bin'], ampT); np.maximum.at(aF, bAF, aAF)
        nR = KRMAX - KRMIN + 1
        MR = up(up((np.bincount(I['R_kr'] - KRMIN, mR, minlength=nR) + np.bincount(bAR - KRMIN, mAR, minlength=nR))*fac))
        aR = np.zeros(nR); np.maximum.at(aR, I['R_kr'] - KRMIN, ampR); np.maximum.at(aR, bAR - KRMIN, aAR)
        # lambda intervals
        kF = np.arange(KF)
        lamF_lo = kF/20.0; lamF_hi = (kF + 1)/20.0                  # exact in binary? k/20 not exact -> round outward
        lamF_lo = dn(lamF_lo); lamF_hi = up(lamF_hi)
        kR = np.arange(KRMIN, KRMAX + 1)
        L1lo = fd(self.L1)
        om = (self.L2 - self.L1) if self.L2 is not None else arb(0)
        om = om + arb(fmpq(SUB, NG))
        L1hi_om = fu(self.L1 + om)
        lamR_lo = dn(dn(L1lo + dn(kR/20.0)))
        lamR_hi = up(up(L1hi_om + up((kR + 1)/20.0)))
        return dict(bd_amp=self.bd_amp, MF=MF, aF=aF, lamF_lo=lamF_lo, lamF_hi=lamF_hi, kF=kF,
                    MR=MR, aR=aR, lamR_lo=lamR_lo, lamR_hi=lamR_hi, kR=kR, omega=fu(om), L1lo=L1lo, L1hi=fu(self.L1))

class PsiTable:
    """Per cell (depends on t_L, t_R only): upper/lower tables for Psi_g(t_R, u) on the coarse grid u_c = c*PSI/NG."""
    def __init__(self, cell):
        tL = aq(Fr(cell, 4000)); tR = aq(Fr(cell + 1, 4000))
        CM = IV0//PSI + 3
        sm = SMOOTH; ns = len(sm)
        step = arb(fmpq(PSI, NG))
        # Y_r(c) = exp(tR/2 * c*step * log r)
        Y_up = np.zeros((ns, CM)); Y_dn = np.zeros((ns, CM))
        for ri, r in enumerate(sm):
            q = (tR/2*step*arb(r).log()).exp(); y = arb(1); vals = []
            for c in range(CM):
                if c % 256 == 0 and c > 0: y = (tR/2*step*c*arb(r).log()).exp()
                vals.append(y); y = y*q
            Y_up[ri] = arr_up(vals); Y_dn[ri] = arr_dn(vals)
        self.Pp = np.zeros((ns, CM)); self.Pn = np.zeros((ns, CM))
        pos = {x: i for i, x in enumerate(sm)}
        for gi, g in enumerate(sm):
            if g == 1: continue
            R = rad(g)
            for d in DIVS:
                if R % d: continue
                r = g//d
                K = (tR/4*(QD[d] + arb(r).log()**2)).exp()
                if MU[d] > 0: self.Pp[gi] = up(self.Pp[gi] + up(fu(K)*Y_up[pos[r]]))
                else: self.Pn[gi] = dn(self.Pn[gi] + dn(fd(K)*Y_dn[pos[r]]))
        cc = np.arange(CM)
        self.Dt = arr_up([((tR - tL)*(step*c)**2/4).exp() for c in range(CM)])
