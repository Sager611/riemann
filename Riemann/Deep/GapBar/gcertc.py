"""gapbar: Arb certificate that BM's criterion (general r) fails at kappa for EVERY COMPLEX amplifier (any degree, the L2
closure, any v, any theta <= 1/4).  Extends GapCeil gcert_r.py (real P, Re Gam only) to the Hermitian form:
  Q_L(phi) = sum_ij Gam_ij conj(I^i phi)^T W_{b+p} (I^j phi),  Gam = ReH + i ImH (ReH symmetric, ImH antisymmetric).
For phi = ar + i ai in the orthonormal Jacobi basis q_k of L2(X^b): M = A_re + i A_im, and
  Q_L = [ar; ai]^T [[A_re, -A_im], [A_im, A_re]] [ar; ai]  (real 2N embedding; orientation R = conjugate, same spectrum).
High modes (w in W_N, complex): as gcert_r with |Gam_ij| <= |ReH_ij| + |ImH_ij|, except (0,1)+(1,0) for p >= 1:
  2 Re(Gam_01 <w, X^(b+p) I w>) >= -|ReH_01| (b+p) th^(p-1) t_1^2 - 2 |ImH_01| th^p t_1  (by parts only for the real part).
Schur step with the doubled window E2 = [[E_re, -E_im], [E_im, E_re]]; tail eps as in gcert_r (majorant bounds |Gam| itself)."""
import time, math
from flint import arb, arb_mat, fmpq
import gcert as GC
import gcert_r as GR
import certc as CC


class OpsC(GR.OpsR):
    def __init__(self, kap, N, D, gam, r, Nw=None):
        # build everything as OpsR but with the Hermitian Gamma arrays
        orig = GR.CR.gamma_arrays_r
        GR.CR.gamma_arrays_r = CC.gamma_arrays_r
        try:
            super().__init__(kap, N, D, gam, r, Nw)
        finally:
            GR.CR.gamma_arrays_r = orig
        self.GamC = self.Gam                       # {'a': (ReH, ImH), ...}
        self.cacheC = {}

    def AC(self, h, p, part):
        key = (h, p, part)
        if key in self.cacheC: return self.cacheC[key]
        G = self.GamC[h][part]; D = self.D; W = self.Wq[p]
        H = None
        for i in range(D - 1, -1, -1):
            Z = arb_mat(self.n2, self.N)
            for j in range(D):
                Z = Z + G[i][j] * self.Y[j]
            WZ = W * Z
            H = WZ if H is None else WZ + self.IsT * H
        AP = GC.sub(H, self.Nw, self.N)
        out = self.QT * AP * self.QTN.transpose()
        self.cacheC[key] = out
        return out


def C_eff_c(Gre, Gim, p, D, gam, r, cfs, N, b):
    th = GC.th
    key = ('t', N, b)
    if key not in GR._TC:
        GR._TC[key] = GR.tvec_b(N, b, 345)
    t = GR._TC[key]
    a = lambda i, j: abs(Gre[i][j]) + abs(Gim[i][j])
    s = arb(0)
    if p == 0:
        s += 2 * a(0, 1) * t[1]
    else:
        s += abs(Gre[0][1]) * (b + p) * th ** (p - 1) * t[1] ** 2 + 2 * abs(Gim[0][1]) * th ** p * t[1]
    for j in range(2, D):
        s += 2 * a(0, j) * th ** p * t[j]
    for i in range(1, D):
        for j in range(1, D):
            s += a(i, j) * th ** p * t[i] * t[j]
    tk = ('tail', N, D, gam, r, b)
    if tk not in GR._TC:
        tt = arb(0)
        for i in range(340):
            for j in range(340):
                if i >= D or j >= D:
                    w = 2 * t[max(i, j)] if (i == 0 or j == 0) else t[i] * t[j]
                    tt += GR.tail_gam(i, j, gam, r) * w
        GR._TC[tk] = tt + arb(10) ** -60
    return s + GR._TC[tk] * cfs


def certify_c(ops, terms, label=''):
    N, Nw, D, gam, r, b = ops.N, ops.Nw, ops.D, ops.gam, ops.r, ops.b
    Mr = arb_mat(Nw, N); Mi = arb_mat(Nw, N)
    for (h, p, cf) in terms:
        Mr = Mr + cf * ops.AC(h, p, 0)
        Mi = Mi + cf * ops.AC(h, p, 1)
    Br = Mr.tolist(); Bi = Mi.tolist()
    eps = sum((abs(cf) * GR.tail_norm_r(D, gam, r, p) for (h, p, cf) in terms), arb(0))
    dD = None; nsub = 256
    for tt in range(nsub):
        lo = GC.A_(GC.TH) * tt / nsub; hi = GC.A_(GC.TH) * (tt + 1) / nsub
        Xb = (lo + hi) / 2 + arb(0, ((hi - lo) / 2).upper())
        dX = sum((cf * ops.GamC[h][0][0][0] * Xb ** p for (h, p, cf) in terms), arb(0))
        l = dX.lower()
        dD = l if dD is None or l < dD else dD
    dD = arb(dD)
    byp = {}
    for (h, p, cf) in terms:
        Gr, Gi = ops.GamC[h]
        if p not in byp:
            byp[p] = ([[cf * Gr[m][n] for n in range(D)] for m in range(D)], [[cf * Gi[m][n] for n in range(D)] for m in range(D)])
        else:
            R0, I0 = byp[p]
            byp[p] = ([[R0[m][n] + cf * Gr[m][n] for n in range(D)] for m in range(D)], [[I0[m][n] + cf * Gi[m][n] for n in range(D)] for m in range(D)])
    C = arb(0)
    for p, (Gr, Gi) in byp.items():
        cfs = sum((abs(cf) for (h, pp, cf) in terms if pp == p), arb(0))
        C += C_eff_c(Gr, Gi, p, D, gam, r, cfs, N, b)
    dW = dD - C
    info = {'label': label, 'deltaD': dD.str(5), 'C': C.str(5), 'deltaW': dW.str(5), 'eps_tail': eps.str(3)}
    if not dW > 0:
        info['fail'] = 'deltaW <= 0'
        return False, info
    n = 2 * N
    def big(k, l):     # doubled real matrix entry of the Nw x N complex block: rows k < 2Nw, cols l < 2N
        kr, kb = k % Nw, k // Nw; lr, lb = l % N, l // N
        if kb == lb: return Br[kr][lr]
        return -Bi[kr][lr] if kb == 0 else Bi[kr][lr]
    # V block, rows/cols ordered (re 0..N-1, im 0..N-1)
    S = [[None] * n for _ in range(n)]
    for k in range(n):
        for l in range(n):
            kk = (k % N) + (k // N) * Nw
            S[k][l] = big(kk, l) + arb(0, eps.upper())
    # window rows k in [N, Nw) for both parts
    Wrows = [k for k in range(N, Nw)] + [Nw + k for k in range(N, Nw)]
    E = [[big(w, l) for l in range(n)] for w in Wrows]
    Em = arb_mat(E)
    EtE = (Em.transpose() * Em).tolist()
    for k in range(n):
        for l in range(n):
            S[k][l] -= (1 + GC.ETA) * EtE[k][l] / dW
    for k in range(n):
        S[k][k] -= (1 + 1 / GC.ETA) * eps ** 2 / dW
    ok, kk, d = GC.cholesky_ok(arb_mat(S))
    info['cholesky'] = 'OK' if ok else ('FAIL at %d (%s)' % (kk, d.str(5)))
    return ok, info
