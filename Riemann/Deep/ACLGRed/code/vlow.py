# vlow.py: QW_lambda on the span of ACLG's conditioning vectors 1_W(e^lam_m * e_n^*) attached to deep modes (delta_m(lam*sigma) <= tau),
# archimedean pair (S = {inf}): e_n = sqrt2 u^{1/2} h_{n,sigma}(u), e^lam_m = sqrt2 u^{1/2} h_{m,lam sigma}(u); for x <= lam,
# (e^lam_m * e_n^*)(x) = 2 x^{1/2} int_0^sigma h_{m,lam sigma}(x w) h_{n,sigma}(w) dw. QW_lambda: Galerkin matrix of WeilSimple/weil_matrix.py.
# Usage: vlow.py c sigma tau nlist N prec    (c = lambda^2; nlist like 0 or 0,2,4)
import sys, json, math, time, os
import numpy as np
sys.path.insert(0, '/Users/adriansager/playground/riemann/Riemann/Deep/WeilSimple')
from flint import arb, arb_mat, ctx, fmpq
from weil_matrix import WeilMat

def legendre_norm_even(z, K):        # rows: even k = 0,2,..,K of sqrt((2k+1)/2) P_k(z)
    P = np.zeros((K + 1, z.size)); P[0] = 1.0
    if K >= 1: P[1] = z
    for k in range(1, K): P[k + 1] = ((2*k + 1)*z*P[k] - k*P[k - 1])/(k + 1)
    ks = np.arange(0, K + 1, 2)
    return ks, P[ks]*np.sqrt((2*ks + 1)/2.0)[:, None]

class Prolates:                       # even prolates at scale t (t^2 = ct): h_{n,t}(x) = t^{-1/2} Psi_n(x/t), n = 0,2,4,...
    def __init__(s, ct, nmax):
        s.t = math.sqrt(ct); g = 2*math.pi*ct; g2 = g*g
        K = int(2*g + 80); K += K % 2; s.K = K
        ks = np.arange(0, K + 1, 2).astype(float)
        dg = ks*(ks + 1) + g2*(2*ks*ks + 2*ks - 1)/((2*ks - 1)*(2*ks + 3))
        kk = ks[:-1]; of = g2*(kk + 1)*(kk + 2)/((2*kk + 3)*np.sqrt((2*kk + 1)*(2*kk + 5)))
        w, V = np.linalg.eigh(np.diag(dg) + np.diag(of, 1) + np.diag(of, -1))
        s.V = V[:, :nmax//2 + 1]
        _, P0 = legendre_norm_even(np.array([0.0]), K)
        phi0 = s.V.T @ P0[:, 0]
        chi = s.t*math.sqrt(2)*s.V[0, :]/phi0
        s.delta = 1 - np.abs(chi); s.sign = np.sign(chi)
        s.V = s.V*np.sign(phi0)[None, :]          # h_n(0) > 0
    def h(s, idx, x):                 # x array, |x| <= t
        _, P = legendre_norm_even(np.clip(x/s.t, -1, 1), s.K)
        return (s.V[:, idx].T @ P)/math.sqrt(s.t)

def main():
    c = float(sys.argv[1]); sig = float(sys.argv[2]); tau = float(sys.argv[3])
    nlist = [int(v) for v in sys.argv[4].split(',')]; N = int(sys.argv[5]); ctx.prec = int(sys.argv[6])
    lam = math.sqrt(c); ell = math.log(lam); L = 2*ell
    t0 = time.time()
    big = Prolates(c*sig*sig, int(4*c*sig*sig) + 40); small = Prolates(sig*sig, max(nlist) + 2)
    mlow = [i for i in range(big.V.shape[1]) if big.delta[i] <= tau]
    Ny, Nw = 600, 600
    yy, wy = np.polynomial.legendre.leggauss(Ny); y = ell*yy; wy = ell*wy; x = np.exp(y)
    ww, www = np.polynomial.legendre.leggauss(Nw); w = 0.5*sig*(ww + 1); www = 0.5*sig*www
    hn = {n: small.h(n//2, w) for n in nlist}
    Vsel = big.V[:, mlow]; Hn = np.array([www*hn[n] for n in nlist]).T      # Nw x len(nlist)
    Cf3 = np.empty((Ny, len(nlist), len(mlow)))
    for a in range(Ny):                                     # row by row: P is (K+1) x Nw (memory cap)
        _, P = legendre_norm_even(np.clip(x[a]*w/big.t, -1, 1), big.K)
        Hrow = (Vsel.T @ P)/math.sqrt(big.t)                # M_low x Nw
        Cf3[a] = 2*math.sqrt(x[a])*(Hrow @ Hn).T            # len(nlist) x M_low
    Cf = Cf3.reshape(Ny, -1)                                # Ny x (len(nlist)*M_low)
    # coefficients on the WeilSimple basis
    js = np.arange(0, N + 1); wj = 2*math.pi*js/L
    PhiE = np.sqrt(2/L)*np.cos(np.outer(y, wj)); PhiE[:, 0] = 1/math.sqrt(L)
    PhiO = np.sqrt(2/L)*np.sin(np.outer(y, wj[1:]))
    CE = PhiE.T @ (wy[:, None]*Cf); CO = PhiO.T @ (wy[:, None]*Cf)
    norm_true = np.sqrt((wy[:, None]*Cf*Cf).sum(0)); norm_N = np.sqrt((CE*CE).sum(0) + (CO*CO).sum(0))
    W = WeilMat(fmpq(int(round(c))) if abs(c - round(c)) < 1e-12 else fmpq(int(round(c*1000)), 1000))
    def tonp(T): return np.array([[float(T[i, j].mid().str(25, radius=False)) for j in range(T.ncols())] for i in range(T.nrows())])
    QE = tonp(W.block(N, 'even')[0]); QO = tonp(W.block(N, 'odd')[0])
    gal = np.sort(np.concatenate([np.linalg.eigvalsh(QE), np.linalg.eigvalsh(QO)]))
    # orthonormal basis of the numerically significant part of span{P_N c}
    Cn_E = CE/norm_N; Cn_O = CO/norm_N
    G = Cn_E.T @ Cn_E + Cn_O.T @ Cn_O
    s, U = np.linalg.eigh(G); s = s[::-1]; U = U[:, ::-1]
    out = {'c': c, 'sigma': sig, 'tau': tau, 'nlist': nlist, 'N': N, 'M_low': len(mlow), 'm_max': 2*max(mlow) if mlow else None,
           'nvec': Cf.shape[1], 'proj_loss_max': float(np.max(1 - (norm_N/norm_true)**2)),
           'gram_eigs': [float(v) for v in s], 'gram_eigs_rank': {str(e): int((s > e*s[0]).sum()) for e in (1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12)},
           'galerkin_QW_low': [float(v) for v in gal[:12]], 'N_QW_below_1e-2': int((gal < 1e-2).sum())}
    for cut in (1e-4, 1e-8, 1e-11):
        r = int((s > cut*s[0]).sum()); B = U[:, :r]/np.sqrt(s[:r])[None, :]
        PE = Cn_E @ B; PO = Cn_O @ B
        Qc = PE.T @ QE @ PE + PO.T @ QO @ PO
        ev = np.sort(np.linalg.eigvalsh(0.5*(Qc + Qc.T)))
        out['span_cut_%g' % cut] = {'dim': r, 'QW_eigs': [float(v) for v in ev], 'n_below_1e-2': int((ev < 1e-2).sum()),
                                    'n_ge_0.05': int((ev >= 0.05).sum())}
    out['sec'] = round(time.time() - t0, 1)
    print(json.dumps(out), flush=True)

main()
