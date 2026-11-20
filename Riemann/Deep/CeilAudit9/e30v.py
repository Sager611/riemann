# ceilaudit9 E30v: integer check of the LP adversary (e30.py) with explore30's exact two-bump kernel.
# f = phi_w(x - a') + e^{i theta} phi_w(x + a'), phi(s) = (1 - s^2)^8, a' = a - w, w = 0.03a, modulated to t0 (block centre);
# K(gamma, d) = Re[fhat(gamma - i d) conj(fhat(gamma + i d))], W_norm = sum_rho K / (D 2 pi ||f||^2).
# Block: M periods of the LP placement, integer zeros by cumulative rounding (k zeros at an atom: floor(k/2) off-line pairs d = 0.37 where
# cos(phase) < 0, the rest on-line); outside: on-line lattice of spacing 1/D. Reports W_norm, max|S|, max|S_1| over the block (S = N - D t).
import sys, math, numpy as np
sys.path.insert(0, 'Riemann/Deep/CeilAudit9')
from e30 import lp, bounds
xg, wg = np.polynomial.legendre.leggauss(96)
def phihat(z):   # int_{-1}^{1} (1-s^2)^8 e^{i z s} ds, z complex array
    return ((1 - xg**2)**8*wg*np.exp(1j*np.multiply.outer(z, xg))).sum(-1)
def run(a, T=3e12 + 1e4, wf=0.03, dmax=0.37, J=240, Lw=40.0, rm=1.0):
    D, Sm0, S1m0 = bounds(T); Sm, S1m = Sm0 - rm, S1m0 - rm     # LP keeps rm extra units for integer rounding; check against Sm0, S1m0
    w = wf*a; ap = a - w; P = math.pi/ap; n = D*P
    ch = math.cosh(2*ap*dmax); best = None
    for psi in np.arange(32)*math.pi/16:      # block starts at phase psi with S = S_1 = 0 (anchored LP)
        phi_ = 2*ap*(np.arange(J)*P/J) + psi
        u_ = np.where(np.cos(phi_) < 0, -ch*np.cos(phi_), -np.cos(phi_))
        r_ = lp(u_, n, D, P, Sm, S1m, J, anchor=True)
        if r_.status == 0 and (best is None or r_.fun < best[0].fun): best = (r_, phi_, psi)
    r, phi, psi = best; mu = r.x[:J]; Wlp = 1 - (-r.fun)/n
    ell = Lw/w; M = int(2*ell/P) + 1; tj = np.arange(J)*P/J
    mass = np.tile(mu, M); tt = (np.arange(M)[:, None]*P + tj[None, :]).ravel() - M*P/2; dT = P/J
    off = np.tile(np.cos(phi) < 0, M)
    # second-order sigma-delta rounding: integer counts k_i with N_i - C_i in (-1, 1) and running integral of (N - C) kept near 0
    nzm = np.nonzero(mass > 1e-9)[0]; k = np.zeros(len(mass), int); C = 0.0; Nn = 0; E = 0.0
    for ii, i in enumerate(nzm):
        C += mass[i]; gap = (tt[nzm[ii+1]] - tt[i]) if ii + 1 < len(nzm) else 0.0
        lo = int(np.floor(C - Nn)); best_ = None
        for kk in (lo, lo + 1):
            if kk < 0: continue
            e = Nn + kk - C
            if abs(e) >= 1: continue
            sc = abs(E + e*gap)
            if best_ is None or sc < best_[0]: best_ = (sc, kk, e)
        _, kk, e = best_; k[i] = kk; Nn += kk; E += e*gap
    gam = []; dl = []
    for t_, kk, o in zip(tt, k, off):
        if kk <= 0: continue
        npair = kk//2 if o else 0; non = kk - 2*npair
        # split the cluster symmetrically around its atom (distinct ordinates, centroid exactly t_): single on-line zero(s) at t_, pairs at +-offsets
        for q_ in range(non): gam += [t_]; dl += [0.0]
        offs = dT*(np.arange(npair) - (npair - 1)/2)/max(npair, 1)
        for o_ in offs: gam += [t_ + o_]*2; dl += [dmax]*2
    gam = np.array(gam); dl = np.array(dl)
    tb0, tb1 = tt[0], tt[-1] + P/J
    # S over the block relative to the lattice baseline, S_1 = int S
    gs = np.sort(gam); t0b = tt[0]                   # exact: S(t) = N(t) - D(t - t0b), S_1(t) = sum_{g<=t}(t - g) - D(t - t0b)^2/2
    Nk = np.arange(1, len(gs)+1); cs = np.cumsum(gs)
    S = np.concatenate([Nk - 1 - D*(gs - t0b), Nk - D*(gs - t0b)])             # just before / after each zero
    tc = t0b + Nk/D; ok = (tc > gs) & (np.concatenate([tc[:-1] < gs[1:], [True]]))  # S = 0 crossings between zeros
    S1z = Nk*gs - cs - D*(gs - t0b)**2/2; S1c = Nk*tc - cs - D*(tc - t0b)**2/2
    S1 = np.concatenate([S1z, S1c[ok]])
    # outside lattice
    h = 1/D; out = np.concatenate([np.arange(tb0 - h/2, tb0 - 2.5*ell, -h), np.arange(tb1 + h/2, tb1 + 2.5*ell, h)])
    allg = np.concatenate([gam, out]); alld = np.concatenate([dl, np.zeros(len(out))])
    th_ = -psi - M*math.pi                       # test-function phase matching the LP phase 2a't + psi at the block start
    def fh(z): return w*phihat(w*z)*(np.exp(1j*z*ap) + np.exp(1j*th_)*np.exp(-1j*z*ap))
    K = np.empty(len(allg))
    for i in range(0, len(allg), 4000):
        g_ = allg[i:i+4000]; d_ = alld[i:i+4000]
        K[i:i+4000] = (fh(g_ - 1j*d_)*np.conj(fh(g_ + 1j*d_))).real
    nf2 = 2*w*((1 - xg**2)**16*wg).sum()        # ||f||^2 = 2 w int phi^2
    Wn = K.sum()/(D*2*math.pi*nf2)
    print(f"a={a:.3f} (c={math.exp(2*a):.2f}) T={T:.0e}: anchored LP (psi={psi:.2f}) W_norm {Wlp:+.3f}; integer block: {len(gam)} zeros ({(dl>0).sum()} off-line), M={M} periods; "
          f"W_norm {Wn:+.4f}; max|S| {np.abs(S).max():.2f} (allowed {Sm0:.2f}); S_1 range [{S1.min():.2f}, {S1.max():.2f}] (allowed |S_1| <= {S1m0:.2f}); "
          f"{'ADMISSIBLE' if (np.abs(S).max() <= Sm0 and np.abs(S1).max() <= S1m0 and Wn < 0) else 'not admissible'}; max zeros at one ordinate {np.max(np.unique(gam, return_counts=True)[1])}", flush=True)
if __name__ == '__main__':
    T = float(sys.argv[2]) if len(sys.argv) > 2 else 3e12 + 1e4
    for a in [float(x) for x in sys.argv[1].split(',')]: run(a, T=T)
