# Finite-T numerics for Lamzouri 2611.14486 (Prop SmallSetZeros / Thm Main2 chain) on Odlyzko's zeros.
# usage: python3 num_pairs.py <zeros.gz> <comma-separated N cut-offs>
# Kernels (xi = (gamma-gamma') log T / 2pi, w(D) = 4/(4+D^2)):
#   f^(xi) = R(xi/lam)                 (Lamzouri's f = lam*R^(lam u), support [-1/lam,1/lam])
#   g0^(xi) = FT of f*1_[-th,th]       (sharp truncation; BGST Lemma 5 needs only L1, even, supp in [-1,1], Lipschitz at 0)
import sys, gzip, time, math
import numpy as np

LAM = 0.7892
TH = 4/7 - 1e-5
TWO_PI = 2*math.pi

def R(x):
    x = np.abs(x)
    den = 1.0 - x*x
    s = np.sinc(x)**2
    near = np.abs(den) < 1e-9
    return np.where(near, -(x - 1.0)/2.0, s/np.where(near, 1.0, den))

def fhat(xi):
    return R(xi/LAM)

def g0hat(xi):
    k = TWO_PI*np.abs(xi); m = TWO_PI*LAM; th = TH
    small = k*th < 1e-3
    ks = np.where(small, 1.0, k)
    I1 = np.where(small, th - k*k*th**3/6, np.sin(ks*th)/ks)
    I2 = np.where(small, th*th/2 - k*k*th**4/8, th*np.sin(ks*th)/ks + (np.cos(ks*th) - 1)/ks**2)
    def h(u):  # (1-cos(u th))/u
        sm = np.abs(u)*th < 1e-4
        us = np.where(sm, 1.0, u)
        return np.where(sm, u*th*th/2, (1 - np.cos(us*th))/us)
    I3 = 0.5*(h(m + k) + h(m - k))
    return 2*LAM*(I1 - LAM*I2 + I3/TWO_PI)

def J_of(th):  # Lamzouri's closed form J_th(f) = f(0) + 2 int_0^th f(a) a da
    l = LAM
    return l + l*th**2 - 2/3*l*l*th**3 - th*math.cos(2*math.pi*l*th)/(2*math.pi**2) + math.sin(2*math.pi*l*th)/(4*math.pi**3*l)

def selftest():
    # g0hat vs direct quadrature; fhat(0)=1; R vs FT of R^ at a few points
    from numpy.polynomial.legendre import leggauss
    xg, wg = leggauss(600)
    a = TH*(xg + 1)/2; wa = wg*TH/2
    fa = LAM*(1 - LAM*a + np.sin(TWO_PI*LAM*a)/TWO_PI)
    err = 0.0
    for xi in (0.0, 0.3, 0.7892, 1.0, 2.5, 7.3, 40.1):
        q = 2*np.sum(wa*fa*np.cos(TWO_PI*xi*a))
        err = max(err, abs(q - g0hat(np.array([xi]))[0]))
    u = (xg + 1)/2/LAM; wu = wg/(2*LAM)
    fu = LAM*(1 - LAM*u + np.sin(TWO_PI*LAM*u)/TWO_PI)
    err2 = 0.0
    for xi in (0.0, 0.5, 0.7892, 1.3, 3.7):
        q = 2*np.sum(wu*fu*np.cos(TWO_PI*xi*u))
        err2 = max(err2, abs(q - fhat(np.array([xi]))[0]))
    J_num = LAM + 2*np.sum(wa*fa*a)
    print(f'selftest: max|g0hat-quad|={err:.2e}  max|fhat-quad|={err2:.2e}  J_closed={J_of(TH):.12f}  J_quad={J_num:.12f}', flush=True)

def pair_sum(gam, L, kern, D):
    """S/N_set where S = sum_{i,j} kern(|g_i-g_j| L) w(g_i-g_j); offsets until all gaps > D."""
    n = len(gam); tot = n*kern(np.array([0.0]))[0]; d = 1
    while d < n:
        dl = gam[d:] - gam[:-d]
        msk = dl <= D
        if not msk.any(): break
        dd = dl[msk]
        for c in range(0, len(dd), 200_000):
            e = dd[c:c+200_000]
            tot += 2*np.sum(kern(e*L)*4/(4 + e*e))
        d += 1
    return tot, d

def tail(L, kern, D):  # continuum estimate of the omitted |Delta|>D part, per zero: 2L int_D^inf kern(DL) w dD
    x = np.linspace(D, D + 4000, 200_001)
    y = kern(x*L)*4/(4 + x*x)
    return 2*L*np.trapezoid(y, x)

def model(gam, L, kern, D):
    # GUE + local density heuristic: per zero at height t, rho=log(t/2pi)/2pi,
    # M(rho) = K(0) + 2 int_0^D K(Delta L) w(Delta) rho (1 - sinc^2(rho Delta)) dDelta ; averaged over zeros
    rho = np.log(gam/TWO_PI)/TWO_PI
    grid = np.linspace(max(rho.min(), 0.05), rho.max(), 120)
    x = np.linspace(0, D, 200_001); kx = kern(x*L)*4/(4 + x*x)
    M = np.array([kern(np.array([0.0]))[0] + 2*np.trapezoid(kx*r*(1 - np.sinc(r*x)**2), x) for r in grid])
    return np.mean(np.interp(rho, grid, M))

def load(path):  # streaming parse (low peak memory)
    from array import array
    a = array('d')
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line.strip(): a.append(float(line))
    return np.frombuffer(a, dtype=float).copy()

def main():
    path = sys.argv[1]; cuts = [int(c) for c in sys.argv[2].split(',')]
    gam_all = load(path)
    print(f'loaded {len(gam_all)} zeros, last {gam_all[-1]:.6f}', flush=True)
    selftest()
    for Ncut in cuts:
        t0 = time.time()
        gam = gam_all[:Ncut]; T = gam[-1] + 1e-9; N = len(gam); L = math.log(T)/TWO_PI
        h = 2*math.pi*LAM/math.log(T)
        gaps = np.diff(gam)
        frac_gap = np.sum(gaps <= h)/N            # D(lambda) estimate (n with gamma_{n+1} <= T)
        left = np.r_[np.inf, gaps]; right = np.r_[gaps, np.inf]
        inE = np.minimum(left, right) <= h; inE[0] = False   # Lamzouri's E_lambda(T): n>=2
        Z0 = gam[~inE]; delta = inE.mean()
        Sg, dg = pair_sum(gam, L, g0hat, 100.0); tg = tail(L, g0hat, 100.0)
        Sf, df = pair_sum(gam, L, fhat, 30.0)
        Sg0, _ = pair_sum(Z0, L, g0hat, 100.0)
        SgE, _ = pair_sum(gam[inE], L, g0hat, 100.0)
        cross = (Sg + SgE - Sg0)/2; cs = cross/math.sqrt(Sg*SgE) if SgE > 0 else float('nan')
        mod = model(gam, L, g0hat, 100.0)
        dscan = []
        if '--scan' in sys.argv:
            for lp in (0.05, 0.1, 0.2, 0.4):
                hp = 2*math.pi*lp/math.log(T); Ep = np.minimum(left, right) <= hp; Ep[0] = False
                Sp, _ = pair_sum(gam[~Ep], L, g0hat, 100.0)
                dscan.append((lp, Ep.mean(), Sp/N - Sg/N))
        Sf0, _ = pair_sum(Z0, L, fhat, 30.0)
        # Lemma DoubleSumZeros moments A_q/N, q=1..6
        s = np.ones(N); d = 1; LT = math.log(T)
        while d < N:
            dl = gam[d:] - gam[:-d]
            if dl.min() > 60: break
            k = 1/(1 + (dl*LT)**2); s[d:] += k; s[:-d] += k; d += 1
        Aq = [np.mean(s**q) for q in range(1, 7)]
        J = J_of(TH)
        print(f'T={T:.1f} N={N} logT={LT:.3f}  D(lam)_T={frac_gap:.4f}  delta=#E/N={delta:.4f}', flush=True)
        print(f'  J_theta(f)={J:.6f}  S_g0(Z)/N={Sg/N:.6f} (tail est {tg:+.1e}, offsets {dg})  diff={Sg/N + tg - J:+.5f}  1/sqrt(logT)={1/math.sqrt(LT):.3f} 1/logT={1/LT:.3f}', flush=True)
        print(f'  S_f(Z)/N={Sf/N:.6f}  S_g0(Z0)/N={Sg0/N:.6f}  S_f(Z0)/N={Sf0/N:.6f}  #Z0/N={len(Z0)/N:.6f}  checks: S_f(Z0)<=#Z0: {Sf0 <= len(Z0)+1e-6}, S_g0(Z0)<=S_f(Z0): {Sg0 <= Sf0 + 1e-3*N}', flush=True)
        print(f'  GUE+local-density model S_g0(Z)/N={mod:.6f} (obs-model={Sg/N-mod:+.5f}); predicted deficit (1+log2pi)(f(0)-2int a f)/logT={(1+math.log(TWO_PI))*(2*LAM-J)/LT:.4f}', flush=True)
        print(f'  S_g0(E)/N={SgE/N:.6f}  cross/sqrt(S(Z)S(E))={cs:.4f} (Cauchy-Schwarz needs <=1)', flush=True)
        for lp, dl_, dS in dscan:
            print(f'  scan lam\'={lp}: delta={dl_:.5f}  S_g0(Z-E)/N - S_g0(Z)/N={dS:+.5f}  ratio/delta={dS/dl_ if dl_>0 else 0:+.3f}', flush=True)
        print(f'  A_q/N q=1..6: ' + ' '.join(f'{a:.4f}' for a in Aq) + f'   [{time.time()-t0:.0f}s]', flush=True)

if __name__ == '__main__':
    main()
