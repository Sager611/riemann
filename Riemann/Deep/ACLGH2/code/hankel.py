# hankel.py: the semilocal pair (P_t, hat P_t) of ACLG 2611.13129 Sec. 8 (l. 1097-1262) on a Galerkin (Laguerre) truncation, and the
# localisation quantity of (H2ii) (eq. H2ii, l. 1257):  T(n0, R) = Tr(M_{chi_R} N_t^(n0)) = sum_{n > n0} kappa_n(t) ||chi_R^{1/2} a_n||^2.
# Model (derivation in ../DERIVATION.md): Mellin variable s = i tau; Cayley z = (tau - i)/(tau + i) maps L^2(R, dtau/2pi) unitarily to
# L^2(S^1) with ONB F_k ~ z^k/(tau + i); functions supported in (0,1] <-> span{F_k, k <= -1}, in [1, inf) <-> span{F_k, k >= 0}.
# Gamma_t = (1 - P_{1/t}) U P_t is unitarily equivalent (conjugation by Mellin multipliers t^{+-i tau}) to the Hankel matrix
# H_{jl} = phi_{j+l+1} (j, l >= 0) of the symbol phi_t(tau) = t^{2 i tau} u_inf(i tau) prod_{p in S} u_p(i tau),
#   u_inf(i tau) = exp(-2 i vartheta(tau)) (Riemann-Siegel theta),  u_p(i tau) = (1 - p^{-1/2 - i tau})/(1 - p^{-1/2 + i tau}),
#   phi_m = (1/2pi) int phi e^{-i m theta} dtheta = (1/pi) int_R phi(tau) ((tau + i)/(tau - i))^m dtau/(1 + tau^2)  (real numbers).
# Coefficients: smooth circle part chi_{T1} phi by FFT on a uniform theta grid (spectrally accurate), tail (1 - chi_{T1}) phi on
# |tau| in [T1/2, T] by 16-point Gauss-Legendre panels in tau (truncation error ~ 1/(T^2 log T)).
# Singular triplets H v = kappa u; e = sum_l v_l F_{-l-1}, g = sum_j u_j F_j, f = conj(phi) g, a = (e + f)/sqrt(2 + 2 kappa) and
#   ||chi^{1/2} a||^2 = (v^T Tc v + u^T Tc u + 2 v^T X u)/(2 + 2 kappa),  Tc_{ab} = chi_{a-b} (Toeplitz), X_{lj} = (chi phi)_{l+j+1}.
# Usage: hankel.py M t primes(comma or -) R1,R2,... [T] [T1]  -> JSON (kappa_n; m_n(R) per R).
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'): os.environ[k] = '1'
import sys, json, math, time, numpy as np
from scipy.special import loggamma

def symbol(tau, t, primes):
    th = np.imag(loggamma(0.25 + 0.5j*tau)) - 0.5*tau*math.log(math.pi)
    ph = np.exp(-2j*th + 2j*tau*math.log(t))
    for p in primes:
        a = p**-0.5
        ph = ph*(1 - a*np.exp(-1j*tau*math.log(p)))/(1 - a*np.exp(1j*tau*math.log(p)))
    return ph

def chi_cut(tau, R):          # C^infty, 1 on [-R, R], 0 outside [-2R, 2R]
    x = (np.abs(tau) - R)/R
    out = np.zeros_like(tau, dtype=float)
    out[x <= 0] = 1.0
    mid = (x > 0) & (x < 1)
    a = np.exp(-1/x[mid]); b = np.exp(-1/(1 - x[mid]))
    out[mid] = b/(a + b)
    return out

def circle_coeffs(fun, Rsupp, mmax, rate_extra):
    # fun(tau) supported in |tau| <= Rsupp; c_m = (1/2pi) int fun(tau(theta)) e^{-i m theta} dtheta, m = 0..mmax-1
    rate = (math.log(Rsupp/(2*math.pi) + 1) + rate_extra)*(1 + Rsupp**2)/2 + mmax
    Nth = 2**int(math.ceil(math.log2(3*rate + 64)))
    th = 2*math.pi*np.arange(1, Nth)/Nth
    tau = -1/np.tan(th/2)
    f = np.zeros(Nth, dtype=complex); f[1:] = fun(tau)
    c = np.fft.fft(f)/Nth
    return c[:mmax], Nth

def tail_coeffs(fun, a, b, mmax, rate_extra):
    # (1/pi) int_{a<=|tau|<=b} fun(tau) z^{-m} dtau/(1+tau^2), Gauss-Legendre panels
    g, gw = np.polynomial.legendre.leggauss(16)
    rate = math.log(b/(2*math.pi) + 1) + rate_extra + 2*mmax/(1 + a*a)
    pw = 12.0/rate
    npan = int(math.ceil((b - a)/pw)); edges = np.linspace(a, b, npan + 1)
    s = ((edges[:-1, None] + edges[1:, None])/2 + (edges[1:, None] - edges[:-1, None])/2*g[None, :]).ravel()
    ws = ((edges[1:, None] - edges[:-1, None])/2*gw[None, :]).ravel()/(math.pi*(1 + s*s))
    c = np.zeros(mmax, dtype=complex)
    for sg in (1.0, -1.0):
        tau = sg*s
        cur = ws*fun(tau)
        zi = (tau + 1j)/(tau - 1j)
        for m in range(mmax):
            c[m] += cur.sum(); cur = cur*zi
    return c, s.size*2

if __name__ == '__main__':
    t0 = time.time()
    M = int(sys.argv[1]); t = float(sys.argv[2])
    primes = [] if sys.argv[3] == '-' else [int(q) for q in sys.argv[3].split(',')]
    Rs = [float(r) for r in sys.argv[4].split(',')]
    T = float(sys.argv[5]) if len(sys.argv) > 5 else 5000.0
    T1 = float(sys.argv[6]) if len(sys.argv) > 6 else 200.0
    rate_extra = 2*abs(math.log(t)) + sum(math.log(p) for p in primes)
    mmax = 2*M + 1
    sym = lambda tau: symbol(tau, t, primes)
    c1, Nth = circle_coeffs(lambda tau: chi_cut(tau, T1/2)*sym(tau), T1, mmax, rate_extra)
    c2, npt = tail_coeffs(lambda tau: (1 - chi_cut(tau, T1/2))*sym(tau), T1/2, T, mmax, rate_extra)
    phis = c1 + c2
    info = {'M': M, 't': t, 'primes': primes, 'T': T, 'T1': T1, 'Nth': Nth, 'tail_pts': npt,
            'imag_max': float(np.max(np.abs(phis.imag))), 'phi0_5': [float(x) for x in phis.real[0:6]]}
    ph = phis.real
    idx = np.arange(M)
    H = ph[idx[:, None] + idx[None, :] + 1]
    U_, S_, Vt = np.linalg.svd(H)
    V = Vt.T
    info['sec_svd'] = round(time.time() - t0, 1)
    res = {'info': info, 'kappa': [float(x) for x in S_], 'mass': {}}
    for R in Rs:
        cchi, _ = circle_coeffs(lambda tau: chi_cut(tau, R).astype(complex), 2*R, mmax, rate_extra)
        cxp, _ = circle_coeffs(lambda tau: chi_cut(tau, R)*sym(tau), 2*R, mmax, rate_extra)
        cchi = cchi.real; cxp = cxp.real
        Tc = cchi[np.abs(idx[:, None] - idx[None, :])]
        X = cxp[idx[:, None] + idx[None, :] + 1]
        ee = (V*(Tc @ V)).sum(0); gg = (U_*(Tc @ U_)).sum(0); eg = (V*(X @ U_)).sum(0)
        m = (ee + gg + 2*eg)/(2 + 2*S_)
        res['mass'][str(R)] = [float(x) for x in m]
        res['info']['chi0_R' + str(R)] = float(cchi[0])
    res['info']['sec'] = round(time.time() - t0, 1)
    print(json.dumps(res))
