# laurent.py: band-limited, exactly unitary model of the semilocal pair of ACLG 2611.13129 Sec. 8 (no Laurent-series cut of u_p).
# x = log u, Mellin M g(i tau) = int g e^{i tau x} dx. Band |tau| <= tmax = pi/dx, samples x_j = (j + 1/2) dx; P_t <-> j <= -1 (x < 0),
# 1 - P_{1/t} <-> i >= 0, after the translation that turns phi = u into phi_t = t^{2 i tau} u (DERIVATION sec. 1 of ACLGFals).
# U = Laurent operator with symbol phi~(theta/dx) on the circle (unimodular => U unitary on l^2(Z), (P, U^*(1-P)U) an exact pair).
# phi~ = exp(i psi~): psi = -2 vartheta(tau) + 2 tau log t + sum_p 2 atan2(a sin(tau log p), 1 - a cos(tau log p)), a = p^{-1/2}
# (= arg u_p(i tau), u_p = (1 - p^{-1/2 - i tau})/(1 - p^{-1/2 + i tau})); psi~ = psi on |tau| <= tau1 = b1 tmax, then a smooth blend
# into the line through (tau1, psi(tau1)) and (tmax, pi K) (negative slope = negative delay: no spurious self-dual modes; smooth
# and odd on the circle). Gamma = Hankel matrix H_{il} = c_{i+l+1}, c_m = (1/2pi) int phi~(theta/dx) e^{-i m theta} dtheta, i, l < N,
# N = depth/dx. H symmetric: H v = mu v, kappa = |mu|, e = v (input l <-> x = -(l+1/2)dx), g = sign(mu) v (output i <-> (i+1/2)dx).
# M e(i tau) = sqrt(dx) sum_l e_l e^{-i tau (l+1/2) dx}, M g = sqrt(dx) sum_i g_i e^{i tau (i+1/2) dx}, M f = conj(phi) M g,
# m_n(R) = (1/2pi) int chi_R |M e + M f|^2 dtau/(2 + 2 kappa), T(n0, R) = sum_{n>n0} kappa_n m_n(R)  (ACLG eq. H2ii).
# Usage: laurent.py primes(comma or -) tmax depth R1,R2.. [t] [b1]   env: LA_EIGH=k pairs per end (400), LA_F32=1
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'): os.environ[k] = '1'
import sys, json, math, time, numpy as np
import scipy.linalg as sla
from scipy.special import loggamma
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hankel import chi_cut
t0 = time.time()
primes = [] if sys.argv[1] == '-' else [int(q) for q in sys.argv[1].split(',')]
tmax = float(sys.argv[2]); depth = float(sys.argv[3]); Rs = [float(r) for r in sys.argv[4].split(',')]
tsc = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0
b1 = float(sys.argv[6]) if len(sys.argv) > 6 else 0.75
DT = np.float32 if os.environ.get('LA_F32', '0') == '1' else np.float64
dx = math.pi/tmax; N = int(round(depth/dx))
def psi(tau):
    th = np.imag(loggamma(0.25 + 0.5j*tau)) - 0.5*tau*math.log(math.pi)
    out = -2*th + 2*tau*math.log(tsc)
    for p in primes:
        a = p**-0.5; L = math.log(p)
        out = out + 2*np.arctan2(a*np.sin(tau*L), 1 - a*np.cos(tau*L))
    return out
tau1 = b1*tmax; tau2 = 0.5*(1 + b1)*tmax
p1 = float(psi(np.array([tau1]))[0]); d1 = float((psi(np.array([tau1 + 1e-4])) - psi(np.array([tau1 - 1e-4])))[0]/2e-4)
K = round((p1 + min(d1, -0.5)*(tmax - tau1))/math.pi)
slope = (math.pi*K - p1)/(tmax - tau1)
def step(x):      # smooth 0 -> 1 on [0, 1]
    x = np.clip(x, 0, 1); a = np.exp(-1/np.maximum(x, 1e-300)); b = np.exp(-1/np.maximum(1 - x, 1e-300))
    return np.where(x <= 0, 0.0, np.where(x >= 1, 1.0, a/(a + b)))
def psit(tau):    # closed phase, odd
    s = np.sign(tau); a = np.abs(tau)
    bl = step((a - tau1)/(tau2 - tau1))
    lin = p1 + slope*(a - tau1)
    return s*((1 - bl)*psi(a) + bl*lin)
Nth = 1 << int(math.ceil(math.log2(4*N + 4096)))
k = np.arange(Nth); th = 2*math.pi*k/Nth; th = np.where(th >= math.pi, th - 2*math.pi, th)
c = (np.fft.fft(np.exp(1j*psit(th/dx)))/Nth)
cim = float(np.max(np.abs(c.imag[:2*N + 2]))); c = c.real
tail = float(np.max(np.abs(c[2*N - 50:2*N + 1])))
def build_H():
    H = np.empty((N, N), dtype=DT, order='F'); ar = np.arange(N)
    for i0 in range(0, N, 256):
        i1 = min(N, i0 + 256)
        H[i0:i1] = c[1 + np.arange(i0, i1)[:, None] + ar[None, :]]
    return H
kE = min(int(os.environ.get('LA_EIGH', '400')), N//2)
H = build_H()
mu1, v1 = sla.eigh(H, subset_by_index=[N - kE, N - 1], driver='evr', overwrite_a=True)
del H; H = build_H()
mu2, v2 = sla.eigh(H, subset_by_index=[0, kE - 1], driver='evr', overwrite_a=True)
del H
cut = max(abs(float(mu1[0])), abs(float(mu2[-1])))
mu = np.concatenate([mu1, mu2]); o = np.argsort(-np.abs(mu)); o = o[np.abs(mu[o]) >= cut]
V = np.concatenate([v1, v2], axis=1)[:, o].astype(np.float64); del v1, v2
mu = mu[o].astype(np.float64); S = np.abs(mu); sgn = np.sign(mu)
# depth profile: fraction of |e_n|^2 beyond depth 10, 20 (x < -10, -20)
xl = (np.arange(N) + 0.5)*dx
deep10 = (V[xl > 10]**2).sum(0); deep20 = (V[xl > 20]**2).sum(0); xmean = (xl[:, None]*V**2).sum(0)
res = {'info': {'model': 'laurent', 'primes': primes, 't': tsc, 'tmax': tmax, 'tau1': tau1, 'depth': depth, 'dx': dx, 'N': N,
                'Nth': Nth, 'K': K, 'slope': slope, 'c_imag_max': cim, 'c_tail': tail, 'kE': kE, 'cut': cut, 'f32': DT == np.float32,
                'sec_eig': round(time.time() - t0, 1)},
       'kappa': [float(x) for x in S], 'mass': {}, 'xmean': [float(x) for x in xmean], 'deep10': [float(x) for x in deep10],
       'deep20': [float(x) for x in deep20]}
gx, gw = np.polynomial.legendre.leggauss(20)
for R in Rs:
    tp = max(4, int(math.ceil(2*R/1.0)))
    te = np.linspace(0, 2*R, tp + 1)
    tau = ((te[:-1, None] + te[1:, None])/2 + (te[1:, None] - te[:-1, None])/2*gx[None, :]).ravel()
    tw = ((te[1:, None] - te[:-1, None])/2*gw[None, :]).ravel()
    chi = chi_cut(tau, R); ph = np.exp(1j*psi(tau))
    m = np.zeros(S.size)
    for c0 in range(0, tau.size, 100):
        c1 = min(tau.size, c0 + 100)
        E = math.sqrt(dx)*np.exp(-1j*np.outer(tau[c0:c1], xl))          # input: x = -(l+1/2) dx
        Er = np.ascontiguousarray(E.real); Ei = np.ascontiguousarray(E.imag); del E
        A1 = Er @ V; A2 = Ei @ V; del Er, Ei
        Me = A1 + 1j*A2; Mg = (A1 - 1j*A2)*sgn[None, :]                    # output: x = +(i+1/2) dx, same samples
        Ma = Me + np.conj(ph[c0:c1])[:, None]*Mg
        m += (2/(2*math.pi))*((tw*chi)[c0:c1][None, :] @ np.abs(Ma)**2).ravel()
        del A1, A2, Me, Mg, Ma
    m /= (2 + 2*S)
    res['mass'][str(R)] = [float(x) for x in m]
res['info']['sec'] = round(time.time() - t0, 1)
print(json.dumps(res))
