# semi_pair.py: the semilocal pair of ACLG 2611.13129 Sec. 8 at t = 1 as an integral operator on L^2(0,1) (derivation: ../DERIVATION.md
# sec. 3). Gamma = (1 - P) U P is unitarily equivalent to (A G)(w) = 2 int_0^1 C_S(u w) G(u) du with
#   C_S(z) = sum_{n} gamma(n) cos(2 pi n z),  n = prod_{p in S} p^{k_p} (k_p >= -1),  gamma(n) = prod_p g_p(k_p),
#   g_p(k) = 1 - 1/p (k >= 0), -1/p (k = -1); truncated at n <= nmax (frequency cut). S = {infinity}: C = cos(2 pi z) (prolate, gamma = 2 pi).
# Input g(u) = u^{1/2} G(u) on (0,1]; output g_out(v) = v^{-1/2} H(1/v) on [1, inf); A G_n = kappa_n H_n.
# Mellin: M e_n(i tau) = int G_n u^{-1/2 + i tau} du, M g_n(i tau) = int H_n w^{-1/2 - i tau} dw, M f_n = conj(u(i tau)) M g_n,
# a_n = (e_n + f_n)/sqrt(2 + 2 kappa_n); m_n(R) = (1/2pi) int chi_R |M a_n|^2 dtau;  T(n0, R) = sum_{n > n0} kappa_n m_n(R) (eq. H2ii).
# Nystrom: composite Gauss-Legendre in s, u = s^2. Usage: semi_pair.py primes(comma or -) nmax R1,R2,... [nodes_per_wavelength] [taper]
# env: SP_EIGH=k (symmetric eigensolver, k pairs per end; else full SVD), SP_F32=1 (float32 kernel), SP_CONSISTENT=0 (exact u in f = U^* g).
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'): os.environ[k] = '1'
import sys, json, math, time, itertools, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hankel import symbol, chi_cut
from scipy.special import loggamma
consistent = os.environ.get('SP_CONSISTENT', '1') == '1'
DT = np.float32 if os.environ.get('SP_F32', '0') == '1' else np.float64   # float32 kernel + eigensolver (memory cap)

t0 = time.time()
primes = [] if sys.argv[1] == '-' else [int(q) for q in sys.argv[1].split(',')]
nmax = float(sys.argv[2]); Rs = [float(r) for r in sys.argv[3].split(',')]
npw = float(sys.argv[4]) if len(sys.argv) > 4 else 7.0
taper = (sys.argv[5] == 'taper') if len(sys.argv) > 5 else False   # smooth cut of gamma(n) on [nmax/4, nmax] (truncation-artifact test)
# frequencies n and weights gamma(n)
terms = {1.0: 1.0}
for p in primes:
    new = {}
    for n, gm in terms.items():
        k = -1
        while True:
            m = n*p**k
            if m > nmax: break
            new[m] = new.get(m, 0.0) + gm*((1 - 1/p) if k >= 0 else -1/p)
            k += 1
    terms = new
freqs = np.array(sorted(terms)); gams = np.array([terms[f] for f in freqs])
if taper:
    xx = np.clip((np.log(freqs) - math.log(nmax/4))/math.log(4), 0, 1)
    gams = gams*np.where(xx <= 0, 1.0, np.where(xx >= 1, 0.0, np.exp(-1/np.maximum(1 - xx, 1e-12))/(np.exp(-1/np.maximum(1 - xx, 1e-12)) + np.exp(-1/np.maximum(xx, 1e-12)))))
# grid: composite GL (q = 20) in s on [0,1], graded panels covering <= 3 wavelengths of cos(2 pi nmax s^2): width <= 1.5/(nmax s)
q = 20
gx, gw = np.polynomial.legendre.leggauss(q)
edges = [0.0]
while edges[-1] < 1.0:
    a = edges[-1]; wdt = min(0.05, 3.0*7.0/(npw*2*max(nmax, 1.0)*max(a, 1e-3)))
    edges.append(min(1.0, a + wdt))
e = np.array(edges); npan = e.size - 1
s = ((e[:-1, None] + e[1:, None])/2 + (e[1:, None] - e[:-1, None])/2*gx[None, :]).ravel()
ws = ((e[1:, None] - e[:-1, None])/2*gw[None, :]).ravel()
u = s*s; om = 2*s*ws                                # du-weights
n = u.size
sq = np.sqrt(om)
def build_K():
    K = np.empty((n, n), dtype=DT, order='F')      # Fortran order: LAPACK works in place (no copy)
    for i0 in range(0, n, 128):                        # row blocks: no n x n temporaries
        i1 = min(n, i0 + 128)
        blk = np.zeros((i1 - i0, n)); arg = 2*math.pi*np.outer(u[i0:i1], u)
        for f, gm in zip(freqs, gams): blk += (2*gm)*np.cos(f*arg)
        K[i0:i1] = blk*(sq[i0:i1, None]*sq[None, :])
    return K
K = build_K()
kE = int(os.environ.get('SP_EIGH', '0'))
if kE > 0:
    # K is symmetric (same grid for input and output): eigenpairs at both ends, singular values |mu|, H = sign(mu) V
    import scipy.linalg as sla
    nK = K.shape[0]; kE = min(kE, nK//2)
    mu1, v1 = sla.eigh(K, subset_by_index=[nK - kE, nK - 1], driver='evr', overwrite_a=True)
    del K; K = build_K()                            # rebuilt instead of copied (memory cap)
    mu2, v2 = sla.eigh(K, subset_by_index=[0, kE - 1], driver='evr', overwrite_a=True)
    del K
    mu = np.concatenate([mu1, mu2]); Vall = np.concatenate([v1, v2], axis=1); del v1, v2
    o = np.argsort(-np.abs(mu)); mu = mu[o]; V = Vall[:, o]; del Vall
    V = V.astype(np.float64); mu = mu.astype(np.float64)
    S = np.abs(mu); Hm = V*np.sign(mu)[None, :]
else:
    Hm, S, Vt = np.linalg.svd(K)                    # K[i, j]: output w_i, input u_j;  K V = S H
    del K
    V = Vt.T
res = {'info': {'taper': taper, 'primes': primes, 'nmax': nmax, 'nterms': int(freqs.size), 'nodes': int(n), 'sec_svd': round(time.time() - t0, 1)},
       'kappa': [float(x) for x in S], 'mass': {}, 'check': {}}
lu = np.log(u)
Gm = V/sq[:, None]; Hn = Hm/sq[:, None]                                # G_n(u_j), H_n(w_i)
del V, Hm
for R in Rs:
    # tau quadrature on [0, 2R]: GL panels of width <= 0.5 (|log u| <= ~25)
    tp = max(4, int(math.ceil(2*R/1.5)))
    te = np.linspace(0, 2*R, tp + 1)
    tau = ((te[:-1, None] + te[1:, None])/2 + (te[1:, None] - te[:-1, None])/2*gx[None, :]).ravel()
    tw = ((te[1:, None] - te[:-1, None])/2*gw[None, :]).ravel()
    chi = chi_cut(tau, R)
    th = np.imag(loggamma(0.25 + 0.5j*tau)) - 0.5*tau*math.log(math.pi)
    uK = np.exp(-2j*th)*sum(gm*f**-0.5*np.exp(1j*tau*math.log(f)) for f, gm in zip(freqs, gams))
    ub = np.conj(uK) if consistent else np.conj(symbol(tau, 1.0, primes))
    # product-integration weights W[tau, j] = int_panel l_j(u) u^{-1/2 + i tau} du (l_j Lagrange basis on the panel's q nodes in s),
    # fine sub-grid in log s per panel (resolves tau log u for tau <= 2R); G entire in u, so degree q-1 interpolation per panel.
    W = np.zeros((tau.size, n), dtype=complex)
    g16, w16 = np.polynomial.legendre.leggauss(16)
    for pi_ in range(npan):
        sa, sb = e[pi_], e[pi_ + 1]; nodes = s[pi_*q:(pi_ + 1)*q]
        ya, yb = math.log(max(sa, 1e-7)), math.log(sb)
        npy = max(1, int(math.ceil((yb - ya)*4*R/6.0)))                # 16-pt GL panels in y = log s, <= 6 rad of 2 tau y
        ye = np.linspace(ya, yb, npy + 1)
        y = ((ye[:-1, None] + ye[1:, None])/2 + (ye[1:, None] - ye[:-1, None])/2*g16[None, :]).ravel()
        wy = ((ye[1:, None] - ye[:-1, None])/2*w16[None, :]).ravel()
        sf = np.exp(y); wf = wy*sf                                     # ds = s dy
        bw = np.array([1.0/np.prod([nodes[j] - nodes[k] for k in range(q) if k != j]) for j in range(q)])
        diff = sf[:, None] - nodes[None, :]
        diff[diff == 0] = 1e-300
        L = (bw[None, :]/diff)/np.sum(bw[None, :]/diff, axis=1, keepdims=True)
        LW = L*(2*wf)[:, None]                                         # u^{-1/2 + i tau} du = 2 s^{2 i tau} ds
        for c0 in range(0, tau.size, 256):
            c1 = min(tau.size, c0 + 256)
            W[c0:c1, pi_*q:(pi_ + 1)*q] = np.exp(2j*np.outer(tau[c0:c1], y)) @ LW
    m = np.zeros(S.size); cross = np.zeros(S.size)

    TB = 300
    for c0 in range(0, tau.size, TB):
        c1 = min(tau.size, c0 + TB)
        Wr = np.ascontiguousarray(W[c0:c1].real); Wi = np.ascontiguousarray(W[c0:c1].imag)   # real products: no complex upcast of Gm, Hn
        Me = (Wr @ Gm) + 1j*(Wi @ Gm); Mg = (Wr @ Hn) - 1j*(Wi @ Hn); del Wr, Wi
        wt = (tw*chi)[c0:c1]
        m += (2/(2*math.pi))*(wt[None, :] @ (np.abs(Me + ub[c0:c1, None]*Mg)**2)).ravel()
        cross += (2/(2*math.pi))*(wt[None, :] @ np.real(np.conj(Me)*ub[c0:c1, None]*Mg)).ravel()
        del Me, Mg
    del W
    m /= (2 + 2*S)
    res['mass'][str(R)] = [float(x) for x in m]
    res['check'][str(R)] = {'cross_over_kappa_first5': [float(cross[i]/S[i]) for i in range(5)]}
res['info']['sec'] = round(time.time() - t0, 1)
print(json.dumps(res))
