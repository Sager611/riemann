# gram_herm.py: the Gram constant g_K^Herm of ACLG 2611.13129 Prop. xipoly (l. 506-510):
#   g_K = inf{ ||E(h)||^2_{L^2(R+*, d*u)} : h in H_K^Herm, ||h||_{L^2(R)} = 1 },  H_K^Herm = {sum_{k=0..K} c_k h_{4k} : int h = 0},
# h_n the L^2(R)-normalised Hermite functions for the character e^{2 pi i x y} (h_{4k} self-dual, int h_{4k} = h_{4k}(0)).
# Mellin-Plancherel (eq. mellinE l. 422): ||E h||^2 = (1/2pi) int_R |zeta(1/2+it)|^2 |Mh(1/2+it)|^2 dt, Mh(w) = int_0^inf h x^{w-1} dx,
# and ||h||^2_{L^2(R)} = 2 (1/2pi) int |Mh(1/2+it)|^2 dt. Mh computed by FFT in v = log x (trapezoid), zeta by mpmath.
# The Hermite basis is orthonormal, so G_2 = I up to quadrature error (printed) and double precision suffices.
# Usage: gram_herm.py Kmax [T]   -> JSON lines per K: g_K, top eigen-direction diagnostics.
import sys, json, math, numpy as np, mpmath as mp

Kmax = int(sys.argv[1]); T = float(sys.argv[2]) if len(sys.argv) > 2 else 340.0
nmax = 4*Kmax
v0, v1, M = -64.0, 3.6, 2**16
dv = (v1 - v0)/M
v = v0 + dv*np.arange(M)
x = np.exp(v); y = math.sqrt(2*math.pi)*x
# Hermite functions psi_n(y), keep n = 0 mod 4
p_prev = np.zeros(M); p = math.pi**-0.25*np.exp(-y*y/2)
rows = [p.copy()]
for n in range(0, nmax):
    p_next = math.sqrt(2.0/(n + 1))*y*p - math.sqrt(n/(n + 1.0))*p_prev
    p_prev, p = p, p_next
    if (n + 1) % 4 == 0: rows.append(p.copy())
H = (2*math.pi)**0.25*np.array(rows)            # h_{4k}(e^v), k = 0..Kmax
F = H*np.exp(v/2)[None, :]
# Mh(1/2 + i t_k) = dv * sum_j F_j e^{i t_k v_j}, t_k = 2 pi k/(M dv): via inverse FFT
tk = 2*math.pi*np.arange(M)/(M*dv)
Mh = dv*M*np.fft.ifft(F, axis=1)*np.exp(1j*tk*v0)[None, :]
keep = tk <= T
tk = tk[keep]; Mh = Mh[:, keep]
dt = tk[1] - tk[0]
w = np.full(tk.size, dt); w[0] = dt/2
mp.mp.dps = 20
z2 = np.array([float(abs(mp.zeta(mp.mpf(0.5) + 1j*mp.mpf(float(t))))**2) for t in tk])
G2 = (2/math.pi)*np.real((np.conj(Mh)*w[None, :]) @ Mh.T)
GE = (1/math.pi)*np.real((np.conj(Mh)*(w*z2)[None, :]) @ Mh.T)
# h_{4k}(0)
h0 = []
for k in range(Kmax + 1):
    m = 2*k
    h0.append((2*math.pi)**0.25*math.pi**-0.25*(-1)**m*math.exp(0.5*math.lgamma(2*m + 1) - m*math.log(2) - math.lgamma(m + 1)))
h0 = np.array(h0)
print(json.dumps({'info': 'G2 - I max', 'err': float(np.max(np.abs(G2 - np.eye(Kmax + 1)))), 'dt': dt, 'T': T, 'npts': int(tk.size)}), flush=True)
for K in range(1, Kmax + 1):
    vK = h0[:K + 1]/np.linalg.norm(h0[:K + 1])
    # orthonormal basis of vK-perp
    Qf, _ = np.linalg.qr(np.column_stack([vK, np.eye(K + 1)[:, :K]]))
    Q = Qf[:, 1:K + 1]
    A = Q.T @ GE[:K + 1, :K + 1] @ Q
    B = Q.T @ G2[:K + 1, :K + 1] @ Q
    L = np.linalg.cholesky(B); Li = np.linalg.inv(L)
    ev, evec = np.linalg.eigh(Li @ A @ Li.T)
    c = Q @ (Li.T @ evec[:, 0])                  # coefficients of the minimiser on h_0..h_{4K}
    dens = np.abs(c @ Mh[:K + 1])**2             # |Mh(1/2+it)|^2 of the minimiser
    tbar = float(np.sum(w*tk*dens)/np.sum(w*dens))
    print(json.dumps({'K': K, 'gK': float(ev[0]), 'g2': float(ev[1]) if K > 1 else None, 'gmax': float(ev[-1]),
                      'tmean_min': round(tbar, 2)}), flush=True)
