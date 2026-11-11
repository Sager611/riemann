# gram_lambda.py: Gram constant "beyond the Hermite regime" (ACLG 2611.13129 l. 1369, 1395; Prop. xipoly l. 506-510 at finite lambda):
#   g_K(lambda) = inf{ ||E(f)||^2_{L^2([1/lambda, lambda], d*u)} : f in H_K(lambda), ||f||_{L^2(R)} = 1 },
#   H_K(lambda) = {sum_{k=0..K} c_k h_{4k,lambda} : int f = 0},  E(f)(u) = u^{1/2} sum_{m >= 1} f(m u)  (restricted to the window).
# Prolate functions from the Bouwkamp matrix (App. A) in double precision: h_{n,lambda}(x) = lambda^{-1/2} phi_n(x/lambda), orthonormal on
# [-lambda, lambda]; int h_{n,lambda} = lambda^{1/2} sqrt2 d_0^(n). Quadrature: Gauss-Legendre in t = log u between the breakpoints
# log(lambda/m). Usage: gram_lambda.py c [Kmax]  -> JSON lines per K.
import sys, json, math, numpy as np

c = float(sys.argv[1]); lam = math.sqrt(c); g = 2*math.pi*c
Kmax = int(sys.argv[2]) if len(sys.argv) > 2 else int(c) + 2
Kl = int(2*g + 60); Kl += Kl % 2
ks = np.arange(0, Kl + 1, 2)
dg = ks*(ks + 1) + g*g*(2*ks*ks + 2*ks - 1)/((2*ks - 1)*(2*ks + 3))
of = g*g*(ks[:-1] + 1)*(ks[:-1] + 2)/((2*ks[:-1] + 3)*np.sqrt((2*ks[:-1] + 1)*(2*ks[:-1] + 5)))
Tm = np.diag(dg) + np.diag(of, 1) + np.diag(of, -1)
w, Vv = np.linalg.eigh(Tm)
nsel = [4*k for k in range(Kmax + 1)]
D = np.array([Vv[:, n//2] for n in nsel])       # rows: Legendre coefficients d^(n) on bar P_k, k even

def phis(z):                                     # phi_n(z) for n in nsel, z array in [-1, 1]
    P = np.zeros((ks.size, z.size)); Pm2 = np.ones_like(z); Pm1 = z.copy()
    # Legendre P_k by recurrence, keep even k, normalise sqrt((2k+1)/2)
    P[0] = Pm2*math.sqrt(0.5)
    for k in range(1, Kl + 1):
        if k == 1: Pk = Pm1
        else:
            Pk = ((2*k - 1)*z*Pm1 - (k - 1)*Pm2)/k
            Pm2, Pm1 = Pm1, Pk
        if k % 2 == 0: P[k//2] = Pk*math.sqrt((2*k + 1)/2)
    return D @ P

# quadrature nodes in t = log u on [-log lam, log lam], breakpoints log(lam/m), sub-panels of width <= 0.04
bps = sorted(set([-math.log(lam), math.log(lam)] + [math.log(lam/m) for m in range(1, int(c) + 1) if -math.log(lam) < math.log(lam/m) < math.log(lam)]))
gx, gw = np.polynomial.legendre.leggauss(20)
T, W = [], []
for a, b in zip(bps[:-1], bps[1:]):
    npan = max(1, int(math.ceil((b - a)/0.04)))
    e = np.linspace(a, b, npan + 1)
    for p, q in zip(e[:-1], e[1:]):
        T.append((p + q)/2 + (q - p)/2*gx); W.append((q - p)/2*gw)
T = np.concatenate(T); W = np.concatenate(W); u = np.exp(T)
# E(h_n)(u) = u^{1/2} sum_{m: m u <= lam} h_n(m u)
EH = np.zeros((len(nsel), T.size))
mmax = int(math.floor(lam/u.min())) + 1
for m in range(1, mmax + 1):
    x = m*u; sel = x <= lam
    if not sel.any(): break
    EH[:, sel] += lam**-0.5*phis(x[sel]/lam)
EH *= np.sqrt(u)[None, :]
GE = (EH*W[None, :]) @ EH.T
I = math.sqrt(lam)*math.sqrt(2)*D[:, 0]          # int h_{4k,lambda}
print(json.dumps({'c': c, 'npts': int(T.size), 'Kmax': Kmax, 'I0': float(I[0]), 'diag0_5': [float(GE[i, i]) for i in range(min(6, len(nsel)))]}), flush=True)
for K in range(1, Kmax + 1):
    v = I[:K + 1]/np.linalg.norm(I[:K + 1])
    Qf, _ = np.linalg.qr(np.column_stack([v, np.eye(K + 1)[:, :K]]))
    Q = Qf[:, 1:K + 1]
    ev = np.linalg.eigvalsh(Q.T @ GE[:K + 1, :K + 1] @ Q)
    # two-mode states f_k = I_0 h_{4k} - I_{4k} h_0, normalised in L^2(R): individual ||E f_k||^2
    fk = np.zeros(K + 1); fk[K] = I[0]; fk[0] = -I[K]; fk /= np.linalg.norm(fk)
    print(json.dumps({'K': K, 'gK': float(ev[0]), 'g2': float(ev[1]) if K > 1 else None, 'gmax': float(ev[-1]),
                      'twomode_K': float(fk @ GE[:K + 1, :K + 1] @ fk)}), flush=True)
