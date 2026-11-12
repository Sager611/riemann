# Finite-T Levinson experiment: f = V*Psi, V = F + (b/L) F', F = zeta or Davenport-Heilbronn.
# usage: python3 levin.py zeta|dh T0 H dt  -> writes levin_<fn>.npz and prints summary
import numpy as np, sys, math, time
from scipy.special import bernoulli
fn = sys.argv[1]; T0 = float(sys.argv[2]); H = float(sys.argv[3]); dt = float(sys.argv[4])
theta, b, R = 0.5, 1.03217, 1.24618
Tm = T0 + H/2; L = math.log(Tm/(2*math.pi)); a = 0.5 - R/L; y = T0**theta
K = 12; B2 = bernoulli(2*K)[2::2]; fac = [math.factorial(2*k) for k in range(1, K+1)]
if fn == 'zeta':
    q = 1; classes = [(1.0, 0.0)]
    N = int((T0+H)/math.pi) + 60; M = N - 1
    coef = np.ones(M)
else:
    q = 5; kap = (math.sqrt(10-2*math.sqrt(5)) - 2)/(math.sqrt(5) - 1)
    per = [1.0, kap, -kap, -1.0, 0.0]          # a_m for m = 1..5 mod 5
    classes = [(per[c-1], c/5) for c in range(1, 5)]
    N = int((T0+H)/math.pi) + 60; M = 5*N - 1
    coef = np.array([per[(m-1) % 5] for m in range(1, M+1)])
ms = np.arange(1, M+1, dtype=float); logm = np.log(ms)
# mollifier coefficients: Dirichlet inverse of F's coefficients, times n^{a-1/2} P(log(y/n)/log y), P(x)=x
ny = int(y); A = np.zeros(ny+1); A[1:] = [coef[n-1] if n <= M else 0 for n in range(1, ny+1)]
inv = np.zeros(ny+1); inv[1] = 1.0
for n in range(2, ny+1):
    s_ = 0.0
    for d in range(1, n):
        if n % d == 0: s_ += inv[d]*A[n//d]
    inv[n] = -s_
nn = np.arange(1, ny+1, dtype=float)
cpsi = inv[1:]*nn**(a-0.5)*(np.log(y/nn)/math.log(y))
xs = np.array([-4, -3, -2.5, -2, -1.6, -R, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 9, 14, 2.5*L])
if len(sys.argv) > 5 and sys.argv[5] == 'dense':
    xs = np.array([-12, -10, -8, -7, -6, -5.5, -5, -4.5, -4, -3.5, -3, -2.75, -2.5, -2.25, -2, -1.8, -1.6, -1.4, -R, -1.1, -1, -0.9, -0.75, -0.6, -0.5, -0.35, -0.25, -0.1, 0, 0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 9, 14, 2.5*L])
sig = 0.5 + xs/L; nl = len(xs)
def Htail(s, X):
    lX = np.log(X); Xs = np.exp(-s*lX)
    h = X*Xs/(s-1) + Xs/2; hs = -lX*X*Xs/(s-1) - X*Xs/(s-1)**2 - lX*Xs/2
    poch = s.copy(); dpoch = np.ones_like(s)    # (s)_{2k-1} and derivative
    for k in range(1, K+1):
        if k > 1:
            for i in (2*k-3, 2*k-2):
                dpoch = dpoch*(s+i) + poch; poch = poch*(s+i)
        term = B2[k-1]/fac[k-1]*X**(-(2*k-1))*Xs
        h = h + term*poch; hs = hs + term*(dpoch - lX*poch)
    return h, hs
def tail(s):
    h = 0; hs = 0
    for w, al in classes:
        hh, hhs = Htail(s, N + al)
        qs = q**(-s) if q > 1 else 1.0
        h = h + w*qs*hh; hs = hs + w*qs*(hhs - (math.log(q)*hh if q > 1 else 0))
    return h, hs
J = int(round(H/dt)); ts = T0 + dt*np.arange(J)
Bk = 100 if q == 1 else 40
Om = np.exp(-1j*dt*np.outer(logm, np.arange(Bk)))
Wm = np.concatenate([coef*np.exp(-np.outer(sig, logm)), -coef*logm*np.exp(-np.outer(sig, logm))])  # 2nl x M
Wp = cpsi*np.exp(-np.outer(sig, np.log(nn)))
logV = np.zeros((nl, J), np.float32); logP = np.zeros((nl, J), np.float32)
Vhalf = np.zeros(J, complex); Phalf = np.zeros(J, complex); Fhalf = np.zeros(J, complex)
i0 = int(np.argmin(abs(xs)))
t_start = time.time()
for j0 in range(0, J, Bk):
    jj = np.arange(j0, min(J, j0+Bk)); tb = ts[jj]; nb = len(jj)
    ph = np.exp(-1j*tb[0]*logm)[:, None]*Om[:, :nb]
    S = Wm @ ph
    s = sig[:, None] + 1j*tb[None, :]
    h, hs = tail(s)
    F = S[:nl] + h; Fp = S[nl:] + hs
    V = F + (b/L)*Fp
    P = Wp @ ph[:ny]
    logV[:, jj] = np.log(np.abs(V)); logP[:, jj] = np.log(np.abs(P))
    Vhalf[jj] = V[i0]; Phalf[jj] = P[i0]; Fhalf[jj] = F[i0]
print('main loop s', round(time.time()-t_start, 1), flush=True)
del Om, Wm, ph, S, F, Fp, V, P, h, hs, s
def evalV_sig(sigs, t):
    ph = np.exp(-1j*t*logm); Vs = []; Ps = []
    for c0 in range(0, len(sigs), 25):
        sc = sigs[c0:c0+25]; E = np.exp(-np.outer(sc, logm))
        S = (coef*E) @ ph; Sp = (-coef*logm*E) @ ph; del E
        s = sc + 1j*t; h, hs = tail(s)
        Vs.append(S + h + (b/L)*(Sp + hs))
        Ps.append((cpsi*np.exp(-np.outer(sc, np.log(nn)))) @ ph[:ny])
    return np.concatenate(Vs), np.concatenate(Ps)
sg = np.linspace(0.5, 0.5 + 2.5, 3001)
Vb, Pb = evalV_sig(sg, T0); Vt, Pt = evalV_sig(sg, T0 + dt*(J-1))
np.savez_compressed(f'levin_{fn}' + ('_dense' if len(sys.argv) > 5 else '') + '.npz', xs=xs, sig=sig, logV=logV, logP=logP, Vhalf=Vhalf, Phalf=Phalf, Fhalf=Fhalf,
    Vb=Vb, Pb=Pb, Vt=Vt, Pt=Pt, sg=sg, T0=T0, H=H, dt=dt, L=L, a=a, b=b, R=R, theta=theta, y=y)
print('saved', fn, 'N', N, 'M', M, 'L', L, 'a', a, 'y', y)
