# ceilaudit9 E45: twisted means and a two-piece split against Levinson's concavity (Jensen) step, explore45's setting
# (zeta, theta = 1/2, P(x) = x, Q(x) = 1 - b x, b = 1.03217, R = 1.24618, y = T0^theta; V = zeta + (b/L) zeta'; line sigma = a = 1/2 - R/L).
# U(Psi) = (1/2) log mean|V Psi|^2 is the Jensen bound for m = mean log|V Psi|; concavity loss J = U - m.
# (i) twist: U_B = (1/2) log mean|V Psi B|^2 - mean log|B|, B = 1 + sum_{k=2..K} beta_k k^{-s} (beta complex, optimised); identity: = Jensen bound of
#     the mollifier Psi*B (length yK) -> compare with the plain Levinson mollifier of length yK. (ii) split: per-piece Jensen on E = {cos(t log p + phi0) > 0}
#     and its complement, and on two halves in t.
import sys, math, numpy as np
from scipy.special import bernoulli
from scipy.optimize import minimize
T0 = float(sys.argv[1]); H = float(sys.argv[2]); dt = float(sys.argv[3])
theta, b, R = 0.5, 1.03217, 1.24618
Tm = T0 + H/2; L = math.log(Tm/(2*math.pi)); a = 0.5 - R/L; y = T0**theta
K = 12; B2 = bernoulli(2*K)[2::2]; fac = [math.factorial(2*k) for k in range(1, K+1)]
N = int((T0+H)/math.pi) + 60; M = N - 1; ms = np.arange(1, M+1, dtype=float); logm = np.log(ms)
def Htail(s, X):
    lX = np.log(X); Xs = np.exp(-s*lX)
    h = X*Xs/(s-1) + Xs/2; hs = -lX*X*Xs/(s-1) - X*Xs/(s-1)**2 - lX*Xs/2
    poch = s.copy(); dpoch = np.ones_like(s)
    for k in range(1, K+1):
        if k > 1:
            for i in (2*k-3, 2*k-2):
                dpoch = dpoch*(s+i) + poch; poch = poch*(s+i)
        term = B2[k-1]/fac[k-1]*X**(-(2*k-1))*Xs
        h = h + term*poch; hs = hs + term*(dpoch - lX*poch)
    return h, hs
def mob(n):
    mu = np.ones(n+1, int); mu[0] = 0; isp = np.ones(n+1, bool); isp[:2] = False
    for p in range(2, n+1):
        if isp[p]:
            isp[2*p::p] = False; mu[p::p] *= -1; mu[p*p::p*p] = 0
    return mu
ylist = [y*k for k in (1, 2, 3, 4)]; nymax = int(max(ylist)); mu = mob(nymax)
cps = []
for yy in ylist:
    nn = np.arange(1, int(yy)+1, dtype=float); cps.append(mu[1:int(yy)+1]*nn**(a-0.5)*(np.log(yy/nn)/math.log(yy)))
J = int(round(H/dt)); ts = T0 + dt*np.arange(J); Bk = 100
Om = np.exp(-1j*dt*np.outer(logm, np.arange(Bk)))
Wm = np.stack([np.exp(-a*logm), -logm*np.exp(-a*logm)])
f = [np.zeros(J, complex) for _ in ylist]
for j0 in range(0, J, Bk):
    nb = min(Bk, J - j0); tb = ts[j0:j0+nb]
    ph = np.exp(-1j*tb[0]*logm)[:, None]*Om[:, :nb]
    S = Wm @ ph; s = a + 1j*tb; h, hs = Htail(s, N)
    V = (S[0] + h) + (b/L)*(S[1] + hs)
    for i, cp in enumerate(cps):
        f[i][j0:j0+nb] = V*((cp*np.exp(-a*logm[:len(cp)])) @ ph[:len(cp)])
U = [0.5*np.log(np.mean(np.abs(g)**2)) for g in f]; m0 = np.mean(np.log(np.abs(f[0])))
print(f"T0={T0:.0e} H={H} dt={dt} L={L:.3f} a={a:.4f} y={y:.0f}: U(Psi_y)={U[0]:.5f}, m={m0:.5f}, concavity loss J={U[0]-m0:.5f}; "
      f"plain Levinson mollifier length yK, K=2,3,4: U={U[1]:.5f},{U[2]:.5f},{U[3]:.5f} (gains {U[0]-U[1]:+.5f},{U[0]-U[2]:+.5f},{U[0]-U[3]:+.5f})", flush=True)
g = f[0]; A2 = np.abs(g)**2
def UB(x, Kt):
    Bv = np.ones(J, complex)
    for k in range(2, Kt+1):
        Bv += (x[2*(k-2)] + 1j*x[2*(k-2)+1])*np.exp(-(a + 1j*ts)*math.log(k))
    return 0.5*np.log(np.mean(A2*np.abs(Bv)**2)) - np.mean(np.log(np.abs(Bv)))
for Kt in (2, 3, 4, 6):
    best = None
    for st in range(4):
        x0 = 0.3*np.random.default_rng(st).standard_normal(2*(Kt-1))
        r = minimize(UB, x0, args=(Kt,), method='Nelder-Mead', options={'maxiter': 3000, 'xatol': 1e-6, 'fatol': 1e-9})
        if best is None or r.fun < best.fun: best = r
    bet = best.x[0::2] + 1j*best.x[1::2]
    print(f"  twist K={Kt}: min U_B={best.fun:.5f}, gain {U[0]-best.fun:+.5f} = {(U[0]-best.fun)/(U[0]-m0):.3%} of J; |beta|={np.round(np.abs(bet),3)}; "
          f"vs plain length-{Kt}y mollifier gain {U[0]-U[min(Kt,4)-1]:+.5f}" + (" (K=6: compared with length 4y)" if Kt == 6 else ""), flush=True)
for p in (2, 3, 5):
    bestg = None
    for phi0 in np.arange(16)*2*math.pi/16:
        E = np.cos(ts*math.log(p) + phi0) > 0; w1 = E.mean()
        Us = w1*0.5*np.log(A2[E].mean()) + (1-w1)*0.5*np.log(A2[~E].mean())
        if bestg is None or Us < bestg[0]: bestg = (Us, phi0)
    print(f"  split cos(t log {p} + phi0) > 0: best U_split={bestg[0]:.5f}, gain {U[0]-bestg[0]:+.5f} = {(U[0]-bestg[0])/(U[0]-m0):.3%} of J", flush=True)
hlf = ts < T0 + H/2; Us = 0.5*(0.5*np.log(A2[hlf].mean()) + 0.5*np.log(A2[~hlf].mean()))
print(f"  split two halves in t: gain {U[0]-Us:+.6f} = {(U[0]-Us)/(U[0]-m0):.3%} of J", flush=True)
c2 = 1 - 2**(a-0.5)*(1 - math.log(2)/math.log(y)) - (b/L)*math.log(2)
print(f"  coefficient of 2^-s in V*Psi: c2 = {c2:.4f} (-> 0 as T -> infinity at fixed theta: twisted first-order gain -> 0)", flush=True)
