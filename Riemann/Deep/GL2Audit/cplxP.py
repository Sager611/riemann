# gl2audit: which P-class does the papers' inf_P formula (B l.796, AT l.996: sqrt(AC - Im(B)^2)) belong to?
# direct double-integral c(P,Q) for a complex Q of Bernard's shape; minimise over real P and over complex P (poly basis)
import numpy as np
from numpy.polynomial import chebyshev as Ch
from scipy.optimize import minimize
xi = 25/256; R = 3.0
x, w = np.polynomial.legendre.leggauss(60); U = (x+1)/2; WU = w/2
def T(n, t): c = np.zeros(n+1); c[n] = 1; return Ch.chebval(t, c), Ch.chebval(t, Ch.chebder(c))
t = 1-2*U
T1, dT1 = T(1, t); T3, dT3 = T(3, t); T2, dT2 = T(2, t)
# Bernard's shape: Q = 1 + i^2 l1[(1-2x)-1] + i^3 l2[(1-2x)^2-1] + i^4 l3[(1-2x)^3-1]
l1, l2, l3 = 0.9, 0.35, 0.2
Q = 1 - l1*(T1-1) - 1j*l2*(t**2-1) + l3*(t**3-1)
dQ = -l1*(-2) - 1j*l2*(2*t)*(-2) + l3*(3*t**2)*(-2)
wv = np.exp(R*U)*Q; dwv = np.exp(R*U)*(R*Q+dQ)
A = np.sum(WU*abs(wv)**2); C = np.sum(WU*abs(dwv)**2); B = np.sum(WU*wv*np.conj(dwv)); w1 = abs(np.exp(R)*(1 - l1*(-2) - 1j*l2*0 + l3*(-2)))**2
NB = 8
def Pv(a, b):
    co = a + 1j*b; s = np.zeros_like(U, dtype=complex); ds = np.zeros_like(U, dtype=complex)
    for j in range(NB):
        Tj, dTj = T(j, 2*U-1); s += co[j]*Tj; ds += co[j]*2*dTj
    g = U*(1-U); dg = 1-2*U
    return U + g*s, 1 + dg*s + g*ds
def cdirect(a, b):
    P, dP = Pv(a, b)
    I = dP[:, None]*Q[None, :] + xi*P[:, None]*(dQ+R*Q)[None, :]
    return 1 + np.sum(WU[:, None]*WU[None, :]*np.exp(2*R*U)[None, :]*abs(I)**2)/xi
orr = minimize(lambda a: cdirect(a, np.zeros(NB)), np.zeros(NB), method='BFGS', options={'gtol': 1e-12})
occ = minimize(lambda v: cdirect(v[:NB], v[NB:]), np.zeros(2*NB) + 1e-3, method='BFGS', options={'gtol': 1e-12})
X = np.sqrt(A*C - B.imag**2)
print(f"Q complex (Im B = {B.imag:.6f}, A = {A:.6f}, C = {C:.6f}), R = {R}, xi = 25/256")
print(f"min over REAL P (direct c, deg {NB+1}):    {orr.fun:.10f}")
print(f"formula sqrt(AC) coth (real-P claim):       {(1+w1)/2 + np.sqrt(A*C)/np.tanh(xi*np.sqrt(C/A)):.10f}")
print(f"min over COMPLEX P (direct c):              {occ.fun:.10f}")
print(f"formula sqrt(AC - ImB^2) (B l.796, AT l.996): {(1+w1)/2 + X/np.tanh(xi*X/A):.10f}")
