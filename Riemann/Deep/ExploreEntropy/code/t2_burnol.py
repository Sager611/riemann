# t2_burnol.py (explore32): check Arch(f) = (1/2) <k, (ln|y| + ln|eta|) k>, k(y) = |y|^{-1/2} f(ln|y|), k^(eta) = int k e^{2 pi i y eta} dy
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import digamma
a = 0.5
f = lambda x: np.where(np.abs(x) < a, (1 - (x/a)**2)**3, 0.0)
D = digamma(0.25) - np.log(np.pi)
def gl(lo, hi, n, m):  # m panels of n Gauss-Legendre nodes
    x0, w0 = leggauss(n); xs = []; ws = []
    e = np.linspace(lo, hi, m + 1)
    for i in range(m):
        h = (e[i+1] - e[i])/2; xs.append(e[i] + h*(x0 + 1)); ws.append(h*w0)
    return np.concatenate(xs), np.concatenate(ws)
# ||f||^2 and Arch via Levy-Khintchine: Arch = D||f||^2 + int_0^inf k(x) (2||f||^2 - 2 C_f(x)) dx
X, WX = gl(-a, a, 40, 40)
nf = np.sum(WX*f(X)**2)
def Cf(x): return np.sum(WX*f(X)*f(X - x))
kk = lambda x: np.exp(x/2)/(np.exp(x) - np.exp(-x))
T, WT = gl(0, 2*a, 40, 60)
lk = np.sum(WT*kk(T)*np.array([2*nf - 2*Cf(t) for t in T])) + 2*nf*np.sum([0.0])
# tail x >= 2a: C_f = 0: int_{2a}^inf k(x) dx * 2||f||^2 ; K(u) = sum_m e^{-b_m u}/b_m
K = lambda u: sum(np.exp(-(2*m + 0.5)*u)/(2*m + 0.5) for m in range(4000))
arch = D*nf + lk + 2*nf*K(2*a)
# Burnol side
Y, WY = gl(np.exp(-a), np.exp(a), 40, 40)
ky = Y**-0.5*f(np.log(Y))
lny = 2*np.sum(WY*np.log(Y)*ky**2)       # both signs of y
nk = 2*np.sum(WY*ky**2)
def khat(eta):
    out = np.empty_like(eta)
    for s in range(0, len(eta), 2000):
        e = eta[s:s+2000]
        out[s:s+2000] = 2*np.cos(2*np.pi*np.outer(e, Y)) @ (WY*ky)
    return out
lne = 0.0
for lo, hi, m in [(0, 1e-3, 4), (1e-3, 1e-2, 4), (1e-2, 0.1, 4), (0.1, 1, 8), (1, 10, 80), (10, 100, 800), (100, 600, 4000)]:
    E, WE = gl(lo, hi, 20, m)
    lne += 2*np.sum(WE*np.log(E)*khat(E)**2)
print('||f||^2 =', nf, ' ||k||^2 =', nk, ' ratio', nk/nf)
print('Arch(f) (Levy-Khintchine, x-space) =', arch, '  Arch/||f||^2 =', arch/nf)
print('<k,ln|y|k> =', lny, ' <k^,ln|eta|k^> =', lne, ' (1/2)<k,Hk> =', (lny + lne)/2)
print('Beckner: <k,Hk>/||k||^2 =', (lny + lne)/nk, ' >= D =', D)
