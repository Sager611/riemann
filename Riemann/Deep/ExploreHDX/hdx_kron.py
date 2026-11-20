# hdx_kron.py (explore64): sanity check of the Kronecker ceiling on actual functions f_T = u(x) cos(T x) (float64, x-space quadrature).
# Natural hard-link localization (1 - m(z) = min(z/2a', 1)). Prediction: R(f_T)/|f_T|^2 -> -sum_n 2 w_n (1-m(ln n)) cos(T ln n) <u,T_n u>/|u|^2
# (smooth parts -> 0 as T grows). Also Q(f_T)/|f_T|^2 = local part + R: the local part carries the archimedean log T growth.
import sys, math, numpy as np
from scipy.special import digamma

def ppow(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n; e = 0
            while m % p == 0: m //= p; e += 1
            return (p, e) if m == 1 else None
    return None

def run(c, cl, Ts=None, M=40001):
    a = math.log(c)/2; ap = math.log(cl)/2
    x = np.linspace(-a, a, M); dx = x[1] - x[0]
    u = (1 - (x/a)**2)**2
    at = [(math.log(n), math.log(ppow(n)[0])/math.sqrt(n)) for n in range(2, int(math.ceil(c))) if n < c and ppow(n)]
    onem = lambda z: np.minimum(np.abs(z)/(2*ap), 1.0)
    def corr(f):   # C_f(z_k) = int f(x) f(x - z_k) dx, z_k = k dx, k = 0..M-1
        L = 1 << int(math.ceil(math.log2(2*M))); Fh = np.fft.rfft(f, L)
        return np.fft.irfft(Fh*np.conj(Fh), L)[:M]*dx
    Cu = corr(u); nu = np.sum(u*u)*dx
    def at_C(C, s):
        k = s/dx; k0 = int(k); t = k - k0; return (1 - t)*C[k0] + t*C[k0 + 1]
    W = [(s, w, 2*w*float(onem(s))*at_C(Cu, s)/nu) for s, w in at]
    if Ts is None:
        Tg = np.linspace(50, 400, 350001); A = sum(Wn*np.cos(Tg*s) for s, w, Wn in W)
        Ts = [float(Tg[np.argmax(A)])]
    z = np.arange(M)*dx
    kap = np.where(z > 0, np.exp(z/2)/(np.exp(z) - np.exp(-z) + 1e-300), 0.0)
    ksm = onem(z)*(2*np.cosh(z/2) - kap); ksm[0] = -1/(4*ap)   # limit (z/2a')(-1/(2z)) at 0
    out = []
    for T in Ts:
        f = u*np.cos(T*x); Cf = corr(f); nf = np.sum(f*f)*dx
        smooth = 2*np.trapezoid(ksm*Cf, z)/nf
        atomsR = -2*sum(w*float(onem(s))*at_C(Cf, s) for s, w in at)/nf
        pred = -sum(Wn*math.cos(T*s) for s, w, Wn in W)
        out.append((T, smooth, atomsR, pred, -sum(Wn for s, w, Wn in W)))
    return out

if __name__ == '__main__':
    print('# c, c_link, T (best alignment in [50,400]), R_smooth/|f|^2, R_atoms/|f|^2, predicted atoms (Kronecker formula at this T), fully aligned value')
    for c, cl in [(10, 10**0.5), (10, 10**0.9), (47, 47**0.5), (47, 47**0.9)]:
        for T, sm, ar, pr, full in run(c, cl):
            print(f'{c:3d} {cl:7.3f} T={T:9.4f}  smooth {sm: .4f}  atoms {ar: .4f}  pred {pr: .4f}  aligned {full: .4f}', flush=True)
