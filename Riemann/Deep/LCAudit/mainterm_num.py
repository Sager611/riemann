#!/usr/bin/env python3
"""lcaudit: numeric sanity check of the leading-order evaluation of Wu's main term (exact finite sums at y = 1e4..1e6).
M(A,B) = L(1+al+be,chi0) S_q(al,be) + e^{-(A+B)} L(1-al-be,chi0) S_q(-be,-al),  al = A/L, be = B/L, L = log(y)/theta,
S_q(al,be) = sum_{h,k<=y, (hk,q)=1} mu(h)mu(k)P[h]P[k](h,k)^{1+al+be} h^{-1-be} k^{-1-al}  (Wu l.400 with the exponent corrected).
Leading order (mainterm_sym.py): c(A,B). Expect M -> c with error O(1/log y), and M_5/M_1 -> 1 (q/phi(q) cancels).
Run: uv run --no-project --with numpy --with mpmath python3 mainterm_num.py
"""
import numpy as np, mpmath as mp, math, sys
mp.mp.dps = 30
def mu_and_primes(y):
    isp = np.ones(y + 1, bool); isp[:2] = False
    for p in range(2, int(y**0.5) + 1):
        if isp[p]: isp[p*p::p] = False
    primes = np.nonzero(isp)[0]
    mu = np.ones(y + 1, np.int8); mu[0] = 0
    for p in primes:
        mu[p::p] *= -1
        if p*p <= y: mu[p*p::p*p] = 0
    return mu, primes
def S(y, primes, base, al, be, q):
    """sum_d J_{1+al+be}(d) U_d(be) U_d(al), U_d(x) = sum_m f_x(dm), f_x(n) = base(n) n^{-1-x}; base = mu P chi0."""
    n = np.arange(y + 1, dtype=float); n[0] = 1.0
    def U(x):
        g = base * n**(-1.0 - x)
        for p in primes:                      # sum over multiples (squarefree support): g[n] += g[n p]
            k = y // p
            g[1:k+1] += g[p:k*p+1:p]
        return g
    Ub, Ua = U(be), U(al)
    J = np.ones(y + 1)
    sexp = 1.0 + al + be
    for p in primes:
        J[p::p] *= (p**sexp - 1.0)
    w = base != 0                              # d squarefree and coprime to q (else U_d = 0)
    return float(np.sum(J[w] * Ub[w] * Ua[w]))
def Lchi0(z, q):
    v = mp.zeta(1 + z)
    for p in [p for p in (2, 3, 5, 7) if q % p == 0]:
        v *= (1 - mp.mpf(p)**(-1 - z))
    return float(v)
def c_formula(A, B, th, P, dP):
    f = lambda e1, e2: float(mp.quad(lambda x: (dP(x) + e1*th*P(x))*(dP(x) + e2*th*P(x)), [0, 1]))
    return (f(A, B) - math.exp(-(A + B))*f(-A, -B))/((A + B)*th)
out = []
def say(m):
    print(m, flush=True); out.append(m)
th = 0.5
rho = 1.2036; Rh = 2*math.pi*rho
say("theta = %.2f, L = log y / theta; columns: y, q, A=B, P, M (exact sums), c(A,B) (leading order), M/c - 1" % th)
for y in (10**4, 10**5, 10**6):
    mu, primes = mu_and_primes(y)
    ly = math.log(y); L = ly/th
    nn = np.arange(y + 1, dtype=float); nn[0] = 1.0
    for Pname in ('x', 'sinh'):
        for A in (-1.0, 1.0, Rh):
            if Pname == 'x':
                Pf = lambda x: x; dPf = lambda x: 1
            else:
                b = th*abs(A)
                Pf = lambda x, b=b: mp.sinh(b*x)/mp.sinh(b); dPf = lambda x, b=b: b*mp.cosh(b*x)/mp.sinh(b)
            Pn = np.array([0.0] + [0.0]*0) # placeholder
            xs = np.log(y/nn)/ly
            Pv = xs if Pname == 'x' else np.sinh(th*abs(A)*xs)/math.sinh(th*abs(A))
            Pv[0] = 0.0
            cval = c_formula(A, A, th, Pf, dPf)
            Ms = {}
            for q in (1, 5):
                base = mu.astype(float) * Pv
                if q == 5: base[::5] = 0.0
                base[0] = 0.0
                al = A/L
                T1 = Lchi0(2*al, q) * S(y, primes, base, al, al, q)
                T2 = math.exp(-2*A) * Lchi0(-2*al, q) * S(y, primes, base, -al, -al, q)
                Ms[q] = T1 + T2
                say("y=1e%d q=%d A=B=%+.4f P=%-4s M=%.6f c=%.6f M/c-1=%+.4f  (T1=%.6f, T2=%.3e)" % (round(math.log10(y)), q, A, Pname, Ms[q], cval, Ms[q]/cval - 1, T1, T2))
            say("   ratio M_5/M_1 - 1 = %+.5f" % (Ms[5]/Ms[1] - 1))
open('logs/mainterm_num.log', 'w').write('\n'.join(out) + '\n')
