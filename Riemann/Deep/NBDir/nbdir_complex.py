"""nbdir: falsification of DERIVATION §8 (complex chi needs no (chi, chibar)-pairing).
chi mod 5 of order 4 (chi(2) = i, odd). Claim: V_chi = W^{-1}(M-1) D_q F_- I with W = tau(chi)/(i^a sqrt q), and
(D_1 a_chi, Y_s) = (J a_chi, Q_lam Phi_s) = W q^{-1/2} int_1^inf A_chibar(t) t^{-1} conj(phi_s(t/q)) dt = L(s,chi)/s  (s on the line).
Also prints the same integral without W and with A_chi in place of A_chibar (must fail), and L(1-s,chibar)/L(s,chi) vs W^-1 gamma_a(s).
"""
import mpmath as mp
import numpy as np
mp.mp.dps = 20
q, a = 5, 1
chil = [0, 1, 1j, -1j, -1]
pi = mp.pi
def L(s, c=chil): return mp.dirichlet(s, c)
chib = [complex(np.conj(c)) for c in chil]
tau = sum(chil[n] * mp.exp(2j * pi * n / q) for n in range(1, q))
W = tau / (1j**a * mp.sqrt(q))
def gam_a(s): return (q / pi)**(s - 0.5) * mp.gamma((s + a) / 2) / mp.gamma((1 - s + a) / 2)
s0 = mp.mpc(0.3, 2.7)
print("|W| =", mp.nstr(abs(W), 15), " W =", mp.nstr(W, 10))
print("FE: |L(1-s,chibar)/L(s,chi) - W^-1 gamma_a(s)| =", mp.nstr(abs(L(1 - s0, chib) / L(s0) - gam_a(s0) / W), 3),
      " without W:", mp.nstr(abs(L(1 - s0, chib) / L(s0) - gam_a(s0)), 3))
def g(x): return (1 - mp.cos(2 * pi * x)) / (pi * x)
def Etail(w, Y, sgn):
    p = -sgn * 1j
    return p**(1 - w) * mp.gammainc(w - 1, p * Y)
def phi(w, x):
    x = mp.mpf(x); Y = 2 * pi * x
    Ttail = Y**(w - 1) / (1 - w) - (Etail(w, Y, 1) + Etail(w, Y, -1)) / 2
    return g(x) + w / (pi * x) * (2 * pi * x)**(1 - w) * Ttail
Ach = np.cumsum(np.array(chil, complex)); Achb = np.conj(Ach)
s = mp.mpc(0.5, 3.0); T = 600
xs, ws = np.polynomial.legendre.leggauss(12)
I_b = mp.mpc(0); I_c = mp.mpc(0)
for n in range(1, T):
    for xi, wi in zip(xs, ws):
        t = n + 0.5 + 0.5 * xi
        cp = mp.conj(phi(s, t / q)) * 0.5 * wi / t
        I_b += complex(Achb[n % q]) * cp; I_c += complex(Ach[n % q]) * cp
# tail beyond T as in nbdir_burnol.py part C: q^{-1/2}(q/pi) mean(A * (1 - cos(2 pi t/q)))/T, times sqrt(q) for the prefactor below
mA = mp.quad(lambda t: complex(Achb[int(mp.floor(t)) % q]) * (1 - mp.cos(2 * pi * t / q)), mp.linspace(0, q, q + 1)) / q
kap = s / (pi * (1 - s))   # phi_s(x) ~ g(x) + kap/x for a = 1 (DERIVATION §2)
tail = q / pi * mA / T + complex(Achb.mean()) * mp.conj(kap) * q / T
fac = 1 / mp.sqrt(q)
print("(D_1 a, Y_s) with W, A_chibar:", mp.nstr(W * fac * (I_b + tail), 10))
print("   without W               :", mp.nstr(fac * (I_b + tail), 10))
print("   with W, A_chi (no bar)  :", mp.nstr(W * fac * I_c, 10))
print("   target L(s,chi)/s       :", mp.nstr(L(s) / s, 10))
