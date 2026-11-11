
# T1: Chandrasekhar-Darboux <-> Maass weight shift. Weight-k archimedean constant-term integral.
from mpmath import mp, mpf, mpc, quad, gamma, sqrt, pi, atan2, exp, inf, fabs, zeta
mp.dps = 30
def I(k, s):
    # I_k(s) = int_R ((u+i)/|u+i|)^(-k) (1+u^2)^(-s) du
    f = lambda u: exp(-1j*k*atan2(1, u)) * (1+u*u)**(-s)
    return quad(f, [-inf, -10, -1, 0, 1, 10, inf])
def closed(k, s):
    return (-1)**(k//2) * sqrt(pi) * gamma(s-mpf(1)/2) * gamma(s) / (gamma(s+k/2) * gamma(s-k/2))
print("k  s  numeric I_k(s)  closed form  rel.err")
for s in [mpc(1.3, 0.7), mpc(0.9, 5.0), mpc(2.1, -3.0)]:
    for k in [0, 2, 4, 6]:
        a = I(k, s); b = closed(k, s)
        print(k, mp.nstr(s, 4), mp.nstr(a, 12), mp.nstr(b, 12), mp.nstr(abs(a-b)/abs(b), 3))
print("Chandrasekhar factor C_k(s) = I_{k+2}/I_k vs (1-s+k/2)/(s+k/2):")
for s in [mpc(1.3, 0.7), mpc(0.9, 5.0)]:
    for k in [0, 2, 4]:
        r = I(k+2, s)/I(k, s); c = (1 - s + k/2)/(s + k/2)
        print(k, mp.nstr(s, 4), mp.nstr(r, 12), mp.nstr(c, 12))
print("|C_k(1/2+it)| on the unitary axis:", [mp.nstr(abs((1-(0.5+1j*t)+k/2)/((0.5+1j*t)+k/2)), 15) for k in (0,2,4) for t in (0.3, 7.0)])
# weight-k scattering 'matrix' phi_k = I_k * zeta(2s-1)/zeta(2s): poles at rho/2 for every k (common arithmetic factor)
rho1 = mpc(0.5, 14.134725141734693790)
for k in (0, 2, 4):
    s = rho1/2 + mpc(1e-8, 0)
    print("k=%d |phi_k| near rho_1/2: %s" % (k, mp.nstr(abs(closed(k, s)*zeta(2*s-1)/zeta(2*s)), 5)))
