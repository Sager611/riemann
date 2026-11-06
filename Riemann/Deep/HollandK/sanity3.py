#!/usr/bin/env python3
# Sanity check 3 (not part of the proof): at n = 1e43 and d on the wedge boundary n^3 log^2(n+2) = 7000 d^5,
# build the finite-free model from the true R_0..R_4 (main term of the saddle, 700 digits), then evaluate the residual
# E_F(z) at points of Omega's boundary (lambda = 1.4, r = 2.6 sqrt(B d)) and compare with the certified bound.
import mpmath as mp
mp.mp.dps = 700
pi = mp.pi
def u0(z):
    u = mp.log(z)/2
    for _ in range(500): u = mp.log((z/u + mp.mpf(9)/4)/pi)/2
    return u
def h(z):
    u = u0(z)
    f = 2*z*mp.log(u) + 9*u/2 - pi*mp.exp(2*u)
    fpp = -2*z/u**2 - 4*pi*mp.exp(2*u)
    return mp.log(2*pi**2) + f + mp.log(2*pi/(-fpp))/2
n = mp.mpf(10)**43
lg = lambda z: h(n+z) - (n+z)*mp.log(4) - mp.loggamma(n+z+mp.mpf(1)/2)    # log gamma(n+z) up to a constant
lg0 = lg(0)
R = [mp.exp(lg(j)-lg0) for j in range(5)]
q = [R[k+1]**2/(R[k]*R[k+2]) for k in range(3)]
L = 2*u0(n); b = n + mp.mpf(1)/2
T = (q[0]-1)/(q[1]-1); AJ = 2*q[0]/(T-q[0]); BJ = 1/(q[0]*(1+1/AJ)-1); D = BJ/2
f = lambda k, U: mp.log((U+k)/(U+k+1))
sol = mp.findroot(lambda a, bb, c: [f(k, a*AJ)-f(k, BJ+bb*n/L)+f(k, D+c*n/L)-f(k, D)-mp.log(q[k]) for k in range(3)],
                  (mp.mpf(6)/5, mp.mpf(1)/3, mp.mpf(1)/24), tol=mp.mpf(10)**(-600))
A = sol[0]*AJ; B = BJ + sol[1]*n/L; C = D + sol[2]*n/L; S = B*R[1]
d = mp.floor((n**3*mp.log(n+2)**2/7000)**(mp.mpf(1)/5))
r = mp.mpf('2.6')*mp.sqrt(B*d); lam = mp.mpf('1.4')
def EF(z):
    return (lg(z) - lg0 - z*mp.log(S) + z*mp.log(A) + z*mp.log(C) - z*mp.log(D)
            + mp.loggamma(B+z) - mp.loggamma(B) - mp.loggamma(A+z) + mp.loggamma(A)
            + mp.loggamma(D+z) - mp.loggamma(D) - mp.loggamma(C+z) + mp.loggamma(C))
print("n = 1e43, L =", mp.nstr(L, 8), " d =", mp.nstr(d, 8), " r =", mp.nstr(r, 6), " solution", [mp.nstr(x, 8) for x in sol])
print("E_F(0..4) =", [mp.nstr(EF(j), 3) for j in range(5)])
pts = [d/2 + 1j*lam*r, -lam*r, d + lam*r, d/2 - 1j*lam*r, d + lam*r*mp.expjpi(mp.mpf(1)/4), -lam*r*mp.expjpi(mp.mpf(-1)/3)]
vals = [abs(EF(z)) for z in pts]
bound_cert = mp.mpf('70.600')*d**mp.mpf(2.5)/(n**mp.mpf(1.5)*L)   # Cfin <= 70.600 (certify.log, L* = 90; n = 1e43 has L > 90)
print("|E_F| on boundary points:", [mp.nstr(v, 5) for v in vals])
print("certified bound Cfin d^{5/2}/(n^{3/2} L) =", mp.nstr(bound_cert, 5), "  (needed <= thetaE, thetaE >= 0.929840)")
