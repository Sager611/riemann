
# A: Basu-Midha mechanism vs the Riemann kernel.
#  R_n(z) = E[(z+iX)^n] = sum_l (-1)^l C(n,2l) m_l z^(n-2l),  X ~ Phi(|u|)du normalized  (z-variable Jensen polys of Xi, reversed)
#  (1) moments growth m_l^(1/2l): Gaussian ~ sqrt(2l/e), compact support -> bounded, Riemann -> ~ (1/4)log l (unbounded, sub-Gaussian)
#  (2) constant-coefficient admissible flows from delta_0 (d_p b_k = -(4 a k^2 + 2 b k) b_k + 2 g k(2k-1) b_{k-1}) give
#      Gaussians (a=0) or compactly supported laws (a>0): b_k^(1/2k) bounded
#  (3) R_n real-rooted for n <= NMAX (sanity), and the dBN-induced generator d_t R_n = -(R_{n+2} - 2 z R_{n+1} + z^2 R_n) leaves degree n
import mpmath as mp, sys
mp.mp.dps = 40
def Phi(u):
    return mp.fsum((2*mp.pi**2*n**4*mp.e**(9*u) - 3*mp.pi*n**2*mp.e**(5*u))*mp.e**(-mp.pi*n**2*mp.e**(4*u)) for n in range(1, 8))
L = 30
pts = [mp.mpf(x)/8 for x in range(0, 17)]
b = [2*mp.quad(lambda u: u**(2*k)*Phi(u), pts) for k in range(L+1)]
m = [bk/b[0] for bk in b]
print('b0 =', mp.nstr(b[0], 12), '(should be Xi(0)/8*2? H_0(0)=Xi(0)/8 -> b0=2*H0(0)=Xi(0)/4 =', mp.nstr(mp.mpf('0.497120778188314')/4, 12), ')')
print(' l   m_l^(1/2l)    gaussian-equiv sqrt((2l-1)!!)^(1/l)*sigma   ratio m_l/m_{l-1}')
sig = mp.sqrt(m[1])
for l in [1,2,3,4,5,8,10,15,20,25,30]:
    g = (mp.fac2(2*l-1))**(mp.mpf(1)/(2*l))*sig
    print('%3d  %.6f   %.6f   %.6f' % (l, float(m[l]**(mp.mpf(1)/(2*l))), float(g), float(m[l]/m[l-1])))
# (2) constant-coefficient flow from delta_0
def flow(a, bb, g, p, K=30, steps=4000):
    bk = [mp.mpf(1)] + [mp.mpf(0)]*K
    h = mp.mpf(p)/steps
    for _ in range(steps):
        nb = bk[:]
        for k in range(1, K+1):
            nb[k] = bk[k] + h*(-(4*a*k*k + 2*bb*k)*bk[k] + 2*g*k*(2*k-1)*bk[k-1])
        bk = nb
    return bk
mp.mp.dps = 20
for (a, bb, g, p) in [(0.05, 0, 1, 2.0), (0.01, 0.1, 1, 3.0), (0.2, -0.3, 1, 1.0)]:
    bk = flow(a, bb, g, p)
    print('flow a=%g b=%g g=%g p=%g: b_k^(1/2k) at k=5,10,20,30:' % (a,bb,g,p), [round(float(bk[k]**(mp.mpf(1)/(2*k))),4) for k in (5,10,20,30)], ' bound 2sqrt(g/a)=%.4f' % (2*(g/a)**0.5))
mp.mp.dps = 60
NMAX = 30
mmp = [mp.mpf(x) for x in m]
worst = 0
for n in range(2, NMAX+1):
    coeffs = [0]*(n+1)  # highest degree first for polyroots
    for l in range(0, n//2+1):
        coeffs[2*l] = (-1)**l*mp.binomial(n, 2*l)*mmp[l]
    r = mp.polyroots(coeffs, maxsteps=400, extraprec=400)
    im = max(abs(mp.im(x)) for x in r)
    worst = max(worst, im)
    if n in (2,5,10,20,30): print('n=%d max|Im root| = %s  largest root = %s  (Hermite-scale 2sqrt(n m1) = %s)' % (n, mp.nstr(im,3), mp.nstr(max(mp.re(x) for x in r),6), mp.nstr(2*mp.sqrt(n*mmp[1]),6)))
print('max |Im| over n<=%d: %s' % (NMAX, mp.nstr(worst, 3)))
