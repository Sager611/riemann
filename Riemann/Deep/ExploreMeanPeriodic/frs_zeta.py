# explore14: FRS boundary function for zeta, Z(s)=Λ(2s)Λ(2s-1)/Λ(s)  (FRS 0803.2821 Sec 5.1, K=Q)
import mpmath as mp, time
mp.mp.dps = 30
pi = mp.pi
def Lam(s): return pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
def Z(s): return Lam(2*s)*Lam(2*s-1)/Lam(s)
t0=time.time()
# 1. Laurent coefficients at s=1/2 via trapezoid on circle
def laurent(k, r=mp.mpf('0.2'), M=96):
    acc = mp.mpc(0)
    for j in range(M):
        th = 2*pi*j/M; u = r*mp.e**(1j*th); s = mp.mpf(1)/2 + u
        acc += Z(s)*u**(k)*u   # (1/2πi)∮ Z u^k du, du = i u dθ
    return acc/M
am2 = laurent(1); am1 = laurent(0)
L12 = Lam(mp.mpf(1)/2)
print('Lambda(1/2) =', mp.nstr(L12, 15))
print('a_{-2} contour =', mp.nstr(am2, 15), ' formula -1/(4Λ(1/2)) =', mp.nstr(-1/(4*L12), 15))
print('a_{-1} contour =', mp.nstr(am1, 5))
a2 = mp.re(-1/(4*L12))
# 2. residues at zeros
B = []
for n in range(1, 11):
    rho = mp.zetazero(n)
    dL = pi**(-rho/2)*mp.gamma(rho/2)*mp.zeta(rho, derivative=1)
    c = Lam(2*rho)*Lam(2*rho-1)/dL
    B.append((mp.im(rho), c))
    print('n=%2d gamma=%s c_rho=%s  (Re/|c|=%s)' % (n, mp.nstr(mp.im(rho),10), mp.nstr(c,6), mp.nstr(mp.re(c)/abs(c),3)))
S1 = sum(abs(c)*g for g,c in B); S0 = sum(abs(c) for g,c in B)
print('sum|b|*gamma (10 zeros) =', mp.nstr(S1,6), ' sum|b| =', mp.nstr(S0,6), ' ratio 2*sum|b|gamma/a_{-2} =', mp.nstr(2*S1/a2,6))
# 3. definition vs Fourier-Carleman expansion
c0 = mp.mpf(2)
def f(y):
    g = lambda t: mp.re(Z(c0+1j*t)*y**(-c0-1j*t))
    return mp.quad(g, mp.linspace(0, 64, 33))/pi
def h_def(x): return f(x) - f(1/x)/x
def h_res(x, N=10):
    # sum of residues of Z(s) x^{-s}: double pole at 1/2 and simple poles at rho, conj(rho)
    lx = mp.log(x)
    tot = x**(-mp.mpf(1)/2)*(am1 - am2*lx)
    for g,c in B[:N]:
        rho = mp.mpf(1)/2 + 1j*g
        tot += c*x**(-rho) + mp.conj(c)*x**(-mp.conj(rho))
    return tot
for x in ['0.02','0.05','0.1','0.3','0.6','0.9','1.5','4']:
    x = mp.mpf(x)
    hd = h_def(x); hr0 = mp.re(x**(-mp.mpf(1)/2)*(-am2*mp.log(x))); hr = h_res(x)
    print('x=%5s h_def=%s  main=%s  def-main=%s  def-(main+10 zeros)=%s' % (mp.nstr(x,3), mp.nstr(hd,18), mp.nstr(hr0,18), mp.nstr(hd-hr0,5), mp.nstr(hd-hr,5)))
print('time %.1f s' % (time.time()-t0))
