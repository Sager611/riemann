# explore14: DH adversary, part 2: on-line zeros up to height 86, residues, definition vs Fourier-Carleman expansion
import mpmath as mp, time
mp.mp.dps = 30
pi = mp.pi; t0 = time.time()
r5 = mp.sqrt(5); kap = (mp.sqrt(10-2*r5)-2)/(r5-1); A = [1, kap, -kap, -1]
def f(s):  return 5**(-s)*sum(A[a-1]*mp.zeta(s, mp.mpf(a)/5) for a in range(1,5))
def fp(s): return -mp.log(5)*f(s) + 5**(-s)*sum(A[a-1]*mp.zeta(s, mp.mpf(a)/5, 1) for a in range(1,5))
def Lf(s): return (5/pi)**(s/2)*mp.gamma((s+1)/2)*f(s)
def Lam(s): return pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
half = mp.mpf(1)/2
def Zf(s): return s*(1-s)*Lam(2*s)*Lam(2*s-1)/Lf(s)
R = lambda t: mp.re(Lf(half+1j*t))
ts = [mp.mpf(k)/10 for k in range(1, 861)]
vals = [R(t) for t in ts]
zeros = []
for i in range(len(ts)-1):
    if vals[i]*vals[i+1] < 0:
        zeros.append(mp.findroot(R, (ts[i], ts[i+1]), solver='anderson'))
print('on-line DH zeros found in (0,86):', len(zeros), ' first:', [mp.nstr(z,8) for z in zeros[:6]], ' last:', [mp.nstr(z,8) for z in zeros[-3:]])
print('max |Im Lf| on the line (should be 0):', mp.nstr(max(abs(mp.im(Lf(half+1j*t))) for t in ts[::50]), 3))
L12 = Lf(half); a2 = -1/(16*L12)
res = []
for g in zeros:
    rho = half + 1j*g
    dLf = (5/pi)**(rho/2)*mp.gamma((rho+1)/2)*fp(rho)
    c = rho*(1-rho)*Lam(2*rho)*Lam(2*rho-1)/dLf
    res.append((g, c))
print('residues (first 6):', [(mp.nstr(g,7), mp.nstr(c,4)) for g,c in res[:6]])
print('max Re(c)/|c| =', mp.nstr(max(abs(mp.re(c))/abs(c) for g,c in res), 3))
sgn = ''.join('+' if mp.im(c) > 0 else '-' for g,c in res)
print('sign pattern of b_gamma = Im c along the zeros:', sgn)
S1 = sum(abs(c)*g for g,c in res)
print('a_{-2}^DH =', mp.nstr(a2,10), ' 2*sum|b|gamma/|a| (all on-line zeros < 86) =', mp.nstr(2*S1/abs(a2), 6))
# definition vs expansion
c0 = mp.mpf(2)
def F(y): return mp.quad(lambda t: mp.re(Zf(c0+1j*t)*y**(-c0-1j*t)), mp.linspace(0, 64, 33))/pi
def hdef(x): return F(x) - F(1/x)/x
def hexp(x):
    tot = x**(-half)*(-a2*mp.log(x))
    for g,c in res:
        rho = half+1j*g; tot += c*x**(-rho) + mp.conj(c)*x**(-mp.conj(rho))
    return tot
for x in ['0.1','0.5','3']:
    x = mp.mpf(x); hd = hdef(x); he = hexp(x); hm = x**(-half)*(-a2*mp.log(x))
    print('x=%s h_def=%s main=%s def-main=%s def-expansion=%s' % (mp.nstr(x,2), mp.nstr(hd,16), mp.nstr(hm,16), mp.nstr(hd-hm,4), mp.nstr(hd-he,4)))
print('time %.1f s' % (time.time()-t0))
