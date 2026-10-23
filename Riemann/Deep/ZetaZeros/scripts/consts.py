
# consts.py -- independent recomputation of the constants of arXiv 2610.08965 (Knausgard)
import json, math
from fractions import Fraction as Fr
import numpy as np
import mpmath as mp
mp.mp.dps = 50
A = 'src/anc/'
cert = {'first': json.load(open(A+'certificate/inputs.json')),
        'second': json.load(open(A+'certificate_line/inputs.json')),
        'third': json.load(open(A+'certificate_third/inputs.json'))}

def Cf(a):  # C(f) = 4/3 + 1/2 sum a^2 + 1/pi^2 sum((-1)^j a_j/j^2 - a_j^2/(4 j^2))
    s1 = sum(Fr(x)**2 for x in a)/2
    s2 = sum(Fr((-1)**j)*a[j-1]/(j*j) - a[j-1]**2/(4*j*j) for j in range(1,len(a)+1))
    return mp.mpf(4)/3 + mp.mpf(s1.numerator)/s1.denominator + (mp.mpf(s2.numerator)/s2.denominator)/mp.pi**2

def Cf_numeric(a, n=200001):  # direct quadrature check of int f^2 + iint |s-t| f f
    u = np.linspace(-0.5,0.5,n); f = 1+sum(a[j]*np.cos(2*np.pi*(j+1)*u) for j in range(len(a)))
    h = u[1]-u[0]; w = np.full(n,h); w[0]=w[-1]=h/2
    I1 = np.sum(w*f*f)
    # iint |s-t| f(s) f(t) = int f(s) [ int |s-t| f(t) dt ] ds ; inner via cumulative sums
    F0 = np.cumsum(w*f); F1 = np.cumsum(w*f*u)
    tot0, tot1 = F0[-1], F1[-1]
    inner = u*F0 - F1 + (tot1 - F1) - u*(tot0 - F0)
    return I1 + np.sum(w*f*inner)

def Kfun(a, t):
    t = np.asarray(t, float)
    s = np.sinc(t)
    for j,aj in enumerate(a,1):
        s = s + aj/2*(np.sinc(t-j)+np.sinc(t+j))
    return s

out = {}
for name,d in cert.items():
    a = [Fr(x, d['denominator']) for x in d['window_coefficients']]
    af = [float(x) for x in a]
    C = Cf(a); cap = Fr(d['energy_cap']); capm = mp.mpf(cap.numerator)/cap.denominator
    th = float(Fr(d['separation'])); g0 = float(Fr(d['majorant_gamma']))
    u = np.linspace(0,0.5,2000001)
    f = 1+sum(af[j]*np.cos(2*np.pi*(j+1)*u) for j in range(len(af)))
    g = g0*(np.sinc(th*(u-0.5))+np.sinc(th*(u+0.5)))**2
    Lam = 2*g0/th*(1+np.sinc(th))
    t = np.linspace(0, th, 400001); K2 = Kfun(af,t)**2
    out[name] = dict(C=mp.nstr(C,20), C_quad=float(Cf_numeric(af)), cap=mp.nstr(capm,20), cap_minus_C=mp.nstr(capm-C,6),
        f_min=float(f.min()), f_half=float(f[-1]), g_minus_f_min=float((g-f).min()), argmin_u=float(u[np.argmin(g-f)]),
        Lambda=Lam, minK2_on_0_theta=float(K2.min()), K_theta=float(Kfun(af,th)))
for k,v in out.items():
    print(k, v)
# certificate formulas
d=cert['first']; D=d['denominator']
x1 = Fr(sum(d['h'][0]),D); x2 = Fr(sum(d['h'][1]),D); b = 6*Fr(d['alpha'],D); C0 = Fr(d['energy_cap']); c = Fr(67,20)
q = (3-C0+x2-x1-b)/(2-2*x1+x2)
mu=2-2*x1+x2; nu=2*c-3+x1-x2; rh=6*c-7-c*c+x1-2*x2; rk=4*c-2-c*c-2*x1
print('first: x1',x1,float(x1),'x2',x2,float(x2),'b',float(b),'q',q,float(q),'== paper', q==Fr(1645064,1960733))
print('  mu,nu,rh,rk', float(mu),float(nu),float(rh),float(rk))
d=cert['second']; D=d['denominator']
x = Fr(sum(d['h'][0]),D); b = 7*Fr(d['alpha'],D); C0=Fr(d['energy_cap'])
s0 = (2-C0-b)/(1-x); print('second: x',x,'b',b,'s0',s0,float(s0),'== paper', s0==Fr(1669159,2478195))

