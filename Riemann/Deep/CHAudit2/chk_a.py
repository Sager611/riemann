# chaudit2: numerics for the technical lemmas of 2511.14736v1 (barbar, artanor, zetinvbound's Stirling input)
import mpmath as mp
import numpy as np
mp.mp.dps = 30

print("== barbar: beta = 6.4 constant ==")
b = mp.mpf('6.4')
cst = mp.ei(b) - mp.e**b*(b-1)/(b*(b-2))
print("Ei(6.4) - e^6.4*5.4/(6.4*4.4) =", mp.nstr(cst, 12), "<= 2.525:", cst <= mp.mpf('2.525'))
print("beta/(beta-2) =", mp.nstr(b/(b-2), 12), "= 16/11:", mp.nstr(mp.mpf(16)/11, 12))
# optimal beta for the u0-free part (minimise Ei(beta) - e^beta (beta-1)/(beta(beta-2))) is irrelevant; check the inequality g(x)>=g(beta)
def g(x, be):
    return mp.e**x/x*(1+be/((be-2)*x)) - mp.ei(x)
mn = min(g(mp.mpf(x)/100, b) for x in range(1, 5000))
print("min_{0<x<=50} g(x) - g(6.4) =", mp.nstr(mn - g(b, b), 6), "(>= 0 expected)")

def kappa(a):
    a = mp.mpf(a)
    return 2*mp.e**(a/2) - a*mp.ei(a/2) + mp.mpf('2.525')*a

# direct test of barbar: int_{u0}^{-1/2} e^{-a u}/u^2 du <= 16/11 e^{-a u0}/(a u0^2) + kappa_a
worst = 0
rng = np.random.default_rng(1)
cases = [(a, u0) for a in [1e-6, 1e-3, 0.05, 1/6, 0.5, 1, 3, 10] for u0 in [-0.5, -0.5001, -0.6, -1, -2, -5, -20, -80]]
cases += [(float(10**rng.uniform(-6, 1.3)), -float(10**rng.uniform(np.log10(0.5), 2))) for _ in range(300)]
for a, u0 in cases:
    a = mp.mpf(a); u0 = mp.mpf(u0)
    lhs = mp.quad(lambda u: mp.e**(-a*u)/u**2, [u0, (u0-0.5)/2, mp.mpf(-0.5)]) if u0 < -0.5 else mp.mpf(0)
    rhs = mp.mpf(16)/11*mp.e**(-a*u0)/(a*u0**2) + kappa(a)
    r = lhs/rhs
    if r > worst:
        worst = r; wc = (a, u0)
print("barbar: max LHS/RHS over", len(cases), "cases =", mp.nstr(worst, 8), "at alpha, u0 =", [mp.nstr(t, 6) for t in wc])

print("== artanor constants ==")
def cs(a):
    a = mp.mpf(a); k = kappa(a); p = mp.pi
    c0 = mp.e**(a/2)/2 + mp.e**(-a/2)/(8*p) + k/(16*p)
    c1 = 2*mp.e**(a/2) + (a/(8*p) + mp.mpf(3)/11)*mp.e**(-a/2) + (a/(16*p) + mp.mpf(3)/22)*k
    c2 = mp.mpf(16)/11*(a/(16*p) + mp.mpf(3)/22)
    return c0, c1, c2
for a in ['0.5', '1/6']:
    c0, c1, c2 = cs(mp.mpf(eval(a)))
    print("alpha =", a, " c0 =", mp.nstr(c0, 10), " c1 =", mp.nstr(c1, 10), " c2 =", mp.nstr(c2, 10))
print("claims: c0(1/2)<16/21=%.7f c1(1/2)<17/5 c2(1/2)<2/9=%.7f; c0(1/6)<=2/3=%.7f c1(1/6)<=2.85" % (16/21, 2/9, 2/3))
# monotonicity on (0,1/2]
prev = None; mono = True
for j in range(1, 501):
    v = cs(mp.mpf(j)/1000)
    if prev is not None and any(v[i] < prev[i] for i in range(3)):
        mono = False; print("non-monotone at", j/1000)
    prev = v
print("c0,c1,c2 increasing on (0,1/2] (step 0.001):", mono, "; inf near 0: c0 ->", mp.nstr(cs(mp.mpf('1e-12'))[0], 8), "c1 ->", mp.nstr(cs(mp.mpf('1e-12'))[1], 8))

print("== zetinvbound: Stirling remainder after 1/(12z), Re z >= 0 ==")
# R2(z) = loggamma(z) - [(z-1/2)log z - z + log sqrt(2 pi) + 1/(12 z)]; claim |R2| <= sqrt2/(180|z|^3) (GR 8.344); Stieltjes gives 1/(90|z|^3) at arg = pi/2
mx = 0; arg_at = None
for th_i in range(0, 91, 3):
    th = mp.pi/2*th_i/90
    for lr in np.linspace(-1.3, 2.0, 120):
        r = mp.mpf(10)**lr
        z = r*mp.e**(1j*th)
        R2 = mp.loggamma(z) - ((z-0.5)*mp.log(z) - z + mp.log(mp.sqrt(2*mp.pi)) + 1/(12*z))
        q = abs(R2)*abs(z)**3
        if q > mx:
            mx = q; arg_at = (float(th_i)*1.0, float(r))
print("sup |R2(z)| |z|^3 over grid =", mp.nstr(mx, 8), "at (arg in degrees, |z|) =", arg_at, "; sqrt2/180 =", mp.nstr(mp.sqrt(2)/180, 8), "; 1/90 =", mp.nstr(mp.mpf(1)/90, 8), "; 1/360 =", mp.nstr(mp.mpf(1)/360, 8))
# what the proof needs: Re z/(12|z|^2) - |R2(z)| > 0 for Re z >= 1/2
need = 1e9
for th_i in range(0, 90):
    th = mp.pi/2*th_i/90
    for lr in np.linspace(np.log10(0.5), 3, 80):
        r = mp.mpf(10)**lr
        z = r*mp.e**(1j*th)
        if z.real < 0.5:
            continue
        R2 = mp.loggamma(z) - ((z-0.5)*mp.log(z) - z + mp.log(mp.sqrt(2*mp.pi)) + 1/(12*z))
        val = (z.real/(12*abs(z)**2) + R2.real)/(z.real/(12*abs(z)**2))
        need = min(need, float(val))
print("min over Re z >= 1/2 of (Re z/(12|z|^2) + Re R2)/(Re z/(12|z|^2)) =", need, "(> 0 needed)")
