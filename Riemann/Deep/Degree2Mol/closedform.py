
# Evaluate 2508.11108's closed form (eqn: SRT), l.476-505, at R = sqrt(3/5)/theta and compare with the shooting solution of (mainODEbdry).
import mpmath as mp
mp.mp.dps = 30
phi = (1 + mp.sqrt(5))/2
Fp  = lambda t: mp.hyp2f1(mp.mpf(1)/2, phi, mp.mpf(1)/2 + phi, -mp.e**(2*t))
Fm  = lambda t: mp.hyp2f1(mp.mpf(1)/2, -1/phi, mp.mpf(1)/2 - 1/phi, -mp.e**(2*t))
F1p = lambda u: mp.hyp2f1(mp.mpf(1)/2, 1 + phi, mp.mpf(1)/2 + phi, -mp.e**(2*u))
F1m = lambda u: mp.hyp2f1(mp.mpf(3)/2, -1/phi, mp.mpf(1)/2 - 1/phi, -mp.e**(2*u))
g1 = lambda t: mp.e**(-t/phi)*Fm(t)
g2 = lambda t: mp.e**(phi*t)*Fp(t)
W  = lambda u: mp.e**u*(2*phi*Fm(u)*F1p(u) - Fp(u)*F1m(u))
v1 = lambda u: g2(u)/(W(u)*(1 + mp.e**(2*u)))
v2 = lambda u: -g1(u)/(W(u)*(1 + mp.e**(2*u)))
f  = lambda t: g1(t) - g1(0)/g2(0)*g2(t)
g0 = lambda t: g2(t)/(2*g2(0))
w1 = lambda t: mp.quad(v1, [0, t]); w2 = lambda t: mp.quad(v2, [0, t])
def C1(R):
    num = -2*Fm(R)*Fp(0)*w1(R) - mp.e**(mp.sqrt(5)*R)*Fp(R)*(1 + 2*Fp(0)*w2(R))
    den = 2*Fm(R)*Fp(0) - 2*mp.e**(mp.sqrt(5)*R)*Fm(0)*Fp(R)
    return num/den
def S(t, R, c1): return c1*f(t) + g0(t) + g1(t)*w1(t) + g2(t)*w2(t)
# Wronskian check: W should equal g1 g2' - g1' g2
t0 = mp.mpf('0.7'); print("Wronskian check:", mp.nstr(W(t0), 12), mp.nstr(g1(t0)*mp.diff(g2, t0) - mp.diff(g1, t0)*g2(t0), 12))
def odeS(R):
    fr = lambda t, y: [y[1], -mp.tanh(t)*y[1] + y[0] - 1/(1+mp.e**(2*t))]
    y0 = mp.odefun(fr, 0, [mp.mpf(1)/2, 0])(R)[0]; y1 = mp.odefun(fr, 0, [mp.mpf(1)/2, 1])(R)[0]
    s = -y0/(y1 - y0); return s, mp.odefun(fr, 0, [mp.mpf(1)/2, s])
for th in ['1/2', '1/8']:
    R = mp.sqrt(mp.mpf(3)/5)/mp.mpf(eval(th)); c1 = C1(R)
    s_ode, sol = odeS(R)
    dS0 = mp.diff(lambda t: S(t, R, c1), 0)
    print(f"theta={th} R={mp.nstr(R,8)}: closed-form C1={mp.nstr(c1,10)}, S(R)={mp.nstr(S(R,R,c1),6)}, S'(0)={mp.nstr(dS0,10)} vs shooting S'(0)={mp.nstr(s_ode,10)}; S(R/2): {mp.nstr(S(R/2,R,c1),10)} vs {mp.nstr(sol(R/2)[0],10)}")
    I1 = mp.quad(lambda t: mp.e**(-t)*S(t, R, c1), [0, R])
    c = mp.mpf(1)/2 + (1/mp.sqrt(15))*(-1 + 2*mp.e**R*(1 - mp.e**(-R) - dS0 - I1))
    print(f"   kappa from closed form + (simplcPRQ): {mp.nstr(1 - mp.log(c)/R, 10)}")
