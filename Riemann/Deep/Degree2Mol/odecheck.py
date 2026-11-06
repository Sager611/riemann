
import mpmath as mp
mp.mp.dps = 30
def kappa_paper(theta):
    R = mp.sqrt(mp.mpf(3)/5)/theta
    # solve S'' + tanh(t) S' - S + 1/(1+e^{2t}) = 0, S(0)=1/2, S(R)=0 by linear shooting
    def f(t, y): return [y[1], -mp.tanh(t)*y[1] + y[0] - 1/(1+mp.e**(2*t))]
    def shoot(s):
        sol = mp.odefun(f, 0, [mp.mpf(1)/2, s]); return sol(R)
    # linear in s: S(R) = A + B s
    y0 = shoot(0); y1 = shoot(1)
    A = y0[0]; B = y1[0] - y0[0]
    s = -A/B
    sol = mp.odefun(f, 0, [mp.mpf(1)/2, s])
    I1 = mp.quad(lambda t: mp.e**(-t)*sol(t)[0], [0, R/2, R])
    c = mp.mpf(1)/2 + (1/mp.sqrt(15))*(-1 + 2*mp.e**R*(1 - mp.e**(-R) - s - I1))
    # direct functional value check: I(S) = int e^t(S^2+S'^2) + e^{-t}((1-S)^2+S'^2)
    I = mp.quad(lambda t: mp.e**t*(sol(t)[0]**2 + sol(t)[1]**2) + mp.e**(-t)*((1-sol(t)[0])**2 + sol(t)[1]**2), [0, R/2, R])
    c2 = mp.mpf(1)/2 - mp.sqrt(mp.mpf(3)/5)/3 + 2/mp.sqrt(15)*mp.e**R*I
    return R, s, 1 - mp.log(c)/R, 1 - mp.log(c2)/R
for th in ['0.5', '0.25', '0.125', '0.6666666666666666']:
    R, s, k1, k2 = kappa_paper(mp.mpf(th))
    print(f"theta={th}: R={mp.nstr(R,8)} S'(0)={mp.nstr(s,10)} kappa(simplcPRQ)={mp.nstr(k1,10)} kappa(direct I)={mp.nstr(k2,10)}")
