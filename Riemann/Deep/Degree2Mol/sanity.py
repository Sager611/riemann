
import mpmath as mp
mp.mp.dps = 30
def c_PQ(Qc, R, th, P='x'):
    # Q polynomial coefficients in y (low->high); P=x: A=1, B=1/3
    Q = lambda y: mp.polyval(Qc[::-1], y)
    dQc = [k*Qc[k] for k in range(1, len(Qc))]
    dQ = lambda y: mp.polyval(dQc[::-1], y) if dQc else mp.mpf(0)
    W0 = mp.quad(lambda y: mp.e**(2*R*y)*Q(y)**2, [0, 1])
    W1 = mp.quad(lambda y: mp.e**(2*R*y)*(R*Q(y) + dQ(y))**2, [0, 1])
    return (1 + mp.e**(2*R)*Q(1)**2)/2 + W0/th + th*W1/3
# Levinson: theta=1/2, P=x, Q=1-y
f = lambda R: -(1 - mp.log(c_PQ([1, -1], R, mp.mpf(1)/2))/R)
Rb = mp.findroot(lambda R: mp.diff(f, R), 1.3)
print("Levinson (theta=1/2,P=x,Q=1-y): R*=", mp.nstr(Rb, 8), "kappa=", mp.nstr(-f(Rb), 8))
# Conrey's Q (2508.11108 l.974): Q(y)=0.492+0.602(1-2y)-0.08(1-2y)^3-0.06(1-2y)^5+0.046(1-2y)^7 ; with optimal P at theta=4/7
import numpy as np
z = np.polynomial.Polynomial([1, -2])
Qp = 0.492 + 0.602*z - 0.08*z**3 - 0.06*z**5 + 0.046*z**7
Qc = [mp.mpf(x) for x in Qp.coef]
def kap_optP(Qc, R, th):
    Q = lambda y: mp.polyval(Qc[::-1], y)
    dQc = [k*Qc[k] for k in range(1, len(Qc))]
    dQ = lambda y: mp.polyval(dQc[::-1], y)
    W0 = mp.quad(lambda y: mp.e**(2*R*y)*Q(y)**2, [0, 1]); W1 = mp.quad(lambda y: mp.e**(2*R*y)*(R*Q(y) + dQ(y))**2, [0, 1])
    c = (1 + mp.e**(2*R)*Q(1)**2)/2 + mp.sqrt(W0*W1)/mp.tanh(th*mp.sqrt(W1/W0))
    return 1 - mp.log(c)/R
for R in [1.2, 1.3, 1.4]:
    print("Conrey Q, optimal P, theta=4/7, R=", R, "kappa=", mp.nstr(kap_optP(Qc, mp.mpf(R), mp.mpf(4)/7), 8), " Q(0)=", float(Qp(0)), "Q(1)=", float(Qp(1)))
