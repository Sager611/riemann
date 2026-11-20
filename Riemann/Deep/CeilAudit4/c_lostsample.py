# C: morirate lost-sample law bookkeeping (floats)
import math
from scipy.integrate import quad
g2 = lambda y: math.log((2+math.sqrt(4-y*y))/y)
I, err = quad(g2, 0, 2, limit=200)
print(f"mean of g2 on (0,2) = {I/2:.12f}  (pi/2 = {math.pi/2:.12f})")
# zeros below 4 pi mu: N(T) = (T/2pi) log(T/(2 pi e)) + 7/8
for mu in [10, 45, 1e3, 1e6, 1e12]:
    T = 4*math.pi*mu
    N = T/(2*math.pi)*math.log(T/(2*math.pi*math.e)) + 7/8
    for eps in [0.01, 0.1]:
        lost = eps/math.log(T)*N          # merge a proportion eps/log T of the zeros (o(1) for every asymptotic statistic)
        # zero density in y = gamma/(2 pi mu) on (0,2): ~ (mu/1) log(mu y) ... weight g2 by local density
        dens = lambda y: mu*math.log(max(mu*y, 1.0000001))   # dN/dy up to constants (dN/dgamma = log(gamma/2pi)/2pi, dgamma = 2 pi mu dy)
        num, _ = quad(lambda y: 4*g2(y)*dens(y), 1e-9, 2, limit=400)
        den, _ = quad(dens, 1e-9, 2, limit=400)
        shift = lost*num/den               # nats lost = sum of 4 g2 over lost samples, density-weighted
        print(f"mu={mu:.0e} eps={eps}: N(4 pi mu)={N:.4g}, lost={lost:.4g} (= {lost/mu:.4f} mu), extra -log lambda = {shift:.4g} = {shift/(4*math.pi*mu):.4f} x 4 pi mu")
print("delta* ~ e^{-3.24 mu + 3.6} (MoriRate l.79): mu = 10 ->", math.exp(-3.24*10+3.6), "; mean spacing at 4 pi mu=125.7:", 2*math.pi/math.log(4*math.pi*10/(2*math.pi)))
