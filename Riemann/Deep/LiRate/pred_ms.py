
import math, numpy as np
def g(a): return np.where(a<=1, a/2, 1-1/(2*a))
for N in (1000,10000,100000,10**6,10**8,10**12):
    t = np.exp(np.linspace(math.log(14.134725), math.log(1e3*math.sqrt(N)), 400000))
    L = np.log(t/(2*math.pi)); rho = L/(2*math.pi); a0 = N/(t**2*L)
    I = 2*np.trapezoid(rho*g(a0)*t, np.log(t))
    print(N, round(I,3), round(I/math.sqrt(N*math.log(N)),4))
