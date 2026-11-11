# Reproduce Lamzouri 2611.14486 l.452-471: Selberg minorant value at lambda=0.7892, theta=4/7-1e-5 (Arb)
from flint import arb, ctx
ctx.prec = 200
pi = arb.pi()
th = arb(4)/7 - arb(1)/10**5
def J(lam):
    lam = arb(lam)
    return lam + lam*th**2 - arb(2)/3*lam**2*th**3 - th*(2*pi*lam*th).cos()/(2*pi**2) + (2*pi*lam*th).sin()/(4*pi**3*lam)
print("J(0.7892) =", J("0.7892"))
# smallest lambda with J > 1 (bisection, Selberg family)
lo, hi = arb("0.7"), arb("0.7892")
for _ in range(60):
    mid = (lo+hi)/2
    if J(mid) > 1: hi = mid
    else: lo = mid
print("Selberg family threshold lambda* ~", hi.mid().str(10))
# cross-check by numerical quadrature of the closed form
import math
def Jq(lam, n=200000):
    t = 4/7-1e-5; s = 0.0
    for i in range(n):
        a = (i+0.5)*t/n
        s += (1-lam*a+math.sin(2*math.pi*lam*a)/(2*math.pi))*a
    return lam + 2*lam*s*t/n
print("quadrature J(0.7892) =", Jq(0.7892))
