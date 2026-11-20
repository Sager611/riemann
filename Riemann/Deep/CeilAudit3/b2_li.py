
import mpmath as mp, math
mp.mp.dps=30
T=mp.mpf(3000175332800); a=mp.mpf("0.8733458561067185768")-mp.mpf(1)/2
EG=mp.euler
def main(n): return n/2*(mp.log(n)-1-mp.log(2*mp.pi)+EG)
def q(n,t,d=a): return 4*mp.sinh(n*d/(2*t**2))**2+4*n*d**2/t**3
def allzeros(c):   # method model: every zero above T at beta=B, density (1/2pi)log(t/2pi)
    n=c*T**2
    f=lambda t: q(n,t)*mp.log(t/(2*mp.pi))/(2*mp.pi)
    return mp.quad(f,[T,T*1.002,T*1.01,T*1.05,T*1.3,T*3,T*100])
def onezero(c):
    n=c*T**2; return 2*q(n,T)   # one quartet at height T (2 conj pairs -> factor; generous)
def cmax(fun):
    lo,hi=50,300
    for _ in range(60):
        mid=(lo+hi)/2
        if fun(mid)<main(mid*T**2): lo=mid
        else: hi=mid
    return lo
c1=cmax(allzeros); print("all zeros at beta=B (method model): c =",round(c1,3))
c2=cmax(onezero); print("one worst quartet at T: c =",round(c2,3))
# KLN-shaped density: #zeros with beta>=sigma up to height ~T(1+eps) <= A T^{8(1-s)/3} (log T)^{5-2s}
LT=mp.log(T)
def D(s,A): return A*T**(8*(1-s)/3)*LT**(5-2*s)
def density(c,A):
    n=c*T**2; B=mp.mpf(1)/2+a
    # worst case: allowed zeros placed at height T, beta distributed by -dD/dsigma, capped by the method's own all-zero sum
    g=lambda s: q(n,T,s-mp.mpf(1)/2)*(-mp.diff(lambda z: D(z,A),s))
    val=mp.quad(g,[B-mp.mpf("0.3"),B-mp.mpf("0.05"),B])+D(B,A)*q(n,T)
    return min(val,allzeros(c))
for A in [1,10,100]:
    print("density A=",A,": D(B)=",mp.nstr(D(mp.mpf(1)/2+a,A),4)," c =",round(cmax(lambda c: density(c,A)),3))
print("zeros in [T,1.0126T]:",mp.nstr(T*mp.mpf("0.0126")*mp.log(T/(2*mp.pi))/(2*mp.pi),4))
