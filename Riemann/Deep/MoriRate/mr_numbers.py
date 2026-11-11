
import math
eps = {5:7.447489666e-18, 10:9.133835606e-44, 16:1.512800554e-75, 23:6.02e-113, 25:1.17e-123, 45:1.1506e-231, 47:1.874e-242}
U = {5:1.222e-17,6:9.575e-23,7:1.006e-27,8:5.211e-33,9:3.229e-38,10:1.944e-43,12:6.271e-54,14:1.642e-64,16:3.743e-75,18:6.775e-86,20:1.657e-96}
print('c  -ln eps - 4pi c   +5ln c  | fit(-5lnc-15.86) | ThmB floor -(23 ln c + ln 1.362e9)')
for c in sorted(eps):
    v = -math.log(eps[c]) - 4*math.pi*c
    print(c, round(v,2), round(v+5*math.log(c),2), round(-5*math.log(c)-15.86,2), round(-(23*math.log(c)+math.log(1.362e9)),2))
def lg10(x): return math.log10(x)
for mu in (5,47,50,100):
    l=math.log(mu)
    b = lg10(1.362e9)+23*lg10(mu)-4*math.pi*mu/math.log(10)
    ga = lg10(1.437e23)+4.5*lg10(mu)+4*lg10(l)-4*math.pi*mu/math.log(10)
    gb = lg10(1.24e23)+5*lg10(mu)+5*lg10(l)-4*math.pi*mu/math.log(10)
    fit = 15.86/math.log(10)+5*lg10(mu)-4*math.pi*mu/math.log(10)
    print('mu',mu,'ThmB 10^%.3f'%b,'G(a) 10^%.3f'%ga,'G(b) 10^%.3f'%gb,'fit 10^%.3f'%fit)
# Table 1
for mu in (50,1e3,1e6,1e12):
    l=math.log(mu)
    print('T1 mu=%g Ga/E=%.3g Gb/E=%.3g'%(mu, 1.437e23*mu**4.5*l**4/(5.046e17*mu**8*l**3), 1.24e23*mu**5*l**5/(5.046e17*mu**8*l**3)))
print('6621(2pi)^4.5/0.1089 =', 6621*(2*math.pi)**4.5/0.1089, ' x6.048e14=%.5g x5.2e14=%.5g'%(237522697*6.048e14, 237522697*5.2e14))
# integral of g2 times -log|y| on (-2,2) (Mori: dual value 2pi)
import mpmath as mp
g2=lambda y: mp.log((2+mp.sqrt(4-y*y))/y)
I = 2*mp.quad(lambda y: g2(y)*(-mp.log(y)), [0,1,2])
print('int g2 (-log|y|) =', I, ' 2pi =', 2*mp.pi, ' int_0^2 g2 =', mp.quad(g2,[0,2]))
