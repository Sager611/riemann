
# Exact check of the weighted Golod-Shafarevich margin in "A surface counterexample to Shafarevich holomorphic convexity"
# (Prop. infinite-quotient): generators P(t)=(1+t)G(t), involution relators <= tG(t), G(t)=3t(1+t)(1+t^2)^2,
# two allowances of 1/100 (deep monodromy, boundary powers), t=49/200 (needs 4t<1).
from fractions import Fraction as Fr
t=Fr(49,200)
G=3*t*(1+t)*(1+t*t)**2
P=(1+t)*G
print('P(1)=',(1+1)*3*1*(1+1)*(1+1)**2,'generators (expect 48);  G(1)=',3*2*4,'(expect 24 involution relators)')
print('4t<1:',4*t<1)
print('G(49/200) =',G,'~',float(G))
print('paper lower bound 3*49*249*53^2/(200^2*50^2) =',Fr(3*49*249*53**2,200**2*50**2),' < G:',Fr(3*49*249*53**2,200**2*50**2)<G,' > 257/250:',Fr(3*49*249*53**2,200**2*50**2)>Fr(257,250))
val=1-P+t*G+Fr(1,100)+Fr(1,100)
print('GS expression 1-P+tG+1/50 =',float(val),' < -1/125:',val<Fr(-1,125))
# weight distribution of the 48 generators from (Z/4)^2 Jennings-type filtration: 3 orbits x t(1+t)^2(1+t^2)^2
import sympy as sp
x=sp.symbols('x')
print('P(x) expanded:',sp.expand((1+x)*3*x*(1+x)*(1+x**2)**2))
