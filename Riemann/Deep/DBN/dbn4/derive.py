
# Independent symbolic derivation of the fixed-kernel polynomial; compare with Planat eq.(15)
import sympy as sp
d,b,y,r,lam=sp.symbols('d b y r lam',real=True)
D,v,B,R=sp.symbols('D v B R',nonnegative=True)
def K(h,dd,bb):  # conjugate-pair kernel K_h(d,b)
    return (h-bb)/(dd**2+(h-bb)**2)+(h+bb)/(dd**2+(h+bb)**2)
# h=1 scaling
lhs=2*K(y,d,b)                 # 2 K_y(d,b)
J=(K(1,d+r,b)+K(1,d-r,b))/2   # J_{1,r}(d,b)
expr=sp.together(lhs/y-lam*J)  # kernel inequality divided by y>0: 2K_y/y - lam*J >= 0
num,den=sp.fraction(sp.factor(expr))
# Planat's formula
C=D+R+1+B; A=C**2+4*D*R-4*B
D4=A**2-16*D*R*C**2; G=A*(C-2*B)-8*D*R*C
Qs=(D+v-B)**2+4*D*B
P=(D+v-B)*D4-lam/2*Qs*G
# check: 2K_y/y == 4(D+v-B)/Qs ... and J == 2G/D4 under d^2=D,b^2=B,y^2=v,r^2=R
sub={d:sp.sqrt(D),b:sp.sqrt(B),y:sp.sqrt(v),r:sp.sqrt(R)}
k1=sp.simplify(sp.together(2*K(y,d,b)/y).subs(sub)-4*(D+v-B)/Qs)
j1=sp.simplify(sp.together(J).subs(sub)-2*G/D4)
print('2K_y/y - 4(D+v-B)/Q* =',k1)
print('J - 2G/D4 =',j1)
# so 2K_y/y - lam J = 4(D+v-B)/Q* - 2 lam G/D4 = (4/(Q* D4)) * [ (D+v-B) D4 - (lam/2) Q* G ] = 4P/(Q* D4)
chk=sp.simplify(4*(D+v-B)/Qs-2*lam*G/D4-4*P/(Qs*D4))
print('identity residual =',chk)
# D4 = product of translated denominators (positive): check D4 == Q*(D+R+2 d r..)... via sympy at numeric points
Dn=lambda u,bb: (u**2+(1-bb)**2)*(u**2+(1+bb)**2)
expr2=sp.expand((Dn(d+r,b)*Dn(d-r,b)).subs({}))
e2=sp.expand(D4.subs({D:d**2,R:r**2,B:b**2})-expr2)
print('D4 - prod of the four translated denominators =',sp.simplify(e2))
e3=sp.expand(Qs.subs({D:d**2,v:y**2,B:b**2})-(d**2+(y-b)**2)*(d**2+(y+b)**2))
print('Q* - (d^2+(y-b)^2)(d^2+(y+b)^2) =',sp.simplify(e3))
