
# Symbolic check of the two-shift mixed-probe kernel polynomial
import sympy as sp
d,b,y,r1,r2,l1,l2=sp.symbols('d b y r1 r2 l1 l2',real=True)
D,v,B,R=sp.symbols('D v B R',nonnegative=True)
def K(h,dd,bb): return (h-bb)/(dd**2+(h-bb)**2)+(h+bb)/(dd**2+(h+bb)**2)
def J(r): return (K(1,d+r,b)+K(1,d-r,b))/2
def parts(Rv):
    C=D+Rv+1+B; A=C**2+4*D*Rv-4*B
    return A**2-16*D*Rv*C**2, A*(C-2*B)-8*D*Rv*C
R1,R2=sp.symbols('R1 R2',positive=True)
D41,G1=parts(R1); D42,G2=parts(R2)
Qs=(D+v-B)**2+4*D*B
Pmix=(D+v-B)*D41*D42-Qs/2*(l1*G1*D42+l2*G2*D41)
# claim: 2K_y/y - l1 J_{r1} - l2 J_{r2} = 4 Pmix/(Q* D41 D42)
lhs=sp.together(2*K(y,d,b)/y-l1*J(r1)-l2*J(r2))
sub={d:sp.sqrt(D),b:sp.sqrt(B),y:sp.sqrt(v),r1:sp.sqrt(R1),r2:sp.sqrt(R2)}
res=sp.simplify(lhs.subs(sub)-4*Pmix/(Qs*D41*D42))
print('mixed identity residual =',res)
print('degree in D:',sp.Poly(sp.expand(Pmix),D).degree(),' in B:',sp.Poly(sp.expand(Pmix),B).degree(),' in v:',sp.Poly(sp.expand(Pmix),v).degree())
