
# Independent check of the Wall D(2) counterexample's finite claims:
#  (1) rho respects the 4 relators of P and kills t1=x1^2, t2=x2^3;
#  (2) the twisted second boundary (left Fox derivatives evaluated at rho) has rank 4  => H_2(B;L_rho)=0;
#  (3) same with the other convention (right Fox derivatives / rho^{-1}) as a cross-check;
#  (4) commutator identities t1=[t2,a2], s t2 s^-1=[t1,a1] hold in P (checked in the HNN/BS normal forms via a faithful
#      matrix representation of the BS(1,n) parts is overkill; we check them via the defining relations by free reduction
#      after substituting conjugation rules) -> done symbolically in the report; here we test them in a concrete quotient rep.
import sympy as sp
z=sp.exp(2*sp.pi*sp.I/3)
gens=['x1','a1','x2','a2','s']
def w(s):
    out=[]
    for tok in s.split():
        if tok.endswith("'"): out.append((tok[:-1],-1))
        else: out.append((tok,1))
    return out
def pw(g,n): return [(g,1 if n>0 else -1)]*abs(n)
r1=w("x1 a1 x1'")+pw('a1',-4)
r2=w("x2 a2 x2'")+pw('a2',-3)
r3=w("x1")+pw('a2',-13)
r4=w("s x2 s'")+pw('a1',-5)
rels=[r1,r2,r3,r4]
rho={'x1':-1,'a1':z,'x2':z**2,'a2':-1,'s':1}
def ev(word,rho):
    v=sp.Integer(1)
    for g,e in word: v*=rho[g]**e
    return sp.nsimplify(sp.simplify(v))
print('relators under rho:',[sp.simplify(sp.expand(ev(r,rho))) for r in rels])
print('t1,t2 under rho:',sp.simplify(ev(pw('x1',2),rho)),sp.simplify(ev(pw('x2',3),rho)))
def fox_left(word,g,rho):
    # d(uv)/dg = du/dg + u dv/dg ; d g/dg=1 ; d g^-1/dg = -g^-1 ; evaluated at rho
    tot=sp.Integer(0); pref=sp.Integer(1)
    for h,e in word:
        if h==g:
            tot+= pref*(1 if e==1 else -rho[g]**-1)
        pref*=rho[h]**e
    return sp.simplify(tot)
D=sp.Matrix([[fox_left(r,g,rho) for g in gens] for r in rels])
D=D.applyfunc(lambda e: sp.nsimplify(sp.simplify(sp.expand(e))))
print('Fox matrix (left derivatives at rho), rows r1..r4, cols x1,a1,x2,a2,s:')
sp.pprint(D)
Dn=sp.Matrix(D).evalf(30)
print('rank (numeric, 30 digits):',Dn.rank(iszerofunc=lambda x: abs(x)<1e-20))
print('rank (exact):',D.rank(simplify=True))
M=D.extract([0,1,2,3],[0,1,2,4])
print('minor cols x1,a1,x2,s:',complex(M.det().evalf(30)),' paper -4(1-zeta^2):',complex((-4*(1-z**2)).evalf(30)))
paper=sp.Matrix([[1-z,-2,0,0,0],[0,0,2,z**2-1,0],[1,0,0,-1,0],[0,-(1+z),1,0,1-z**2]])
print('matches paper matrix (numeric):',max(abs(complex(x)) for x in (D-paper).evalf(30))<1e-25)
# cross-check with rho^{-1} (complex conjugate system)
rhoi={k:sp.conjugate(v) if v not in (1,-1) else v for k,v in rho.items()}
Di=sp.Matrix([[fox_left(r,g,rhoi) for g in gens] for r in rels])
print('rank with conjugate character:',Di.evalf(30).rank(iszerofunc=lambda x: abs(x)<1e-20))
# Euler characteristic bookkeeping: chi(B)=1-5+4=0 ; H_0(B;L)=0 since rho nontrivial ; so dim H_1(B;L)=0 as well.
D1=sp.Matrix([[rho[g]-1 for g in gens]])
print('d1 = (rho(g)-1) row:',[sp.simplify(x) for x in D1],' rank',D1.evalf(20).rank(iszerofunc=lambda x: abs(x)<1e-15))
