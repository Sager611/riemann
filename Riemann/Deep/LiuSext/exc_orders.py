# exc_orders.py — exceptional-row exponent bookkeeping of Liu 2610.26835v2 §12.8 (l. 7146-7327) for residue symbols of order n.
# (1) derive the exceptional excess from count + raw volume - child allowance and compare with the printed F1, F2 (n = 6);
# (2) local F2 table >= kappa*b2 for n = 6 and n = 3; (3) centering budget for the balanced input A = M.
import sympy as sp
A,M,m,q,c,d,R,E,w,wo,l,g,g2,t2,V,f,b2,p2,Bc,D0,n = sp.symbols('A M m q c d R E w w_o ell g g2 t2 V f b2 p2 B_c D0 n')
qt = q+R+E                              # q~ = q+R+E (l. 6656)
K = 2*A-c+R+E-m-D0                      # K = K0 + d - D0, K0 = 2A-c-d+R+E-m (l. 6656)
J = d-c+(2*A-c-d+R+E-m-K)-2*w+wo        # l. 6772
a0 = A-c-w                              # l. 6766
mp = 2*a0-K-g-g2-V                      # m' (l. 6923)
Mp = M+J-g-g2+t2                        # M' (l. 6925)
Dp = b2-p2+w+Bc+l                       # Delta' (l. 6968)
count = (mp-f)/n                        # #exceptional rows <= Z^{(m'-f)/n} (l. 7189, n = 6)
vol = a0-b2                             # raw paired volume (l. 7196)
Ex = sp.expand((count+vol-(Mp+Dp)).subs(M,m+q))
F1 = c/n+(1-1/n)*D0+2*w/n+qt/n+wo+Bc-(1-1/n)*g+l
F2 = 2*b2-(1-1/n)*g2-p2+t2+V/n+f/n
claim = sp.expand((A-(1-1/n)*M-F1-F2).subs(M,m+q))
print('(1) Ex - [A-(1-1/n)M-F1-F2] =', sp.simplify(Ex-claim))
F1p = c/6+5*D0/6+w/3+qt/6+wo+Bc-5*g/6+l; F2p = 2*b2-5*g2/6-p2+t2+V/6+f/6
print('    n=6 matches printed F1,F2 (l. 7203-7207):', sp.simplify(F1.subs(n,6)-F1p)==0 and sp.simplify(F2.subs(n,6)-F2p)==0)
# (2) local F2 at one common prime; cases (g2,p2,t2,V,f,b2) as in l. 7230-7236; 'n | i' replaces '6 | i'
def table(N):
    rows=[]
    for i in range(1,61):
        if i%N: rows.append(('unit',i,dict(g2=i,p2=1,t2=1,V=0,f=0,b2=i)))
        if i%N: rows.append(('nonunit',i,dict(g2=i,p2=1,t2=0,V=1,f=(2 if N==6 else 1) if i==1 else 0,b2=i)))   # i=1: h' needs p^4 (n=6) or p^1 (n=3)
        if i%N==0: rows.append(('n|i',i,dict(g2=i,p2=1,t2=0,V=0,f=0,b2=i)))
        for j0 in range(N,i,N): rows.append(('i>j0,n|j0',(i,j0),dict(g2=j0,p2=1,t2=0,V=0,f=0,b2=sp.Rational(i+j0,2))))
    return rows
for N in (6,3):
    worst=None
    for name,i,s in table(N):
        val=F2.subs(n,N).subs({g2:s['g2'],p2:s['p2'],t2:s['t2'],V:s['V'],f:s['f'],b2:s['b2']})
        r=sp.nsimplify(val/s['b2'])
        if worst is None or r<worst[0]: worst=(r,name,i,val)
    print(f'(2) n={N}: min F2/b2 = {worst[0]} at {worst[1]} i={worst[2]} (F2={worst[3]})')
# F1 lower bound: J>=0 (g=J_+): F1 = c+2w+...; J<0 (g=0): min over D0 in [0,c] of c/n+(1-1/n)D0+(3c-5D0)_+/6 (+R/n), w=0, B_c=(3c-5d-R)_+/6 with d<=D0
x=sp.symbols('x',nonnegative=True)
for N in (6,3):
    h=lambda D: sp.Rational(1,N)+(1-sp.Rational(1,N))*D+sp.Max(0,3-5*D)/6
    vals=[h(sp.Rational(k,600)) for k in range(0,601)]
    print(f'    n={N}: F1/c lower bound (J<0, w=0, sextic B_c) = {min(vals)}')
# (3) centering at A = M: needed saving M/n; available min(kappa,1)*L; comparison total L+(M-(A-L)) must be < (1-1/n)M
for N,kap in ((6,sp.Rational(2,3)),(3,sp.Rational(11,15)),(3,sp.Integer(1))):
    Lmax=(1-sp.Rational(1,N))/2; need=sp.Rational(1,N); have=min(kap,1)*Lmax
    print(f'(3) n={N}, kappa={kap}: need {need}M, best saving {have}M (L<{Lmax}M) -> slack {have-need}M; paper L=M/4 gives {min(kap,1)/4-need}M' if N==6 else f'(3) n={N}, kappa={kap}: need {need}M, best saving {have}M (L<{Lmax}M) -> slack {have-need}M')
