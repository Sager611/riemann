
import math, random, sys
# ---------- 1. exact r, I=sum i, W=sum(i+j+k) over i+8(j+k)<=D ----------
def sums(D):
    r=I=JK=0
    for i in range(D+1):
        m=(D-i)//8
        c=(m+1)*(m+2)//2          # pairs (j,k), j+k<=m
        s=m*(m+1)*(m+2)//3        # sum of (j+k) over those pairs
        r+=c; I+=i*c; JK+=s
    return r,I,I+JK
print("1) exact sums vs claimed leading terms")
for D in (10**3,10**4,10**5):
    r,I,W=sums(D)
    print(f"  D={D}: r/(D^3/384)={r/(D**3/384):.6f}  W/(5D^4/6144)={W/(5*D**4/6144):.6f}  I/(D^4/1536)={I/(D**4/1536):.6f}  (W-5I/4)/D^3={(W-1.25*I)/D**3:.6f}")
# E_p = sum floor(i/p) vs D^4/(1536 p): check error/D^3 bounded uniformly in p<=D
D=4000
cnt=[((D-i)//8+1)*((D-i)//8+2)//2 for i in range(D+1)]
worst=0
for p in (11,13,101,997,3989):
    Ep=sum(c*(i//p) for i,c in enumerate(cnt))
    err=(Ep-D**4/(1536*p))/D**3
    worst=max(worst,abs(err))
print(f"  E_p: max |E_p - D^4/(1536p)|/D^3 over sample p (D={D}) = {worst:.4f}")
# ---------- 2. biquadratic arithmetic: elements as (c0,c1,c2,c3) on basis 1,a,b,ab ; a^2=q, b^2=2 ----------
def mul(x,y,q,mod=None):
    x0,x1,x2,x3=x; y0,y1,y2,y3=y
    z0=x0*y0+q*x1*y1+2*x2*y2+2*q*x3*y3
    z1=x0*y1+x1*y0+2*(x2*y3+x3*y2)
    z2=x0*y2+x2*y0+q*(x1*y3+x3*y1)
    z3=x0*y3+x3*y0+x1*y2+x2*y1
    z=(z0,z1,z2,z3)
    return tuple(t%mod for t in z) if mod else z
def pw(x,e,q,mod):
    res=(1,0,0,0); x=tuple(t%mod for t in x)
    while e:
        if e&1: res=mul(res,x,q,mod)
        x=mul(x,x,q,mod); e>>=1
    return res
def leg(a,p):
    t=pow(a%p,(p-1)//2,p); return -1 if t==p-1 else t
def kron(q,p): return leg(q,p)  # odd p not dividing q
sig=lambda x:(x[0],-x[1],x[2],-x[3])
tau=lambda x:(x[0],x[1],-x[2],-x[3])
def primes(n):
    s=bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]
print("2) Frobenius congruences theta^p = sigma(theta) / sigma tau(theta) mod p for chi(p)=-1")
random.seed(1)
bad=0; tested=0
for q in (-3,-4,5,-7,-8,12,13,-15,-23,-163,21,-1155):
    for p in primes(400):
        if p<=8 or (2*q)%p==0: continue
        if kron(q,p)!=-1: continue
        g=sig if leg(2,p)==1 else (lambda x:sig(tau(x)))
        for _ in range(3):
            th=tuple(random.randrange(0,50) for _ in range(4))
            lhs=pw(th,p,q,p); rhs=tuple(t%p for t in g(th))
            tested+=1; bad+= (lhs!=rhs)
print(f"  tested {tested} cases, failures {bad}")
# ---------- 3. det(gamma1,gamma2,gamma3) and kernel vector, exact in K ----------
print("3) det and kernel (exact in Q(a,b))")
for q in (5,-23,12,-8,8):
    one=(1,0,0,0); A=(0,1,0,0); Bv=(0,0,1,0); AB=(0,0,0,1)
    neg=lambda x:tuple(-t for t in x)
    g=[(one,one,one),(A,neg(A),neg(A)),(Bv,Bv,neg(Bv)),(AB,neg(AB),AB)]
    v=[AB,Bv,neg(A),neg(one)]
    add=lambda x,y:tuple(s+t for s,t in zip(x,y))
    comp=[]
    for row in range(3):
        s=(0,0,0,0)
        for j in range(4): s=add(s,mul(v[j],g[j][row],q))
        comp.append(s)
    # det of 3x3 with entries g[col][row]
    M=[[g[c][r] for c in range(3)] for r in range(3)]
    def m(x,y): return mul(x,y,q)
    def sub(x,y): return tuple(s-t for s,t in zip(x,y))
    det=add(add(m(M[0][0],sub(m(M[1][1],M[2][2]),m(M[1][2],M[2][1]))),
                neg(m(M[0][1],sub(m(M[1][0],M[2][2]),m(M[1][2],M[2][0]))))),
                m(M[0][2],sub(m(M[1][0],M[2][1]),m(M[1][1],M[2][0]))))
    print(f"  q={q}: pi(v)={comp}  det={det} (expect 4ab=(0,0,0,4))")
# ---------- 4. Lemma 3 / box constant ----------
ok=all(1+math.sqrt(2)+(1+math.sqrt(2))*math.sqrt(Q)<=8*math.sqrt(Q) for Q in range(3,10**5))
print("4) 1+sqrt2+(1+sqrt2)sqrt|q| <= 8 sqrt|q| for 3<=|q|<1e5:",ok)
# ---------- 5. (PS) => count condition, rigorous lower bound E_p >= I/p - r, theta(Q)<=1.04Q ----------
print("5) slack in (PS) => (count-condition) using exact r,I,W and E_p >= I/p - r")
for (N,qa) in ((10**3,10),(10**4,10),(10**6,10**3)):
    D=int(0.5*((N-1)/3)**(4/3)); Q=D/math.log(N)**2; B=8*N*math.sqrt(qa)
    if D>2*10**6:
        r=D**3/384; I=D**4/1536; W=5*D**4/6144+0.02*D**3
    else:
        r,I,W=sums(D)
    lhs=I*1.3*math.log(B)-r*1.04*Q
    rhs=W*math.log(B)+r*math.log(32*r*math.sqrt(qa))
    print(f"  N={N},|q|={qa}: D={D}, lhs/rhs={lhs/rhs:.4f}")

