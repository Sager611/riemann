
import math, random, signal
signal.alarm(100)
random.seed(7)
def reduced_forms(Dabs):
    # primitive reduced positive definite forms of discriminant -Dabs
    F=[]
    a=1
    while 3*a*a<=Dabs:
        for b in range(-a+1,a+1):
            num=b*b+Dabs
            if num%(4*a): continue
            c=num//(4*a)
            if c<a: continue
            if b<0 and a==c: continue
            if math.gcd(math.gcd(a,abs(b)),c)!=1: continue
            F.append((a,b,c))
        a+=1
    return F
def S_values(Dabs, bdict, f):
    nmax=max(bdict)
    res=[]
    for (a,b,c) in reduced_forms(Dabs):
        s=0.0
        ymax=int(math.sqrt(4*a*nmax/Dabs))+2
        for y in range(-ymax,ymax+1):
            disc=b*b*y*y-4*a*(c*y*y-nmax)
            if disc<0: continue
            sq=math.sqrt(disc)
            x0=math.floor((-b*y-sq)/(2*a))-1; x1=math.ceil((-b*y+sq)/(2*a))+1
            for x in range(x0,x1+1):
                m=a*x*x+b*x*y+c*y*y
                if 1<=m<=nmax and m in bdict and math.gcd(m,f)==1: s+=bdict[m]
        res.append(s/2.0)   # w=2 for |D|>4: each ideal counted w times
    return res
def is_fund(D):  # D<0
    Dabs=-D
    if D%4==1 or D%4==-3:
        m=Dabs
        p=2
        while p*p<=m:
            if m%(p*p)==0: return False
            p+=1
        return True
    if D%4==0:
        m=Dabs//4
        if (m%4) not in (1,2): return False
        p=2
        while p*p<=m:
            if m%(p*p)==0: return False
            p+=1
        return True
    return False
def sqfree(n):
    p=2
    while p*p<=n:
        if n%(p*p)==0: return False
        p+=1
    return True
tests=0; worst=1e9
for D in [-23,-31,-47,-56,-71,-79,-84,-104,-131,-143,-191,-239]:
    if not is_fund(D): continue
    for p0 in [5,7,11,13]:
        if (-D)%p0==0: continue
        adm=[n for n in range(1,400) if n%2==1 and sqfree(n) and n%p0!=0 and n%3!=0]
        b={n:random.gauss(0,1) for n in adm}
        SD=S_values(-D,b,1); SL=S_values(-D*p0*p0,b,p0)
        lhs=len(SD)*sum(s*s for s in SD); rhs=len(SL)*sum(s*s for s in SL)
        # projection identity: sum over fibres; check total sums agree
        tests+=1; worst=min(worst, rhs/lhs)
        print(D,p0,'h(D)=',len(SD),'h(p0^2 D)=',len(SL),' sumS_D=%.6f sumS_lift=%.6f'%(sum(SD),sum(SL)),' rhs/lhs=%.4f'%(rhs/lhs))
print('tests',tests,'min rhs/lhs',worst)
