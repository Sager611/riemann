
# gauss_lemmas.py -- finite checks for "Bounded-Step Walks on Gaussian Primes"
import random, math, itertools, sys

def lemma_differences(trials=3000, seed=1):
    random.seed(seed)
    worst_ratio=1e9; worst_area=1e9; bad=0; cnt=0
    for D2 in [1,2,4,5,9]:      # D^2 values: D=1, sqrt2, 2, sqrt5, 3
        D=math.sqrt(D2)
        steps=[(a,b) for a in range(-3,4) for b in range(-3,4) if 0<a*a+b*b<=D2]
        for t in range(trials//5):
            n=random.choice([5,10,20,40,80,160,320])
            mode=random.random()
            pts=[(0,0)]; seen={(0,0)}
            tries=0
            while len(pts)<n+1 and tries<50*n:
                tries+=1
                x,y=pts[-1]
                if mode<0.3:   # persistent direction (line-like)
                    st=steps[(len(pts)*7)%len(steps)] if random.random()<0.1 else max(steps,key=lambda s:s[0]+0.01*s[1])
                    if random.random()<0.3: st=random.choice(steps)
                else:
                    st=random.choice(steps)
                q=(x+st[0],y+st[1])
                if q in seen: continue
                pts.append(q); seen.add(q)
            if len(pts)<2: continue
            n=len(pts)-1
            E=pts
            # diameter
            R2=-1; P0=P1=None
            for i in range(len(E)):
                for j in range(i+1,len(E)):
                    d=(E[i][0]-E[j][0])**2+(E[i][1]-E[j][1])**2
                    if d>R2: R2=d; P0,P1=E[i],E[j]
            R=math.sqrt(R2)
            ax,ay=P1[0]-P0[0],P1[1]-P0[1]
            W=max([abs(ax*(p[1]-P0[1])-ay*(p[0]-P0[0]))/R for p in E])
            W=max(1.0,W)
            A=R*W
            diffs=set((p[0]-q[0],p[1]-q[1]) for p in E for q in E)
            cnt+=1
            if A>4*math.pi*D2*len(diffs)+1e-9: bad+=1; print('VIOLATION A<=4piD^2|E-E|',D2,n,A,len(diffs))
            if (n+1)>9/math.pi*(R+2/3)*(2*W+2/3)+1e-9: bad+=1; print('VIOLATION n+1 area',D2,n,R,W)
            if A>D2*n*n+1e-9: bad+=1; print('VIOLATION A<=D^2n^2')
            worst_ratio=min(worst_ratio,len(diffs)/(A))
    return cnt,bad,worst_ratio

def cube_isoperimetry():
    bad=0; tot=0
    for r in range(1,5):
        V=list(range(2**r))
        for mask in range(1,2**(2**r)):
            S=[v for v in V if mask>>v&1]
            Sset=set(S)
            edges=sum(1 for v in S for b in range(r) if (v^(1<<b)) not in Sset)
            tot+=1
            if edges < len(S)*math.log2(2**r/len(S))-1e-9: bad+=1
    random.seed(2)
    for r in range(5,11):
        for _ in range(2000):
            k=random.randrange(1,2**r+1)
            if random.random()<0.5:
                S=set(random.sample(range(2**r),k))
            else:   # subcube-like sets (extremal)
                fixed=random.randrange(0,r+1)
                S=set(v for v in range(2**r) if (v & ((1<<fixed)-1))==0)
            edges=sum(1 for v in S for b in range(r) if (v^(1<<b)) not in S)
            tot+=1
            if edges < len(S)*math.log2(2**r/len(S))-1e-9: bad+=1
    return tot,bad

def H(ps):
    return -sum(p*math.log(p) for p in ps if p>0)
def h2(u):
    return 0.0 if u<=0 or u>=1 else -u*math.log(u)-(1-u)*math.log(1-u)

def list_lemma(trials=20000,seed=3):
    random.seed(seed); bad=0
    for _ in range(trials):
        a=random.randrange(2,9); ny=random.randrange(1,6)
        P={(x,y):random.random()**3 for x in range(a) for y in range(ny)}
        Z=sum(P.values()); P={k:v/Z for k,v in P.items()}
        b=random.randrange(1,a+1)
        L={y:set(random.sample(range(a),random.randrange(0,a+1))) for y in range(ny)}
        rho=sum(v for (x,y),v in P.items() if x in L[y] and len(L[y])<=b)
        Py={y:sum(P[(x,y)] for x in range(a)) for y in range(ny)}
        HXgY=sum(Py[y]*H([P[(x,y)]/Py[y] for x in range(a)]) for y in range(ny) if Py[y]>0)
        lhs=math.log(a)-HXgY; rhs=rho*math.log(a/b)-h2(rho)
        if lhs<rhs-1e-9: bad+=1
    return trials,bad

def continuity_lemma(trials=20000,seed=4):
    random.seed(seed); bad=0
    for _ in range(trials):
        ax=random.randrange(1,6); ay=random.randrange(1,6); a=ax*ay
        P=[random.random()**2 for _ in range(a)]; s=sum(P); P=[x/s for x in P]
        Q=[max(0.0,x+random.uniform(-0.2,0.2)*random.random()) for x in P]; s=sum(Q); Q=[x/s for x in Q]
        eps=0.5*sum(abs(x-y) for x,y in zip(P,Q))
        if eps>0.5: continue
        def cond(R):
            tot=0
            for y in range(ay):
                col=[R[x*ay+y] for x in range(ax)]; py=sum(col)
                if py>0: tot+=py*H([c/py for c in col])
            return tot
        lhs=abs(cond(P)-cond(Q)); rhs=2*h2(eps)+2*eps*math.log(a)
        if lhs>rhs+1e-9: bad+=1
    return trials,bad

if __name__=='__main__':
    print('Lemma differences (cnt, violations, min |E-E|/A):',lemma_differences(),flush=True)
    print('cube edge isoperimetry (sets, violations):',cube_isoperimetry(),flush=True)
    print('short-list lemma (trials, violations):',list_lemma(),flush=True)
    print('continuity lemma (trials, violations):',continuity_lemma(),flush=True)

