# Is there a map H(p,q) -> G(>= dc^2) of points7.txt / points8.txt sending each H vertex to a point whose nearest
# family point carries its label (then the NOTE section 3 lemma gives the lower bound)?  Standard library only.
import os, sys, itertools
from fractions import Fraction
sys.setrecursionlimit(10000)
def read(path):
    pts=[]
    for line in open(path):
        line=line.split('#')[0].strip()
        if not line: continue
        pts.append(tuple(Fraction(s) for s in line.replace(',',' ').split()))
    return pts
def sqd(x,y): return sum((a-b)**2 for a,b in zip(x,y))
def fl(p): return [float(c) for c in p]
base=os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../Gale7/')
def label_family(S, m1):
    # P: second part zero; Q: first part zero
    P=[k for k,s in enumerate(S) if all(c==0 for c in s[m1:])]
    Q=[k for k,s in enumerate(S) if all(c==0 for c in s[:m1])]
    import math
    def unit(v):
        nr=math.sqrt(sum(c*c for c in v)); return [c/nr for c in v], nr
    U=[unit(fl(S[k][:m1]))[0] for k in P]; W=[unit(fl(S[k][m1:]))[0] for k in Q]
    lab={}
    for i,k in enumerate(P): lab[k]=('P',i,None)
    for j,k in enumerate(Q): lab[k]=('Q',None,j)
    for k,s in enumerate(S):
        if k in lab: continue
        x,nx=unit(fl(s[:m1])); y,ny=unit(fl(s[m1:]))
        # best i by |dot|
        di=[sum(a*b for a,b in zip(x,u)) for u in U]; dj=[sum(a*b for a,b in zip(y,w)) for w in W]
        i=max(range(len(U)), key=lambda t: abs(di[t])); j=max(range(len(W)), key=lambda t: abs(dj[t]))
        si=1 if di[i]>0 else -1; sj=1 if dj[j]>0 else -1
        if si>0 and sj<0: t='A'
        elif si<0 and sj>0: t='B'
        else: t='C'
        lab[k]=(t,i,j)
    return lab, len(P), len(Q)
EDGES = [('P','P','n','-'), ('Q','Q','-','n'), ('P','A','n','-'), ('Q','B','-','n'), ('A','A','n','e'), ('A','A','n','n'),
         ('B','B','e','n'), ('B','B','n','n'), ('P','C','n','-'), ('A','C','n','e'), ('A','C','n','n'), ('Q','C','-','n'),
         ('B','C','e','n'), ('B','C','n','n'), ('C','C','n','n'), ('P','B','e','-'), ('Q','A','-','e')]
ORB='PQABC'
def rel(x,y): return '-' if x is None or y is None else ('e' if x==y else 'n')
def isH(a,b):
    X,Y=a,b
    if ORB.index(X[0])>ORB.index(Y[0]): X,Y=Y,X
    return (X[0],Y[0],rel(X[1],Y[1]),rel(X[2],Y[2])) in EDGES
def analyse(fn, symfn, m1):
    P=read(base+fn); N=len(P)
    S=read(base+symfn) if symfn else P
    lab,p,q=label_family(S,m1)
    from collections import Counter
    print(fn, 'labels', Counter(v[0] for v in lab.values()), 'p,q', p, q)
    D2=max(sqd(P[i],P[j]) for i in range(N) for j in range(i))
    # nearest S for each P
    near={}
    for i,x in enumerate(P):
        ds=sorted((sqd(x,s),k) for k,s in enumerate(S))
        near[i]=ds[0][1]
    import math
    nd=sorted(math.sqrt(float(sqd(P[i],S[near[i]]))) for i in range(N))
    print(' nearest dists: median %.4f, 90%% %.4f, max %.4f' % (nd[N//2], nd[int(N*0.9)], nd[-1]))
    # candidates for each H vertex (S index): P points whose nearest is it
    cand={k:[i for i in range(N) if near[i]==k] for k in range(len(S))}
    Hv=list(range(len(S)))
    Hedges=[(a,b) for a in Hv for b in Hv if a<b and isH(lab[a],lab[b])]
    print(' H edges', len(Hedges))
    # for given threshold s, test if choice exists; first: for each H edge, min over best choices
    d2=[[sqd(P[i],P[j]) for j in range(N)] for i in range(N)]
    return P,S,lab,cand,Hedges,d2,D2
for fn,sym,m1,dcnum,dcden in [('points7.txt','points7_sym69.txt',3,8900856812318928,10**16), ('points8.txt',None,4,9023994457279939,10**16)]:
    P,S,lab,cand,Hedges,d2,D2=analyse(fn,sym,m1)
    dc2=Fraction(dcnum,dcden)
    print(' D2*1e16', D2*10**16, ' dc2', float(dc2), 'ratio', float(dc2/D2))
    N=len(P)
    # count edges in G(>=dc2)
    E=sum(1 for i in range(N) for j in range(i) if d2[i][j]>=dc2)
    print(' edges of G(>=dc2):', E, ' edges >dc2:', sum(1 for i in range(N) for j in range(i) if d2[i][j]>dc2))
    # CSP: phi(k) in cand[k] (allow all candidates); check H edges
    order=sorted(range(len(S)), key=lambda k: len(cand[k]))
    adjH={k:set() for k in range(len(S))}
    for a,b in Hedges: adjH[a].add(b); adjH[b].add(a)
    phi={}
    nodes=[0]
    def ok(k,i):
        for w in adjH[k]:
            if w in phi and d2[i][phi[w]]<dc2: return False
        return True
    def rec(t):
        nodes[0]+=1
        if t==len(order): return True
        k=order[t]
        for i in cand[k]:
            if ok(k,i):
                phi[k]=i
                if rec(t+1): return True
                del phi[k]
        return False
    res=rec(0)
    print(' homomorphism H -> G(>=dc2) with nearest-candidates:', res, 'nodes', nodes[0])
    if not res:
        # which H edges fail for all candidate pairs
        bad=[(a,b) for a,b in Hedges if all(d2[i][j]<dc2 for i in cand[a] for j in cand[b])]
        print(' H edges impossible under candidates:', len(bad), [ (lab[a],lab[b]) for a,b in bad[:10]])
    else:
        mn=min(d2[phi[a]][phi[b]] for a,b in Hedges)
        print(' min sqd over mapped H edges', float(mn), ' equals dc2?', mn==dc2)
