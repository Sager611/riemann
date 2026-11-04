
# Small-case check of the label-type design (Prop. "Balanced types and parity") of
# "A Torsion-Free Group Algebra That Is Not Directly Finite" (Oct 4 2026), with PG(2,4) in place of PG(2,128)
# (q=4 has q+1=5 = 1 mod 4 like q+1=129). Also checks the parity criterion's degree argument on the product graph abstractly.
import itertools, random
# GF(4) = {0,1,w,w^2} encoded 0..3 with w^2 = w+1
def add(a,b): return a^b
def mul(a,b):
    if a==0 or b==0: return 0
    lg={1:0,2:1,3:2}; ex=[1,2,3]
    return ex[(lg[a]+lg[b])%3]
F=range(4)
vecs=[v for v in itertools.product(F,repeat=3) if any(v)]
def normalize(v):
    for c in v:
        if c:
            inv={1:1,2:3,3:2}[c]; return tuple(mul(inv,x) for x in v)
pts=sorted({normalize(v) for v in vecs}); pid={p:i for i,p in enumerate(pts)}
lines=[]
for l in pts:
    S=frozenset(pid[p] for p in pts if (mul(l[0],p[0])^mul(l[1],p[1])^mul(l[2],p[2]))==0)
    lines.append(S)
q=4; v=len(pts); assert v==q*q+q+1 and all(len(L)==q+1 for L in lines)
assert all(len(L1&L2)==1 for L1,L2 in itertools.combinations(lines,2))
# Fano plane on extra letters E = 0..6 ; lines of PG(2,2)
fano=[frozenset(s) for s in [(0,1,2),(0,3,4),(0,5,6),(1,3,5),(1,4,6),(2,3,6),(2,4,5)]]
assert all(len(a&b)==1 for a,b in itertools.combinations(fano,2))
D=[frozenset(range(7))-L for L in fano]
print('|D|=',[len(d) for d in D],' pairwise intersections:',sorted({len(a&b) for a,b in itertools.combinations(D,2)}))
print('each extra letter in',[sum(e in d for d in D) for e in range(7)],'members; each pair in',sorted({sum(a in d and b in d for d in D) for a,b in itertools.combinations(range(7),2)}))
m=5; assert m%4==1 and (q+1)%4==1
am=((q+1)*m-1)//4; bm=(q+1)*(m-1)//4
# build vertices: A classes indexed by lines with m vertices, B with m-1
A=[(li,k) for li in range(v) for k in range(m)]; B=[(li,k) for li in range(v) for k in range(m-1)]
random.seed(1)
def assign(V,counts,full_vertex=None):
    extra={x:frozenset() for x in V}
    if full_vertex is not None: extra[full_vertex]=frozenset(range(7))
    free=[x for x in V if x!=full_vertex]; random.shuffle(free)
    it=iter(free)
    for di,d in enumerate(D):
        for _ in range(counts): extra[next(it)]=d
    return extra
xA=A[0]
exA=assign(A,am,xA); exB=assign(B,bm)
S={('A',x):(frozenset(('P',p) for p in lines[x[0]])|frozenset(('E',e) for e in exA[x])) for x in A}
S.update({('B',y):(frozenset(('P',p) for p in lines[y[0]])|frozenset(('E',e) for e in exB[y])) for y in B})
T=[('P',p) for p in range(v)]+[('E',e) for e in range(7)]
n=v*m; p=(q+1)/v
bal=all(sum(t in S[('A',x)] for x in A)==round(n*p) and sum(t in S[('B',y)] for y in B)==round((n-v)*p) for t in T)
viol=0
for x,x2 in itertools.product(A,A):
    k=len(S[('A',x)]&S[('A',x2)])
    want=0 if (x==xA and x2==xA) else 1
    if k%2!=want: viol+=1
for x,y in itertools.product(A,B):
    if len(S[('A',x)]&S[('B',y)])%2!=1: viol+=1
print(f'q={q}: |A|={len(A)} |B|={len(B)} label balance (|Y_t| = |Y|p for all t): {bal}; parity violations over all pairs: {viol}')
print('sizes |S_x|:',sorted({len(s) for s in S.values()}),' |S_xA| =',len(S[('A',xA)]))
