
# Independent check of the explicit seed complex in lean/ComparatorChallenges/EilenbergGanea.lean
# (paper: height kernel of the RAAG on a flag triangulation of the presentation complex <x,y | x^2=y^5, x^2=(xy^-1)^3>).
import itertools, sympy
from sympy.combinatorics.free_groups import free_group
from sympy.combinatorics.fp_groups import FpGroup
seedWord={False:[(False,True),(False,True)]+[(True,False)]*5,
          True:[(False,True),(False,True),(True,True),(False,False),(True,True),(False,False),(True,True),(False,False)]}
def word_str(r):
    return ''.join(('x' if not g else 'y')+('' if s else '^-1') for g,s in seedWord[r])
print('relator words:',word_str(False),'|',word_str(True))
# group-theoretic sanity: abelianisation, order of the presented group
F,x,y=free_group('x,y')
def tow(r):
    w=F.identity
    for g,s in seedWord[r]: w=w*((x if not g else y)**(1 if s else -1))
    return w
G=FpGroup(F,[tow(False),tow(True)])
print('order of <x,y | relators> =',G.order())
M=sympy.Matrix([[sum((1 if s else -1) for g,s in seedWord[r] if g==gg) for gg in (False,True)] for r in (False,True)])
print('abelianisation relation matrix',M.tolist(),'det',M.det())
# build the complex exactly as in the Lean file
rose=[None]+[(g,i) for g in (False,True) for i in range(2)]
SV=[('r',v) for v in rose]+[('c',b) for b in (False,True)]
CE=[(g,k) for g in (False,True) for k in range(3)]
def cstart(e): return [None,(e[0],0),(e[0],1)][e[1]]
def cend(e):   return [(e[0],0),(e[0],1),None][e[1]]
L={r:3*len(seedWord[r]) for r in (False,True)}
SEC=[(r,i) for r in (False,True) for i in range(L[r])]
def letter(s): return seedWord[s[0]][s[1]//3]
def secedge(s):
    g,sg=letter(s); return (g, s[1]%3 if sg else 2-s[1]%3)
def secstart(s):
    g,sg=letter(s); e=secedge(s); return cstart(e) if sg else cend(e)
def nxt(s): return (s[0],(s[1]+1)%L[s[0]])
SE=[('ce',e) for e in CE]+[('sp',s) for s in SEC]
def evs(E):
    if E[0]=='ce': return {('r',cstart(E[1])),('r',cend(E[1]))}
    s=E[1]; return {('c',s[0]),('r',secstart(s))}
def sedges(s): return {('ce',secedge(s)),('sp',s),('sp',nxt(s))}
def svs(s): return set().union(*[evs(E) for E in sedges(s)])
print('counts: vertices',len(SV),'edges',len(SE),'sectors',len(SEC))
# each sector a genuine triangle with 3 distinct vertices, 3 distinct edges, closing up
ok=all(len(svs(s))==3 and len(sedges(s))==3 for s in SEC)
ok2=all(evs(('ce',secedge(s)))=={('r',secstart(s)),('r',secstart(nxt(s)))} for s in SEC)
print('every sector is a triangle (3 vertices,3 edges):',ok,' circle edge joins consecutive spoke feet:',ok2)
# cellular chain complex of the CW complex (orientations: edges start->end; spokes cone->foot; sector boundary cycle)
idxV={v:i for i,v in enumerate(SV)}; idxE={E:i for i,E in enumerate(SE)}
def tail_head(E):
    if E[0]=='ce': return ('r',cstart(E[1])),('r',cend(E[1]))
    s=E[1]; return ('c',s[0]),('r',secstart(s))
d1=sympy.zeros(len(SV),len(SE))
for E in SE:
    t,h=tail_head(E); d1[idxV[h],idxE[E]]+=1; d1[idxV[t],idxE[E]]-=1
d2=sympy.zeros(len(SE),len(SEC))
for j,s in enumerate(SEC):
    # boundary: cone -> foot(s) -> foot(next) -> cone  : +spoke(s) + circle edge (signed by reading direction) - spoke(next)
    g,sg=letter(s); e=('ce',secedge(s))
    d2[idxE[('sp',s)],j]+=1; d2[idxE[e],j]+= (1 if sg else -1); d2[idxE[('sp',nxt(s))],j]-=1
print('d1*d2 == 0:',(d1*d2).is_zero_matrix)
from sympy.matrices.normalforms import smith_normal_form
def homology_Z(dk,dk1,nk):
    # H_k = ker dk / im dk1 ; return (free rank, torsion invariants)
    rk=dk.rank() if dk.shape[1] else 0
    snf=smith_normal_form(dk1,domain=sympy.ZZ) if dk1.shape[1] else None
    inv=[abs(snf[i,i]) for i in range(min(snf.shape)) if snf[i,i]!=0] if snf is not None else []
    free=nk-rk-len(inv)
    return free,[a for a in inv if a!=1]
zero=sympy.zeros(1,len(SV))[:, :0]
print('CW homology H0,H1,H2 over Z (free rank, torsion):',
      homology_Z(sympy.zeros(1,len(SV)),d1,len(SV)), homology_Z(d1,d2,len(SE)), homology_Z(d2,sympy.zeros(len(SEC),0),len(SEC)))
# the flag complex of the Lean comparability graph = order complex; check flagness (cliques=chains) and acyclicity over Q,F2,F3,F5
cells=[('v',v) for v in SV]+[('e',E) for E in SE]+[('s',s) for s in SEC]
def less(a,b):
    if a[0]=='v' and b[0]=='e': return a[1] in evs(b[1])
    if a[0]=='v' and b[0]=='s': return a[1] in svs(b[1])
    if a[0]=='e' and b[0]=='s': return a[1] in sedges(b[1])
    return False
adj={(i,j) for i,a in enumerate(cells) for j,b in enumerate(cells) if less(a,b) or less(b,a)}
n=len(cells); nb={i:{j for j in range(n) if (i,j) in adj} for i in range(n)}
edges=sorted({tuple(sorted(p)) for p in adj})
tris=sorted({tuple(sorted((i,j,k))) for (i,j) in edges for k in nb[i]&nb[j]})
quads=[q for q in itertools.combinations(range(n),4) if False]
k4=sum(1 for (i,j,k) in tris for l in (nb[i]&nb[j]&nb[k]) if l>k)
print('flag complex: vertices',n,'edges',len(edges),'triangles',len(tris),'4-cliques',k4,' Euler char',n-len(edges)+len(tris))
eidx={e:i for i,e in enumerate(edges)}
D1=sympy.zeros(n,len(edges)); D2=sympy.zeros(len(edges),len(tris))
for c,(i,j) in enumerate(edges): D1[j,c]+=1; D1[i,c]-=1
for c,(i,j,k) in enumerate(tris):
    D2[eidx[(j,k)],c]+=1; D2[eidx[(i,k)],c]-=1; D2[eidx[(i,j)],c]+=1
print('order complex homology over Z: H0',homology_Z(sympy.zeros(1,n),D1,n),' H1',homology_Z(D1,D2,len(edges)),' H2',homology_Z(D2,sympy.zeros(len(tris),0),len(tris)))
