
import re, itertools, collections, sys
B='/Users/adriansager/playground/riemann/.scratch/oaicx/pp/A-counterexample-to-Sidorenkos-conjecture-September-23-2026/build/sections/complex.tex'
t=open(B).read()
# Table 1 rows: j & i1&i2&i3 & j12&j13&j23 & bits
rows=re.findall(r'^(\d+)\s*&(\d+)&(\d+)&(\d+)\s*&(\d+)&(\d+)&(\d+)\s*&([01]{3})\\\\',t,re.M)
assert len(rows)==22, len(rows)
faces=[tuple(map(int,r[1:4])) for r in rows]; opp=[tuple(map(int,r[4:7])) for r in rows]; bits=[r[7] for r in rows]
# Lean faces
L=open('/Users/adriansager/playground/riemann/.scratch/oaicx/files/lean/ComparatorChallenges/SidorenkoCounterexample.lean').read()
lf=[tuple(map(int,x.split(','))) for x in re.findall(r'\{(\d+,\d+,\d+)\}',L)]
print('lean faces == table faces:', lf==faces)
# (1) pair in exactly two faces, |E|=33, neighbour columns correct, connected
pc=collections.defaultdict(list)
for j,f in enumerate(faces):
    assert len(set(f))==3
    for a,b in itertools.combinations(f,2): pc[frozenset((a,b))].append(j)
print('|E|=',len(pc),' all pairs in exactly 2 faces:',all(len(v)==2 for v in pc.values()))
ok=True
for j,f in enumerate(faces):
    for (a,b),o in zip([(0,1),(0,2),(1,2)],opp[j]):
        e=frozenset((f[a],f[b])); other=[x for x in pc[e] if x!=j]
        if other!=[o]: ok=False; print('neighbour mismatch face',j,'pair',sorted(e),'table',o,'actual',other)
print('neighbour columns correct:',ok)
adj=collections.defaultdict(set)
for e,v in pc.items(): adj[v[0]].add(v[1]); adj[v[1]].add(v[0])
seen={0}; st=[0]
while st:
    x=st.pop()
    for y in adj[x]:
        if y not in seen: seen.add(y); st.append(y)
print('face-neighbour graph connected:',len(seen)==22)
V=13;E=len(pc);F=22
print('Euler char V-E+F =',V-E+F)
# vertex links: for each point, faces containing it form a cycle?
islink=True
for i in range(13):
    fs=[f for f in faces if i in f]
    g=collections.defaultdict(list)
    for f in fs:
        a,b=[x for x in f if x!=i]; g[a].append(b); g[b].append(a)
    if not all(len(v)==2 for v in g.values()): islink=False; continue
    # connected single cycle
    s=next(iter(g)); seen2={s}; st=[s]
    while st:
        x=st.pop()
        for y in g[x]:
            if y not in seen2: seen2.add(y); st.append(y)
    if len(seen2)!=len(g): islink=False
print('all vertex links are single cycles (closed surface):',islink)
# (2) bit classes
for pos in range(3):
    for val in '01':
        cls=[j for j in range(22) if bits[j][pos]==val]
        pts=set().union(*[faces[j] for j in cls])
        print('bit',pos+1,'val',val,'size',len(cls),'covers all 13:',len(pts)==13)
# Table 2 orders
t2=re.findall(r'^(\d)&(\d)&\$([\d,]+)\$\s*\n\s*&\$([\d,]+)\$',t,re.M)
print('table2 rows',len(t2))
for pos,val,order,newp in t2:
    order=list(map(int,order.split(','))); newp=list(map(int,newp.split(',')))
    cls=sorted(j for j in range(22) if bits[j][int(pos)-1]==val)
    good = sorted(order)==cls and len(newp)==10
    covered=set(faces[order[0]])
    for j,np_ in zip(order[1:],newp):
        f=set(faces[j])
        shares=any(len(f & set(faces[k]))==2 for k in order[:order.index(j)])
        new=f-covered
        if not(shares and new=={np_}): good=False; print('  order fail pos',pos,val,'face',j,'new',new,'claimed',np_,'shares',shares)
        covered|=f
    print('order pos',pos,'val',val,'valid:',good and len(covered)==13)
# (3) a_e
ae={e:sum(bits[v[0]][k]!=bits[v[1]][k] for k in range(3))/3 for e,v in pc.items()}
print('min a_e',min(ae.values()),'max a_e',max(ae.values()),'a_{0,2}=',ae[frozenset((0,2))])
# point degree table
deg={i:sum(1 for e in pc if i in e) for i in range(13)}
nb={i:[next(iter(e-{i})) for e in pc if i in e] for i in range(13)}
tab={}
for m in re.finditer(r'(\d+)&(\d)&(\d)&(\d)&(\d)',t[t.find('eq:complex-colors')-2000:t.find('eq:complex-colors')]):
    tab[int(m.group(1))]=tuple(map(int,m.groups()[1:]))
okd=True
for i in range(13):
    act=(deg[i],sum(deg[k]==4 for k in nb[i]),sum(deg[k]==5 for k in nb[i]),sum(deg[k]==6 for k in nb[i]))
    if tab.get(i)!=act: okd=False; print('deg table mismatch',i,tab.get(i),act)
print('degree table rows parsed',len(tab),'correct:',okd)
# (4) colour refinement on H with separate colours for the two parts
Vh=[('p',i) for i in range(13)]+[('f',j) for j in range(22)]
nbH={('p',i):[('f',j) for j in range(22) if i in faces[j]] for i in range(13)}
nbH.update({('f',j):[('p',i) for i in faces[j]] for j in range(22)})
col={v:v[0] for v in Vh}
for it in range(40):
    sig={v:(col[v],tuple(sorted(collections.Counter(col[u] for u in nbH[v]).items()))) for v in Vh}
    m={s:k for k,s in enumerate(sorted(set(sig.values()),key=repr))}
    new={v:m[sig[v]] for v in Vh}
    if len(set(new.values()))==len(set(col.values())): break
    col=new
print('colour refinement classes:',len(set(col.values())),'(35 = discrete)')
# extra: H bipartite, degrees
print('H: vertices',35,'edges',sum(len(f) for f in faces),'point degrees in H',sorted(sum(i in f for f in faces) for i in range(13)))
