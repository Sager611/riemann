
# Independent check of the explicit 10-dim algebra C and T = C x| DC used in BOTH the Auslander-Reiten and the
# Tachikawa counterexample preprints, and of the claimed  Ext_T^*(s,s) = k[tau], |tau| = 3.
# Parameter q is specialised to a primitive element of GF(2^20) (multiplicative order 2^20-1 = 1048575),
# a generic point for these low-degree computations (ranks over F_2(q) equal generic ranks).
import numpy as np, re, itertools, sys
NBIT=20; POLY=(1<<20)|(1<<3)|1          # x^20 + x^3 + 1
N=(1<<NBIT)-1
exp=np.zeros(2*N+2,dtype=np.int64); log=np.full(1<<NBIT,-1,dtype=np.int64)
a=1
for i in range(N):
    exp[i]=a; log[a]=i
    a<<=1
    if a>>NBIT: a^=POLY
assert a==1 and (log[1:]>=0).all(), 'modulus not primitive'
exp[N:2*N]=exp[:N]
def mul(a,b):
    a=np.asarray(a,dtype=np.int64); b=np.asarray(b,dtype=np.int64)
    r=exp[(log[a]+log[b]) % N]
    return np.where((a==0)|(b==0),0,r)
def inv(a): return int(exp[(N-log[a])%N])
def smul(c,a): return int(mul(c,a))
q=2  # the class of x: primitive
one=1
# ---------------- C: basis and products (AR paper rules) ----------------
CB=['e','x','y','z','u','v','t','j','f','n']; ci={b:i for i,b in enumerate(CB)}
corner={'e':('e','e'),'x':('e','e'),'y':('e','e'),'z':('e','e'),'u':('e','f'),'v':('e','f'),'t':('f','e'),'j':('f','e'),'f':('f','f'),'n':('f','f')}
q2=smul(q,q); q1=q^1
rules={('x','y'):{'z':q},('y','x'):{'z':1},('x','u'):{'v':1},('y','u'):{'v':1},('t','x'):{'j':1},('t','y'):{'j':q2},
       ('u','t'):{'y':1,'x':q},('v','t'):{'z':q},('u','j'):{'z':1},('t','u'):{'n':q1},('n','t'):{'j':q},('u','n'):{'v':1}}
def cprod(a,b):
    out=np.zeros(10,dtype=np.int64)
    if a in 'ef':
        if corner[b][0]==a: out[ci[b]]=1
        return out
    if b in 'ef':
        if corner[a][1]==b: out[ci[a]]=1
        return out
    for k,v in rules.get((a,b),{}).items(): out[ci[k]]^=v
    return out
CM=np.array([[cprod(a,b) for b in CB] for a in CB])   # CM[a,b,:]
# cross-check with the full multiplication table printed in the Tachikawa preprint
tex=open('/Users/adriansager/playground/riemann/.scratch/oaicx/pp/A-counterexample-to-Tachikawas-second-conjecture-September-23-2026/build/sections/algebra.tex').read()
st=tex.index(r'\begin{array}{c|cccccccccc}'); blk=tex[st:tex.index(r'\end{array}',st)]
rows=[r for r in blk.split('\\\\') if '&' in r][1:]
cols=['e','f','x','y','z','u','v','t','j','n']
def parse_entry(s):
    s=s.strip(); out=np.zeros(10,dtype=np.int64)
    if s=='0': return out
    terms=re.findall(r'(\(1\+q\)|q\^2|q)?([a-z])',s.replace(' ',''))
    assert '+'.join(a+b for a,b in terms)==s.replace(' ',''), s
    for c,v in terms:
        coef={'':1,'q':q,'q^2':q2,'(1+q)':q1}[c]
        out[ci[v]]^=coef
    return out
mism=0
for r in rows:
    ents=[x for x in r.replace(r'\hline','').split('&')]
    a=ents[0].strip()
    for b,s in zip(cols,ents[1:]):
        if not np.array_equal(parse_entry(s),CM[ci[a],ci[b]]): mism+=1; print('table mismatch',a,b,s)
print('Tachikawa full table rows parsed:',len(rows),' mismatches vs AR rules:',mism)
def cmul_vec(u,v):
    out=np.zeros(10,dtype=np.int64)
    for i in np.nonzero(u)[0]:
        for j in np.nonzero(v)[0]:
            out^=mul(mul(u[i],v[j]),CM[i,j])
    return out
E=np.eye(10,dtype=np.int64)
assoc=all(np.array_equal(cmul_vec(cmul_vec(E[i],E[j]),E[k]),cmul_vec(E[i],cmul_vec(E[j],E[k]))) for i in range(10) for j in range(10) for k in range(10))
unit=E[ci['e']]^E[ci['f']]
unital=all(np.array_equal(cmul_vec(unit,E[i]),E[i]) and np.array_equal(cmul_vec(E[i],unit),E[i]) for i in range(10))
utut=cmul_vec(cmul_vec(cmul_vec(E[ci['u']],E[ci['t']]),E[ci['u']]),E[ci['t']])
print('C associative:',assoc,' unital:',unital,' utut == q(1+q) z:',np.array_equal(utut,mul(smul(q,q1),E[ci['z']])))
# ---------------- T = C x| DC ----------------
TB=CB+[b.upper() for b in CB]; ti={b:i for i,b in enumerate(TB)}
# product structure constants TM[a,b,:] over the 20-dim basis
TM=np.zeros((20,20,20),dtype=np.int64)
for i in range(10):
    for j in range(10):
        TM[i,j,:10]=CM[i,j]
        # a * b^*  : (a psi)(c) = psi(c a), psi=b^*  -> coefficient of c^* is b^*(c a) = CM[c,a][b]
        for c in range(10): TM[i,10+j,10+c]=CM[c,i,j]
        # a^* * b  : (phi b)(c) = phi(b c), phi=a^* -> coefficient of c^* is CM[b,c][a]
        for c in range(10): TM[10+i,j,10+c]=CM[j,c,i]
def tmul(u,v):
    out=np.zeros(20,dtype=np.int64)
    for i in np.nonzero(u)[0]:
        for j in np.nonzero(v)[0]:
            out^=mul(mul(u[i],v[j]),TM[i,j])
    return out
ET=np.eye(20,dtype=np.int64)
assocT=all(np.array_equal(tmul(tmul(ET[i],ET[j]),ET[k]),tmul(ET[i],tmul(ET[j],ET[k]))) for i in range(20) for j in range(20) for k in range(20))
lam=lambda w: w[10+ci['e']]^w[10+ci['f']]            # lambda(a,phi)=phi(1)=phi(e)+phi(f)
G=np.array([[lam(TM[i,j]) for j in range(20)] for i in range(20)])
def rank(Mat):
    A=np.array(Mat,dtype=np.int64).copy(); r=0; rows,cols=A.shape
    for c in range(cols):
        p=next((i for i in range(r,rows) if A[i,c]),None)
        if p is None: continue
        A[[r,p]]=A[[p,r]]; A[r]=mul(inv(A[r,c]),A[r])
        f=A[:,c].copy(); f[r]=0
        nz=np.nonzero(f)[0]
        if len(nz): A[nz]^=mul(f[nz][:,None],A[r][None,:])
        r+=1
        if r==rows: break
    return r
print('T associative:',assocT,' trace form symmetric:',np.array_equal(G,G.T),' nondegenerate (rank 20):',rank(G)==20)
# ---------------- minimal projective resolution of the simple s at f ----------------
def rref(A):
    A=np.array(A,dtype=np.int64).copy(); r=0; piv=[]
    rows,cols=A.shape
    for c in range(cols):
        if r==rows: break
        p=next((i for i in range(r,rows) if A[i,c]),None)
        if p is None: continue
        A[[r,p]]=A[[p,r]]; A[r]=mul(inv(A[r,c]),A[r])
        f=A[:,c].copy(); f[r]=0
        nz=np.nonzero(f)[0]
        if len(nz): A[nz]^=mul(f[nz][:,None],A[r][None,:])
        piv.append(c); r+=1
    return A[:r],piv
def nullspace_rows(Phi):
    # vectors c with c @ Phi = 0  ->  nullspace of Phi^T
    A,piv=rref(Phi.T); n=Phi.shape[0]; free=[c for c in range(n) if c not in piv]; out=[]
    for fc in free:
        v=np.zeros(n,dtype=np.int64); v[fc]=1
        for r,pc in enumerate(piv): v[pc]=A[r,fc]   # char 2: -x = x
        out.append(v)
    return np.array(out,dtype=np.int64).reshape(len(out),n)
idem={'e':ti['e'],'f':ti['f']}
def right_idem(i):
    return 'e' if tmul(ET[i],ET[idem['e']]).any() else 'f'
PB={k:[i for i in range(20) if right_idem(i)==k and tmul(ET[i],ET[idem[k]])[i]] for k in 'ef'}
print('dim Te =',len(PB['e']),' dim Tf =',len(PB['f']))
# left multiplication matrices on Te, Tf (coordinates in PB lists)
def Lmat(w,k):
    B=PB[k]; M=np.zeros((len(B),len(B)),dtype=np.int64)
    for c,i in enumerate(B):
        img=tmul(ET[w],ET[i])
        assert all(img[j]==0 for j in range(20) if j not in B)
        for r,j in enumerate(B): M[r,c]=img[j]
    return M
L={(w,k):Lmat(w,k) for w in range(20) for k in 'ef'}
radT=[i for i in range(20) if TB[i] not in ('e','f')]
def act(w,blocks,vecs):
    # left action of basis element w on row vectors (ambient = concatenation of blocks)
    out=np.zeros_like(vecs); o=0
    for k in blocks:
        d=len(PB[k]); Mw=L[(w,k)]
        seg=vecs[:,o:o+d]
        # out_seg = (Mw @ seg^T)^T  over GF(2^20)
        res=np.zeros_like(seg)
        for c in range(d):
            col=Mw[:,c]
            nz=np.nonzero(col)[0]
            for r in nz: res[:,r]^=mul(col[r],seg[:,c])
        out[:,o:o+d]=res; o+=d
    return out
def span(vecs):
    if len(vecs)==0: return vecs
    A,_=rref(vecs); return A
# Omega^1(s) = rad(Tf) inside ambient ['f']
blocks=['f']; d=len(PB['f'])
M=np.array([np.eye(d,dtype=np.int64)[r] for r,i in enumerate(PB['f']) if TB[i]!='f'],dtype=np.int64)
ext_ss=[1]; ext_se=[0]; dims=[1]
for step in range(1,13):
    M=span(M); dims.append(len(M))
    R=span(np.vstack([act(w,blocks,M) for w in radT])) if len(M) else M
    eM=span(act(idem['e'],blocks,M)); fM=span(act(idem['f'],blocks,M))
    eR=span(act(idem['e'],blocks,R)) if len(R) else R; fR=span(act(idem['f'],blocks,R)) if len(R) else R
    gens=[]
    for kk,(XM,XR) in (('e',(eM,eR)),('f',(fM,fR))):
        base=XR.copy() if len(XR) else np.zeros((0,M.shape[1]),dtype=np.int64)
        cur=rank(base) if len(base) else 0
        for v in XM:
            t=np.vstack([base,v[None,:]]); rk=rank(t)
            if rk>cur: base=t; cur=rk; gens.append((kk,v))
    a=sum(1 for g in gens if g[0]=='e'); b=sum(1 for g in gens if g[0]=='f')
    ext_se.append(a); ext_ss.append(b)
    # projective cover phi: domain blocks = types of gens
    rowsPhi=[]
    for kk,v in gens:
        for i in PB[kk]: rowsPhi.append(act(i,blocks,v[None,:])[0])
    Phi=np.array(rowsPhi,dtype=np.int64)
    assert rank(Phi)==len(M), 'cover not surjective'
    K=nullspace_rows(Phi)
    blocks=[g[0] for g in gens]; M=K
    print(f'i={step}: dim Omega^{step}={dims[-1]}  top: e-part {a}, f-part {b}  -> dim Ext^{step}(s,s)={b}, Ext^{step}(s,s_e)={a}',flush=True)
print('Ext^i(s,s), i=0..12 :',ext_ss)
print('claimed k[tau],|tau|=3:',[1 if i%3==0 else 0 for i in range(13)])
print('MATCH:',ext_ss==[1 if i%3==0 else 0 for i in range(13)])
