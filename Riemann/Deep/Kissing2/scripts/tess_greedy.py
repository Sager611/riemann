"""tess_greedy.py MODE SEED SECONDS -- randomized greedy packing in the 28D PackingStar family (tails = 2*D4 roots).
MODE sq: squares (4 deleted, 8 caps, +4); MODE tess: tesseracts (16 deleted, 64 caps, +48, one tail frame per
tesseract). Exact compatibility: same tail line 2b.2b' <= 32, 60-degree lines <= 48 (orthogonal lines need a shared
endpoint to conflict, excluded by vertex-disjointness). Measured (600 s): sq 202 squares (+808); tess 9-13.
"""
import numpy as np, sys, time, itertools, random
import os
from decode_ps import leech_min
mode=sys.argv[1]; seed=int(sys.argv[2]); TLIM=float(sys.argv[3])
rng=np.random.default_rng(seed)
HERE=os.path.dirname(os.path.abspath(__file__))
rows=[l.strip() for l in open(os.path.join(HERE,'..','..','Kissing','golay_basis.txt')) if l.strip()]
words=[0]
for r in rows:
    b=sum(1<<i for i,c in enumerate(r) if c=='1'); words+=[x^b for x in words]
U=leech_min(words).astype(np.int32)
N=len(U); key={tuple(u):i for i,u in enumerate(U.tolist())}
# 28D tails: 2*D4 roots, 12 lines, 3 frames
tails=[]
for i in range(4):
    for j in range(i+1,4):
        for si in (2,-2):
            for sj in (2,-2):
                t=np.zeros(4,dtype=np.int32); t[i]=si; t[j]=sj; tails.append(t)
def cl(v):
    v=tuple(int(x) for x in v); nz=next(x for x in v if x); return v if nz>0 else tuple(-x for x in v)
lines=sorted(set(cl(t) for t in tails)); L=np.array(lines); LL=L@L.T
nl=len(L)
orth=[[b for b in range(nl) if LL[a,b]==0 and a!=b] for a in range(nl)]
sixty=[[b for b in range(nl) if abs(LL[a,b])==4] for a in range(nl)]
frames=[]
for a in range(nl):
    fr=tuple(sorted([a]+orth[a]))
    if fr not in frames: frames.append(fr)
CAP=20000
M=[np.zeros((CAP,24),dtype=np.int32) for _ in range(nl)]; cnt=[0]*nl
D=np.zeros(N,dtype=bool)
def compatible(newm):
    # newm: dict line -> array of 2b
    for l,arr in newm.items():
        if cnt[l]:
            if (arr@M[l][:cnt[l]].T).max()>32: return False
        for l2 in sixty[l]:
            if cnt[l2] and (arr@M[l2][:cnt[l2]].T).max()>48: return False
        # also within new unit: same line pairs must be <=32 (excluding self), 60 deg lines <=48
    ls=list(newm)
    for i,l in enumerate(ls):
        G=newm[l]@newm[l].T; np.fill_diagonal(G,-999)
        if G.max()>32: return False
        for l2 in ls[i+1:]:
            if abs(LL[l,l2])==4 and (newm[l]@newm[l2].T).max()>48: return False
    return True
def commit(newm,verts):
    for l,arr in newm.items():
        M[l][cnt[l]:cnt[l]+len(arr)]=arr; cnt[l]+=len(arr)
    D[verts]=True
units=0; caps=0; t0=time.time(); tries=0
while time.time()-t0<TLIM:
    a=int(rng.integers(N))
    if D[a]: continue
    av=U[a]; G=U@av; nb=np.flatnonzero((G==16)&(~D))
    if len(nb)<3: continue
    NB=U[nb]
    for _ in range(10):
        tries+=1
        a1=NB[int(rng.integers(len(nb)))]; g1=NB@a1; c2=np.flatnonzero(g1==0)
        if len(c2)==0: continue
        a2=NB[c2[int(rng.integers(len(c2)))]]
        if mode=='sq':
            a4=a1+a2-av; i4=key.get(tuple(a4.tolist()))
            if i4 is None or D[i4]: continue
            u=av-a1; v=av-a2
            vs=[a,key[tuple(a1.tolist())],key[tuple(a2.tolist())],i4]
            if D[vs].any(): continue
            # midpoints: direction u edges: {a,a1},{a2,a4}: 2b = a+a1, a2+a4 ; direction v: {a,a2},{a1,a4}
            mu=np.array([av+a1,a2+a4]); mv=np.array([av+a2,a1+a4])
            opts=[(x,y) for x in range(nl) for y in orth[x]]
            rng.shuffle(opts)
            for (x,y) in opts:
                newm={x:mu,y:mv}
                if compatible(newm):
                    commit(newm,vs); units+=1; caps+=8; break
            else: continue
            break
        else:
            g2=NB@a2; c3=np.flatnonzero((g1==0)&(g2==0))
            if len(c3)==0: continue
            a3=NB[c3[int(rng.integers(len(c3)))]]
            a4=2*av-a1-a2-a3
            us=[av-a1,av-a2,av-a3,av-a4]
            verts=[]; vv={}
            ok=True
            for S in range(16):
                x=av.copy()
                for i in range(4):
                    if S>>i&1: x=x-us[i]
                ix=key.get(tuple(x.tolist()))
                if ix is None or D[ix]: ok=False; break
                verts.append(ix); vv[S]=x
            if not ok: continue
            grp=[np.array([2*vv[S]-us[i] for S in range(16) if not S>>i&1]) for i in range(4)]
            done=False
            fr_order=list(range(len(frames))); rng.shuffle(fr_order)
            for fi in fr_order:
                fr=frames[fi]
                for p in itertools.permutations(fr):
                    newm={p[i]:grp[i] for i in range(4)}
                    if compatible(newm):
                        commit(newm,verts); units+=1; caps+=64; done=True; break
                if done: break
            if done: break
print(mode,'seed',seed,'units',units,'caps',caps,'deleted',int(D.sum()),'gain',caps-int(D.sum()),'per line',cnt,'tries',tries,'time %.0f'%(time.time()-t0))
