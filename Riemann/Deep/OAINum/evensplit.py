
# evensplit.py -- Prop ff:even-split (pair-orientation tournament) over A (x) A, masks, digit logs
import sys, random, itertools
import numpy as np
from ffpoly import *

def two_primary_gen(p):
    s=0;t=p-1
    while t%2==0: t//=2; s+=1
    for g in range(2,p):
        if pow(g,(p-1)//2,p)==p-1: return pow(g,t,p),s,t

class TA:
    # R = F_p[T1,T2]/(F(T1),F(T2)), elements N x N int64 arrays
    def __init__(self,F,p):
        self.F=F; self.p=p; self.N=deg(F)
        self.red=np.array([(-c)%p for c in F[:-1]],dtype=np.int64)  # T^N = sum red[k] T^k
    def reduce_axis(self,A,axis):
        p=self.p; N=self.N
        A=np.moveaxis(A,axis,0).copy()
        for i in range(A.shape[0]-1,N-1,-1):
            row=A[i]
            if row.any():
                A[i-N:i]=(A[i-N:i]+self.red[:,None]*row[None,:])%p
        A=A[:N]
        return np.moveaxis(A,0,axis)
    def mul(self,X,Y):
        p=self.p; N=self.N
        Rm=np.zeros((2*N-1,2*N-1),dtype=np.int64)
        for i in range(N):
            for j in range(N):
                c=int(X[i,j])
                if c: Rm[i:i+N,j:j+N]=(Rm[i:i+N,j:j+N]+c*Y)%p
        Rm=self.reduce_axis(Rm,0); Rm=self.reduce_axis(Rm,1)
        return Rm%p
    def pow(self,X,e):
        Rr=self.one()
        while e:
            if e&1: Rr=self.mul(Rr,X)
            e>>=1
            if e: X=self.mul(X,X)
        return Rr
    def one(self):
        o=np.zeros((self.N,self.N),dtype=np.int64); o[0,0]=1; return o
    def scal(self,c,X): return (c*X)%self.p

def power_sums(F,p):
    # s_j = sum of j-th powers of roots, j=0..N-1, via Newton identities (F monic)
    N=deg(F); e=[1]
    # elementary symmetric: F = X^N - e1 X^{N-1} + e2 ... -> coefficient a_{N-k} = (-1)^k e_k
    el=[1]+[((-1)**k*F[N-k])%p for k in range(1,N+1)]
    s=[N%p]
    for k in range(1,N):
        v=(-1)**(k-1)*k*el[k]
        for i in range(1,k):
            v+=(-1)**(i-1)*el[i]*s[k-i]
        s.append(v%p)
    return s

def even_split(F,p,check_roots=None):
    N=deg(F); R=TA(F,p)
    u,s,t=two_primary_gen(p)
    D=np.zeros((N,N),dtype=np.int64); D[1,0]=1; D[0,1]=(-1)%p
    e=R.pow(D,p-1)
    z0=R.pow(D,t)
    z0=R.mul(e,z0)
    pieces=[(e,0)]
    uinv=pow(u,p-2,p)
    for j in range(s):
        newp=[]
        for (m,ell) in pieces:
            z=R.mul(m,R.scal(pow(uinv,ell,p),z0))
            z=R.pow(z,2**(s-j-1))
            z=R.mul(m,z)
            for c in (0,1):
                diff=(z-(1 if c==0 else p-1)*m)%p
                mask=(m-R.pow(R.mul(m,diff),p-1))%p
                if mask.any(): newp.append((mask,ell+c*2**j))
        pieces=newp
    estar=np.zeros((N,N),dtype=np.int64)
    for (m,ell) in pieces:
        if ell<2**(s-1): estar=(estar+m)%p
    sj=power_sums(F,p)
    cvec=[int(sum(int(estar[i,j])*sj[j] for j in range(N))%p) for i in range(N)]
    c=norm(cvec)
    fac=None
    for a in range(N):
        g=pgcd(F,psub(c,[a],p),p)
        if 0<deg(g)<N: fac=g; break
    orient_ok=None
    if check_roots is not None:
        xs=check_roots
        def ev(X,a,b):
            return int(sum(int(X[i,j])*pow(a,i,p)*pow(b,j,p) for i in range(N) for j in range(N))%p)
        orient_ok=all(ev(estar,xs[i],xs[j])+ev(estar,xs[j],xs[i])==1 for i in range(N) for j in range(N) if i!=j) and all(ev(estar,x,x)==0 for x in xs)
    return fac,orient_ok

if __name__=='__main__':
    random.seed(1)
    tot=0; fails=0; orient_bad=0
    for p in [5,13,17,29,37,41,7,11,19,23,31,43,97,101,257,10009,65537]:
        for N in [2,4,6,8]:
            if N>=p: continue
            if p<=17 and N<=4:
                subsets=list(itertools.combinations(range(p),N))
            else:
                subsets=[random.sample(range(p),N) for _ in range(12)]
            for xs in subsets:
                F=[1]
                for x in xs: F=pmul(F,[(-x)%p,1],p)
                fac,ok=even_split(F,p,list(xs))
                tot+=1
                if fac is None: fails+=1; print('FAIL',p,xs)
                if ok is False: orient_bad+=1; print('ORIENT',p,xs)
    print('even split: cases',tot,'no factor',fails,'orientation violations',orient_bad)

