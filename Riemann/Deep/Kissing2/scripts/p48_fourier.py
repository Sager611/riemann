"""p48_fourier.py -- exact inner-product distribution of the 52 416 000 minimal vectors of an extremal even unimodular
48-dim lattice (11-design) against a vector y of norm 6, 8, 10 (minimal in its class mod 2P), and
fhat(y) = sum_u (-1)^(u.y): +1833472, +380416, +64000.  Pure fractions, < 1 s.
"""
from fractions import Fraction as Fr
import itertools
N=52416000; n=48
def moments(normy, vals):
    # unknown counts a_v for v in vals (v>=0), symmetric; known: none
    # equations: sum_v c_v a_v v^{2m} = N * (6*normy)^m * (2m-1)!! / prod_{i<m}(n+2i)  for m=0..5
    rows=[];rhs=[]
    for m in range(6):
        row=[]
        for v in vals:
            mult = 1 if v==0 else 2
            row.append(Fr(mult*(v**(2*m))))
        dfact=1
        for i in range(1,2*m,2): dfact*=i
        den=1
        for i in range(m): den*=(n+2*i)
        rows.append(row); rhs.append(Fr(N)*Fr(6*normy)**m*dfact/den)
    return rows,rhs
def solve(rows,rhs):
    import copy
    A=[r[:]+[b] for r,b in zip(rows,rhs)]
    m=len(A); k=len(A[0])-1
    piv=[]; r=0
    for c in range(k):
        p=next((i for i in range(r,m) if A[i][c]!=0),None)
        if p is None: continue
        A[r],A[p]=A[p],A[r]
        A[r]=[x/A[r][c] for x in A[r]]
        for i in range(m):
            if i!=r and A[i][c]!=0:
                f=A[i][c]; A[i]=[x-f*y for x,y in zip(A[i],A[r])]
        piv.append(c); r+=1
    # check consistency
    for i in range(r,m):
        assert A[i][-1]==0, ('inconsistent',A[i][-1])
    sol=[Fr(0)]*k
    for i,c in enumerate(piv): sol[c]=A[i][-1]
    return sol, len(piv)
# y min (norm 6): values 0,1,2,3 plus known a_6=1 -> move to rhs
for normy, vals, known in [(6,[0,1,2,3],{6:1}),(8,[0,1,2,3,4],{}),(10,[0,1,2,3,4,5],{})]:
    rows,rhs=moments(normy,vals)
    for v,c in known.items():
        for m in range(6): rhs[m]-=2*c*Fr(v)**(2*m)
    sol,rk=solve(rows,rhs)
    d={v:s for v,s in zip(vals,sol)}; d.update(known)
    tot=sum((1 if v==0 else 2)*a for v,a in d.items())
    fh=sum((1 if v==0 else 2)*a*(-1)**v for v,a in d.items())
    print('norm',normy,'rank',rk,'dist',{v:str(a) for v,a in d.items()},'total',tot,'fhat',fh)
