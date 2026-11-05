"""p48_fhat_range.py -- for classes of P/2P with minimal norm s = 12..28 the moment equations leave one free count
(a_6); this prints the exact feasible range of fhat(y).  Norm 12: [-4608, 396800]; norm 24: [-51.25M, 38.1M].
"""
from fractions import Fraction as Fr
N=52416000; n=48
def rref(A):
    A=[r[:] for r in A]; m=len(A); k=len(A[0]); r=0; piv=[]
    for c in range(k-1):
        p=next((i for i in range(r,m) if A[i][c]!=0),None)
        if p is None: continue
        A[r],A[p]=A[p],A[r]; A[r]=[x/A[r][c] for x in A[r]]
        for i in range(m):
            if i!=r and A[i][c]!=0:
                f=A[i][c]; A[i]=[x-f*y for x,y in zip(A[i],A[r])]
        piv.append(c); r+=1
    return A[:r],piv
for s in (12,14,16,18,20,22,24,26,28,30,32,36,40,48):
    vals=list(range(0,7))
    A=[]
    for m in range(6):
        row=[Fr((1 if v==0 else 2)*v**(2*m)) for v in vals]
        dfact=1
        for i in range(1,2*m,2): dfact*=i
        den=1
        for i in range(m): den*=(n+2*i)
        A.append(row+[Fr(N)*Fr(6*s)**m*dfact/den])
    R,piv=rref(A)
    free=[c for c in range(7) if c not in piv]
    assert len(free)==1; f=free[0]
    # a_piv = rhs - coef*t, a_f = t
    def sol(t):
        a=[None]*7; a[f]=t
        for i,c in enumerate(piv): a[c]=R[i][-1]-R[i][f]*t
        return a
    # feasible t: all a>=0
    lo=Fr(0); hi=None
    for i,c in enumerate(piv):
        coef=R[i][f]; rhs=R[i][-1]
        # rhs - coef t >= 0
        if coef>0:
            ub=rhs/coef; hi=ub if hi is None or ub<hi else hi
        elif coef<0:
            lb=rhs/coef; lo=lb if lb>lo else lo
        else:
            assert rhs>=0
    fh=lambda a: sum((1 if v==0 else 2)*x*(-1)**v for v,x in zip(vals,a))
    if hi is None or lo>hi: print('norm',s,'infeasible'); continue
    print('norm',s,'free a_%d in [%s, %s]'%(f,float(lo),float(hi)),'fhat range',float(min(fh(sol(lo)),fh(sol(hi)))),float(max(fh(sol(lo)),fh(sol(hi)))))
