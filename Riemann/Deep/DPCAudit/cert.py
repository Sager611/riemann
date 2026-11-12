
# Rigorous Arb check of DPC large-gap margin (TeX 950-966) with OUR coefficients b (from top eigenvector).
import sys, time, itertools, math
import numpy as np
from flint import arb, arb_mat, ctx
t0=time.time()
ctx.prec=int(sys.argv[1]) if len(sys.argv)>1 else 1000
RMAX,DMAX,KT=12,5,40
ell=arb(39327)/20000; th=ell*ell
phi=arb(363846723)/100000000; c=arb.pi()*phi/2
E=[arb(5)/4, arb(3)/2, arb(2), arb(3)]
def edot(n): return sum((E[i]*n[i] for i in range(4)), arb(0))
def vecs(maxdeg):
    out=[]
    for d in range(maxdeg+1):
        for n in itertools.product(range(d+1),repeat=4):
            if sum(n)==d: out.append(n)
    return out
fact=[math.factorial(i) for i in range(120)]
def mfact(n):
    p=1
    for x in n: p*=fact[x]
    return p
allv=vecs(2*DMAX)
G={n: th*edot(n).gamma()/mfact(n) for n in allv if sum(n)>0}
Fh=[{(0,0,0,0):arb(1)}]; Gh=[{n:G[n] for n in G if sum(n)==d} for d in range(2*DMAX+1)]
for d in range(1,2*DMAX+1):
    acc={}
    for k in range(1,d+1):
        for m,gv in Gh[k].items():
            for n,fv in Fh[d-k].items():
                key=tuple(m[i]+n[i] for i in range(4)); acc[key]=acc.get(key,arb(0))+k*gv*fv
    Fh.append({key:val/d for key,val in acc.items()})
M={}
for d in range(2*DMAX+1):
    for n,fv in Fh[d].items(): M[n]=th.gamma()/(th+edot(n)).gamma()*fv*mfact(n)
idt=vecs(DMAX); NT=len(idt); NV=RMAX+1
At=arb_mat(NT,NT)
for a,j in enumerate(idt):
    for b_,k in enumerate(idt): At[a,b_]=M[tuple(j[i]+k[i] for i in range(4))]
Av=arb_mat(NV,NV)
for r in range(NV):
    for s in range(NV): Av[r,s]=1/(th+r+s)
def chol(A):
    n=A.nrows(); L=arb_mat(n,n)
    for i in range(n):
        s=A[i,i]
        for k in range(i): s-=L[i,k]*L[i,k]
        L[i,i]=s.sqrt()
        for j in range(i+1,n):
            s=A[j,i]
            for k in range(i): s-=L[j,k]*L[i,k]
            L[j,i]=s/L[i,i]
    return L
Lv=chol(Av); Lt=chol(At); LviT=Lv.inv().transpose(); LtiT=Lt.inv().transpose()
Js=np.load('/Users/adriansager/playground/riemann/.scratch/dpcaudit/Jsym_363846723_100000000.npy')
w,V=np.linalg.eigh(Js); y=V[:,-1]; print("float lam_max=%.15f"%w[-1], "t=%.1f"%(time.time()-t0))
# b = (Lv^{-T} kron Lt^{-T}) y ; take exact midpoints
Y=arb_mat(NV,NT)
for r in range(NV):
    for a in range(NT): Y[r,a]=arb(float(y[r*NT+a]))
Bm=LviT*Y*LtiT.transpose()   # b[r,:] = sum_r' LviT[r,r'] (LtiT y_r')  -> matrix form B = LviT * Y * LtiT^T
B=arb_mat(NV,NT)
for r in range(NV):
    for a in range(NT): B[r,a]=arb(Bm[r,a].mid())   # exact dyadic coefficients b_{r,j}
sumabs=sum((abs(B[r,a]) for r in range(NV) for a in range(NT)), arb(0))
# I(b) = sum_{r,s} Av[r,s] b_r^T At b_s
U=B*At            # U[r,l] = sum_j b_{r,j} M(j+l)
Ival=arb(0)
for r in range(NV):
    for s in range(NV):
        Ival+=Av[r,s]*sum((U[r,a]*B[s,a] for a in range(NT)), arb(0))
# J_K(b) streaming over sparse pairs
cpow=[c**(2*kk+1)/fact[2*kk+1] for kk in range(KT+1)]
Jval=arb(0); np_=0
for b_,k in enumerate(idt):
    be=edot(k)
    for a,l in enumerate(idt):
        if not all(l[i]<=k[i] for i in range(4)): continue
        np_+=1
        bc=1
        for i in range(4): bc*=math.comb(k[i],l[i])
        al=edot(l)
        g2=[(be-al+2*kk+1).gamma() for kk in range(KT+1)]
        for r in range(NV):
            g1=(th+r+al).gamma()
            coef=[(-1)**kk*cpow[kk]*g1*g2[kk]/(th+r+be+2*kk+1).gamma() for kk in range(KT+1)]
            ur=U[r,a]
            for s in range(NV):
                Wv=sum((coef[kk]/(th+r+s+2*kk+1) for kk in range(KT+1)), arb(0))
                Jval+=bc*ur*Wv*B[s,b_]
rem=sumabs**2/(th*(th+1))*c**(2*KT+3)/fact[2*KT+3]
lam=Jval/Ival; fl=4*ell/arb.pi()
marg=1-phi+fl*lam; marg_rig=1-phi+fl*(Jval-rem)/Ival
print("pairs",np_,"sum|b| =",sumabs.str(8)," I =",Ival.str(20)," J_K =",Jval.str(20))
print("J_K/I =",lam.str(20))
print("Taylor remainder bound rem =",rem.str(5)," rem/I =",(rem/Ival).str(5))
print("margin (J_K)        =",marg.str(15))
print("rigorous margin >=  =",marg_rig.str(15), " positive:", bool(marg_rig>0))
print("t=%.1f"%(time.time()-t0))
