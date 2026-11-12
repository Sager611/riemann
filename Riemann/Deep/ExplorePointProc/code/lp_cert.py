
# explore40: LP with safety margins for the augmented criterion at one lambda; writes rational certificate data.
# f(a) = sum_n B_n H_{2n}(sqrt(2pi) a/s) e^{-pi a^2/s^2}  (physicists' Hermite, B_n rational), fhat(u) = s sum_n (-1)^n B_n H_{2n}(sqrt(2pi) s u) e^{-pi s^2 u^2}
import numpy as np, sys, time, threading, os, resource, json
from fractions import Fraction as Fr
def _wd(cap=285e6):
    while True:
        if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss > cap:
            print("MEMCAP exceeded, abort", flush=True); os._exit(3)
        time.sleep(0.1)
threading.Thread(target=_wd, daemon=True).start()
from scipy.optimize import linprog
import scipy.sparse as sp
from math import lgamma, log, pi, sqrt, exp
import importlib.util
spec=importlib.util.spec_from_file_location('pg','primal_2gap.py'); pg=importlib.util.module_from_spec(spec); spec.loader.exec_module(pg)

lam=Fr(sys.argv[1]); theta=Fr(sys.argv[2]); D=int(sys.argv[3]); s=Fr(sys.argv[4]); U=Fr(sys.argv[5]); dg=Fr(sys.argv[6])
tau1=float(sys.argv[7]); tau2=float(sys.argv[8]); tau3=float(sys.argv[9]); out=sys.argv[10]
assert U<=3*lam and D%2==1
lamf=float(lam); thf=float(theta); sf=float(s); Uf=float(U)
nc=D+1
al=np.unique(np.r_[np.arange(0,3.0,0.0005), np.arange(3.0,float(os.environ.get('AMAX','12')),0.004)]); Bf=pg.basis_f(D,sf,al)
x,wq=np.polynomial.legendre.leggauss(240); aq=(x+1)*thf/2; wq=wq*thf/2
Jvec=pg.basis_f(D,sf,[0.0])[0]+2*(pg.basis_f(D,sf,aq)*(wq*aq)[:,None]).sum(0)
fh0=pg.basis_fh(D,sf,[0.0])[0]
uu=np.arange(Uf,float(os.environ.get('UMAX','12')),0.002); Bu=pg.basis_fh(D,sf,uu)
ng=int((U-lam)/dg)+1; g=np.array([float(lam+k*dg) for k in range(ng)]); assert abs(g[-1]-Uf)<1e-12, (g[-1],Uf)
nv=nc+ng+2; iA=nc; ib=nc+ng; icc=nc+ng+1
def pad(M): return sp.hstack([sp.csr_matrix(M), sp.csr_matrix((M.shape[0], nv-nc))]).tocsr()
sf_=np.abs(Bf).max(1); su_=np.abs(Bu).max(1)
ok_f=sf_>1e-250; ok_u=su_>1e-250
Bf=Bf[ok_f]/sf_[ok_f,None]; Bu=Bu[ok_u]/su_[ok_u,None]
A_ub=[pad(-Bf), pad(Bu)]; b_ub=[-tau1*np.ones(Bf.shape[0]), -tau2*np.ones(Bu.shape[0])]
print('sign rows', Bf.shape[0], Bu.shape[0])
Bg=pg.basis_fh(D,sf,g)
dense=np.zeros((ng*ng,nc)); ri=[]; ci=[]; vv=[]
for i1 in range(ng):
    s12=g[i1]+g; B12=pg.basis_fh(D,sf,s12)          # keep fhat(g1+g2) for all (<=0 beyond U anyway)
    dense[i1*ng:(i1+1)*ng]=Bg[i1][None,:]+B12
    for i2 in range(ng):
        rr=i1*ng+i2
        if i2!=i1: ri+=[rr,rr]; ci+=[i1,i2]; vv+=[1.0,-1.0]
        ri+=[rr,rr]; ci+=[ng,ng+1]; vv+=[1-g[i1],-1.0]
A_ub.append(sp.hstack([sp.csr_matrix(dense), sp.csr_matrix((vv,(ri,ci)),shape=(ng*ng,ng+2))]).tocsr()); b_ub.append(-tau3*np.ones(ng*ng))
# 1D: g2>=U (A(g2)=A(U), fhat(g1+g2) dropped):  fhat(g1)+A(g1)-A(U)+beta(1-g1) <= c - tau3
M=np.zeros((ng,nv)); M[:,:nc]=Bg; M[np.arange(ng),iA+np.arange(ng)]+=1; M[:,iA+ng-1]-=1; M[:,ib]=1-g; M[:,icc]=-1
A_ub.append(sp.csr_matrix(M)); b_ub.append(-tau3*np.ones(ng))
# 1D: g1>=U: A(U)-A(g2)+beta(1-U) <= c - tau3
M=np.zeros((ng,nv)); M[:,iA+ng-1]+=1; M[np.arange(ng),iA+np.arange(ng)]-=1; M[:,ib]=1-Uf; M[:,icc]=-1
A_ub.append(sp.csr_matrix(M)); b_ub.append(-tau3*np.ones(ng))
A_ub=sp.vstack(A_ub).tocsr(); b_ub=np.concatenate(b_ub)
A_eq=np.zeros((2,nv)); A_eq[0,:nc]=fh0; A_eq[1,iA+ng-1]=1.0; b_eq=[1.0,0.0]
cobj=np.zeros(nv); cobj[:nc]=-Jvec; cobj[icc]=2.0
bounds=[(None,None)]*nc+[(None,None)]*ng+[(0,None),(None,None)]
r=linprog(cobj,A_ub=A_ub,b_ub=b_ub,A_eq=A_eq,b_eq=b_eq,bounds=bounds,method='highs')
print('LP status',r.status,r.message)
cn=r.x[:nc]; J=Jvec@cn; c=r.x[icc]; beta=r.x[ib]; A=r.x[iA:iA+ng]
print('float: J=%.6f c=%.6f beta=%.5f margin J-1-2c=%.6f'%(J,c,beta,J-1-2*c))
# convert normalized coefficients c_n (of h_{2n}) to B_n (of H_{2n} e^{-y^2/2}): h_m = (2^m m! sqrt(pi))^{-1/2} H_m e^{-y^2/2}
Bn=[cn[n]*exp(-0.5*(2*n*log(2)+lgamma(2*n+1)+0.5*log(pi))) for n in range(nc)]
Bq=[Fr(repr(float('%.17g'%b))) for b in Bn]   # exact rationals from doubles
json.dump({'lam':str(lam),'theta':str(theta),'D':D,'s':str(s),'U':str(U),'dg':str(dg),
           'B':[str(b) for b in Bq],'A':[str(Fr(repr(float(a))).limit_denominator(10**12)) for a in A],
           'beta':str(Fr(repr(float(beta))).limit_denominator(10**12)),'c_float':c,'J_float':J}, open(out,'w'))
print('saved',out)
