
import numpy as np, math, sys, time, json
from scipy.optimize import linprog
from scipy import sparse
from numpy.polynomial import chebyshev as C
x=6e12; L=math.log(x/(4*math.pi)); pr=L/4; us=L/2
th=0.8733458561067185768; y0=2*th-1
ROW=2*pr*math.cosh(pr*y0)   # row fundamental amplitude per unit v(u*)
d=np.load("/tmp/ceil3_lambda.npz"); u=d["u"]; mu=d["mu"]
UEX=math.log(1e7)
def points(U_max, ucut=math.log(3000.), bw=0.002, du=0.01):
    m1=u<ucut
    pu=[u[m1]]; pw=[mu[m1]]
    m2=~m1
    b=np.floor((u[m2]-ucut)/bw).astype(np.int64)
    nb=b.max()+1
    wsum=np.bincount(b,weights=mu[m2],minlength=nb); usum=np.bincount(b,weights=mu[m2]*u[m2],minlength=nb)
    k=wsum>0
    pu.append(usum[k]/wsum[k]); pw.append(wsum[k])
    uc=np.arange(UEX+du/2,U_max,du)
    pu.append(uc); pw.append(np.exp(uc/2)*du)
    return np.concatenate(pu), np.concatenate(pw)
def solve(h0,m,U_max=80.0):
    pu,pw=points(U_max)
    om=0.5*pw*np.exp(-h0*(pu-us)/2)
    keep=om>1e-14*om.max()
    pu=pu[keep]; om=om[keep]
    N=len(pu)
    xi=2*pu/U_max-1
    Phi=C.chebvander(xi,m)   # N x (m+1)
    phs=C.chebvander(np.array([2*us/U_max-1]),m)[0]
    I=sparse.identity(N,format="csr")
    Ph=sparse.csr_matrix(Phi)
    A=sparse.vstack([sparse.hstack([Ph,-I]),sparse.hstack([-Ph,-I])]).tocsr()
    b=np.zeros(2*N)
    Aeq=sparse.csr_matrix(np.concatenate([phs,np.zeros(N)])[None,:])
    c=np.concatenate([np.zeros(m+1),om])
    bounds=[(None,None)]*(m+1)+[(0,None)]*N
    res=linprog(c,A_ub=A,b_ub=b,A_eq=Aeq,b_eq=[1.0],bounds=bounds,method="highs")
    a=res.x[:m+1]
    return res, a, N
def exactB(h0,a,U_max):
    # exact over n<=1e7, tail integral to 200
    v=np.exp(-h0*(u-us)/2)*C.chebval(2*u/U_max-1,a)
    B1=0.5*np.sum(mu*np.abs(v))
    du=0.002; uc=np.arange(UEX+du/2,200,du)
    vc=np.exp(-h0*(uc-us)/2)*C.chebval(2*uc/U_max-1,a)
    B2=0.5*np.sum(np.exp(uc/2)*np.abs(vc))*du
    vs=C.chebval(2*us/U_max-1,a)
    return B1+B2, B1, B2, vs
out=[]
hs=[float(z) for z in sys.argv[1].split(",")]; ms=[int(z) for z in sys.argv[2].split(",")]
U_max=float(sys.argv[3]) if len(sys.argv)>3 else 80.0
for h0 in hs:
    # single-height reference
    v1=np.exp(-h0*(u-us)/2); Bs=0.5*np.sum(mu*v1)+0.5*np.sum(np.exp(np.arange(UEX,200,0.002)/2)*np.exp(-h0*(np.arange(UEX,200,0.002)-us)/2))*0.002
    print(f"h0={h0} single-height B={Bs:.4g} ratio={ROW/Bs:.4g}",flush=True)
    for m in ms:
        t=time.time()
        res,a,N=solve(h0,m,U_max)
        B,B1,B2,vs=exactB(h0,a,U_max)
        r={"h0":h0,"m":m,"N":N,"status":res.status,"lp":res.fun,"Bexact":B,"Btail":B2,"vstar":vs,"ratio":ROW/B,"sec":round(time.time()-t,1)}
        out.append(r); print(json.dumps(r),flush=True)
json.dump(out,open(f"b1_lp_{sys.argv[1]}_{sys.argv[2]}.json","w"))
print("ROW amplitude", ROW)
