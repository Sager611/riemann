
import numpy as np, math, sys, time, json
from scipy.optimize import linprog
from numpy.polynomial import chebyshev as C
x=6e12; L=math.log(x/(4*math.pi)); pr=L/4; us=L/2
th=0.8733458561067185768; y0=2*th-1
ROW=2*pr*math.cosh(pr*y0)
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
def solve_dual(h0,m,U_max):
    pu,pw=points(U_max)
    om=0.5*pw*np.exp(-h0*(pu-us)/2)
    keep=om>KEEP*om.max(); pu=pu[keep]; om=om[keep]; N=len(pu)
    Phi=C.chebvander(2*pu/U_max-1,m)          # N x (m+1)
    phs=C.chebvander(np.array([2*us/U_max-1]),m)[0]
    # y_i = om_i z_i, |z_i|<=1 ; maximize lam
    Aeq=np.hstack([(Phi*om[:,None]).T, -phs[:,None]])
    rs=np.abs(Aeq).max(axis=1); Aeq=Aeq/rs[:,None]
    c=np.zeros(N+1); c[-1]=-1
    bounds=[(-1,1)]*N+[(None,None)]
    res=linprog(c,A_eq=Aeq,b_eq=np.zeros(m+1),bounds=bounds,method="highs")
    if res.status!=0:
        print("LP status",res.status,res.message,flush=True); return res,None,N
    a=-res.eqlin.marginals/rs
    return res, a, N
def exactB(h0,a,U_max):
    vs=C.chebval(2*us/U_max-1,a); a=a/vs
    v=np.exp(-h0*(u-us)/2)*C.chebval(2*u/U_max-1,a)
    B1=0.5*np.sum(mu*np.abs(v))
    du=0.002; uc=np.arange(UEX+du/2,250,du)
    vc=np.exp(-h0*(uc-us)/2)*C.chebval(2*uc/U_max-1,a)
    B2=0.5*np.sum(np.exp(uc/2)*np.abs(vc))*du
    return B1+B2, B2, a
hs=[float(z) for z in sys.argv[1].split(",")]; ms=[int(z) for z in sys.argv[2].split(",")]
U_max=float(sys.argv[3]); tag=sys.argv[4]; KEEP=float(sys.argv[5]) if len(sys.argv)>5 else 1e-12
out=[]
for h0 in hs:
    for m in ms:
        t=time.time()
        res,a,N=solve_dual(h0,m,U_max)
        if a is None: continue
        lam=-res.fun
        B,B2,a=exactB(h0,a,U_max)
        r={"h0":h0,"m":m,"N":N,"st":res.status,"Bmin_dual":lam,"Bexact":B,"Btail":B2,"ratio_dual":ROW/lam,"ratio_exact":ROW/B,"sec":round(time.time()-t,1)}
        out.append(r); print(json.dumps(r),flush=True)
        np.save(f"/tmp/ceil3_a_{h0}_{m}.npy",a)
json.dump(out,open(f"b1_dual_{tag}.json","w"))
