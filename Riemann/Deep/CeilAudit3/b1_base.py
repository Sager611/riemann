
import numpy as np, math, sys, time, json
from scipy.optimize import linprog
from numpy.polynomial import chebyshev as C
x=6e12; L=math.log(x/(4*math.pi)); pr=L/4; us=L/2
th=0.8733458561067185768; y0=2*th-1
ROW=2*pr*math.cosh(pr*y0)
d=np.load("/tmp/ceil3_lambda.npz"); u=d["u"]; mu=d["mu"]
UEX=math.log(1e7)
def lbase(uu,K,a):  # log of base(u) = K log(u/u*) - a (u-u*)
    return K*np.log(uu/us)-a*(uu-us)
def points(ucut=math.log(3000.), bw=0.002, du=0.005, U=250.):
    m1=u<ucut
    pu=[u[m1]]; pw=[mu[m1]]
    m2=~m1
    b=np.floor((u[m2]-ucut)/bw).astype(np.int64)
    wsum=np.bincount(b,weights=mu[m2]); usum=np.bincount(b,weights=mu[m2]*u[m2])
    k=wsum>0
    pu.append(usum[k]/wsum[k]); pw.append(wsum[k])
    uc=np.arange(UEX+du/2,40,du); pu.append(uc); pw.append(np.exp(uc/2)*du)
    uc=np.arange(40+0.025,U,0.05); pu.append(uc); pw.append(np.exp(uc/2)*0.05)
    return np.concatenate(pu), np.concatenate(pw)
PU,PW=points(bw=0.003,du=0.01)
def solve(K,m,eps=0.05,KEEP=1e-11):
    a=max(K/us,eps)+0.5
    lw=np.log(0.5*PW)+lbase(PU,K,a)
    om=np.exp(lw-lw.max()); sc=math.exp(lw.max())
    keep=om>KEEP; pu=PU[keep]; omk=om[keep]
    lo,hi=pu.min(),pu.max()
    xi=lambda z:(2*z-lo-hi)/(hi-lo)
    Phi=C.chebvander(xi(pu),m); phs=C.chebvander(np.array([xi(us)]),m)[0]
    N=len(pu)
    Aeq=np.hstack([(Phi*omk[:,None]).T,-phs[:,None]]); rs=np.abs(Aeq).max(axis=1); Aeq=Aeq/rs[:,None]
    c=np.zeros(N+1); c[-1]=-1
    res=linprog(c,A_eq=Aeq,b_eq=np.zeros(m+1),bounds=[(-1,1)]*N+[(None,None)],method="highs")
    if res.status!=0: return None
    q=-res.eqlin.marginals/rs; q=q/C.chebval(xi(us),q)
    lam=-res.fun*sc
    # exact re-evaluation over all n<=1e7 and the continuous tail to 250 (step 0.001)
    vq=C.chebval(xi(u),q)*np.exp(lbase(u,K,a))
    B1=0.5*np.sum(mu*np.abs(vq))
    uc=np.arange(UEX+0.0005,250,0.001)
    vc=C.chebval(xi(uc),q)*np.exp(lbase(uc,K,a))
    B2=0.5*np.sum(np.exp(uc/2)*np.abs(vc))*0.001
    return dict(K=K,m=m,h0=2*a,N=N,lo=round(lo,2),hi=round(hi,2),Bdual=lam,Bexact=B1+B2,Btail=B2,ratio=ROW/(B1+B2),ratio_dual=ROW/lam)
Ks=[int(z) for z in sys.argv[1].split(",")]; ms=[int(z) for z in sys.argv[2].split(",")]; tag=sys.argv[3]
out=[]
for K in Ks:
    for m in ms:
        t=time.time(); r=solve(K,m)
        if r is None: print(json.dumps({"K":K,"m":m,"fail":1}),flush=True); continue
        r["sec"]=round(time.time()-t,1); out.append(r)
        print(json.dumps({k:(round(v,5) if isinstance(v,float) else v) for k,v in r.items()}),flush=True)
json.dump(out,open(f"b1_base_{tag}.json","w"))
