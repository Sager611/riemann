
import numpy as np, math, sys, json, time
from scipy.optimize import linprog
from numpy.polynomial import chebyshev as C
X=6000000185827.0; LOGT=math.log(X/(4*math.pi)); RHO=LOGT/(4*math.pi); S=1/RHO
TH=0.8733458561067185768; Y0=2*TH-1
pr=math.pi*RHO; us=2*pr; A=2*S; M=200; PT=M*A
d=np.load("/tmp/ceil3_lambda.npz"); u=d["u"]; mu=d["mu"]; UEX=math.log(1e7)
def points(ucut=math.log(3000.), bw=0.003, du=0.01, U=250.):
    m1=u<ucut; pu=[u[m1]]; pw=[mu[m1]]; m2=~m1
    b=np.floor((u[m2]-ucut)/bw).astype(np.int64)
    wsum=np.bincount(b,weights=mu[m2]); usum=np.bincount(b,weights=mu[m2]*u[m2]); k=wsum>0
    pu.append(usum[k]/wsum[k]); pw.append(wsum[k])
    uc=np.arange(UEX+du/2,40,du); pu.append(uc); pw.append(np.exp(uc/2)*du)
    uc=np.arange(40.025,U,0.05); pu.append(uc); pw.append(np.exp(uc/2)*0.05)
    return np.concatenate(pu), np.concatenate(pw)
PU,PW=points()
def design(K,m,uc):
    a=max(K/uc,0.05)+0.5
    lb=lambda z: K*np.log(z/uc)-a*(z-uc)
    lw=np.log(0.5*PW)+lb(PU); om=np.exp(lw-lw.max()); sc=math.exp(lw.max())
    keep=om>1e-11; pu=PU[keep]; omk=om[keep]; lo,hi=pu.min(),pu.max()
    xi=lambda z:(2*z-lo-hi)/(hi-lo)
    Phi=C.chebvander(xi(pu),m); phs=C.chebvander(np.array([xi(uc)]),m)[0]; N=len(pu)
    Aeq=np.hstack([(Phi*omk[:,None]).T,-phs[:,None]]); rs=np.abs(Aeq).max(axis=1); Aeq=Aeq/rs[:,None]
    c=np.zeros(N+1); c[-1]=-1
    res=linprog(c,A_eq=Aeq,b_eq=np.zeros(m+1),bounds=[(-1,1)]*N+[(None,None)],method="highs")
    q=-res.eqlin.marginals/rs; q=q/C.chebval(xi(uc),q)
    v=lambda z: C.chebval(xi(z),q)*np.exp(lb(z))
    B1=0.5*np.sum(mu*np.abs(v(u))); zz=np.arange(UEX+0.0005,250,0.001)
    B2=0.5*np.sum(np.exp(zz/2)*np.abs(v(zz)))*0.001
    return v,B1+B2,(lo,hi)
def env_cue(N, seed):
    rng=np.random.default_rng(seed)
    Am=(rng.standard_normal((N,N))+1j*rng.standard_normal((N,N)))/math.sqrt(2)
    Q,R=np.linalg.qr(Am); Q=Q*(np.diag(R)/abs(np.diag(R)))
    th=np.sort(np.angle(np.linalg.eigvals(Q))%(2*np.pi)); return th/(2*np.pi)*N*S
def config(seed,lc):
    base=(np.arange(M)+0.5)*A
    if seed is None: return base
    env=env_cue(2*M,seed); mids=(env[0::2]+env[1::2])/2
    uu=mids-base; uu-=uu.mean(); U=np.fft.rfft(uu); k=np.arange(len(U))
    U[~((k>0)&(M/np.maximum(k,1)<lc))]=0
    return base+np.fft.irfft(U,n=M)
def response(pos,v,lohi,y=Y0):
    lo,hi=lohi; mlo=int(lo*PT/(4*math.pi)); mhi=int(hi*PT/(4*math.pi))+2
    ms=np.arange(max(mlo,1),mhi); um=4*math.pi*ms/PT
    G=2*np.cosh(um*y/2)*np.exp(-2j*math.pi*np.outer(ms,pos)/PT).sum(axis=1)
    NG=2**int(math.ceil(math.log2(4*mhi))); spec=np.zeros(NG,dtype=complex); spec[ms]=v(um)*G
    D=(2*math.pi/PT)*np.real(np.fft.ifft(spec)*NG)
    return np.max(np.abs(D))
designs=[(60,8),(120,16),(190,16)]
offs=[float(z) for z in sys.argv[1].split(",")]
cfgs=[(None,None)]+[(s,lc) for lc in (8,12,16) for s in (1,2,3)]
out=[]
for (K,m) in designs:
    for off in offs:
        uc=us+off; v,B,lohi=design(K,m,uc)
        row={"K":K,"m":m,"uc_off":off,"Band":round(B,2)}
        for (s,lc) in cfgs:
            pos=config(s,lc); R=response(pos,v,lohi)
            row[f"{s}_{lc}"]=round(R/B,3)
        out.append(row); print(json.dumps(row),flush=True)
json.dump(out,open("b1_cue.json","w"))
