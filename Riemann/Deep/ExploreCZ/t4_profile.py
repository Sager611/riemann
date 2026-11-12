
import numpy as np
lam=np.log(2.0); c=3*np.pi
def F(x):
    y=lam*x/2
    return np.where(np.abs(y)<1e-12,1.0,(np.sin(y)/np.where(np.abs(y)<1e-12,1.0,y))**2)
u=np.linspace(-3000,3000,1200001); du=u[1]-u[0]
IF={k:float(np.sum(F(u)**k)*du) for k in (1,2,3)}
rho={k:1.13*IF[k] for k in IF}
def run(N,tau,s,K=3,dt=0.05):
    T=N**tau; V2=N**(2*s); Gk={k:float(N**k) for k in range(1,K+1)}   # minimal G_k (hardest case for the profile)
    Bk={k:(T+c*2**k*N**k)*Gk[k]/(rho[k]*V2**k) for k in Gk}
    B=min(Bk.values()); Lp=T/B; R=int(B)
    ok_single=all(2*rho[k]*V2**k<=c*2**k*N**k*Gk[k] for k in Gk)
    t=np.arange(-200,T+200,dt); w=np.zeros_like(t)
    for r in range(1,R+1): w+=V2*F(t-r*Lp)
    worst=0.0
    for k in Gk:
        cs=np.concatenate([[0],np.cumsum(w**k)*dt])
        for H in np.unique(np.round(np.geomspace(1,T+400,60)/dt).astype(int)):
            if H>=len(t): continue
            mass=(cs[H:]-cs[:-H]).max()
            worst=max(worst, mass/(H*dt*Gk[k]+c*2**k*N**k*Gk[k]))
    # 1-separated points with w >= V2
    idx=np.where(w>=V2*(1-1e-9))[0]; cnt=0; last=-1e9
    for i in idx:
        if t[i]-last>=1.0: cnt+=1; last=t[i]
    hux=N**(2-2*s)+T*N**(4-6*s); gm=N**(2-2*s)+N**(3.6-4*s)+T*N**(2.4-4*s); mvt=(T+N)*N**(1-2*s)
    print(f"N={N} tau={tau} s={s}: B_k={ {k:round(v,1) for k,v in Bk.items()} } L'={Lp:.1f} single-spike ok={ok_single} "
          f"max(window mass/(|I|G_k+c2^kN^kG_k))={worst:.3f} (<=1 admissible) R(w)={cnt} | MVT~{mvt:.0f} Hux~{hux:.0f} GM~{gm:.0f} (no constants)")
print("int F^k:",{k:round(v,4) for k,v in IF.items()})
for (N,tau,s) in [(10,3.0,0.8),(30,2.0,0.8),(100,1.5,0.8),(300,1.2,0.8)]:
    run(N,tau,s)
