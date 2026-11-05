# omega_np.py: fast double-precision check of Thm 4.1 Part 2 (Omega_r strictly increasing at m=d-1), k=0
import sys, json, numpy as np
from numpy.polynomial.hermite_e import hermeroots
pi=np.pi
def Phi_derivs(u):
    # Phi, Phi', Phi'' for Chung kernel, vectorized
    P=np.zeros_like(u); P1=np.zeros_like(u); P2=np.zeros_like(u)
    for n in range(1,8):
        a=pi*n*n*np.exp(2*u)
        e=np.exp(-a+u/2)
        # phi=2 a(2a-3) e^{u/2} e^{-a}; write phi=2*h(u)*e with h=a(2a-3)=2a^2-3a, a'=2a
        h=2*a*a-3*a; h1=(4*a-3)*2*a; h2=(4*2*a)*2*a+(4*a-3)*4*a
        L1=0.5-2*a; L2=-4*a   # (ln e)' and ''
        phi=2*h*e; phi1=2*e*(h1+h*L1); phi2=2*e*(h2+2*h1*L1+h*(L1*L1+L2))
        P+=phi; P1+=phi1; P2+=phi2
    return P,P1,P2
def kappa0p(u):
    # full-kernel formula for u<=1 (no underflow there); n=1 rational identity for u>1 (n>=2 modes < e^{-3 pi e^2} there)
    out=np.empty_like(u)
    lo=u<=1.0
    P,P1,P2=Phi_derivs(u[lo])
    out[lo]=-(P2*u[lo]*P-P1*(P+u[lo]*P1))/((u[lo]*P)**2)
    uu=u[~lo]; a=pi*np.exp(2*uu)
    out[~lo]=(2*a*(2*uu-1)+4.5+6*(4*a*uu+2*a-3)/(2*a-3)**2)/uu**2
    return out
def He(n,y):
    a=np.ones_like(y); b=y.copy()
    for k in range(1,n): a,b=b,y*b-k*a
    return b
if __name__=='__main__':
    R=json.load(open(sys.argv[1]))
    res=[]
    for d in range(4,31):
        key=f"{d},0"
        if key not in R: continue
        for m in [d-1]:
            z=float(R[key]['z'][m-1]); n=2*d-1
            U=np.linspace(1e-6,4.0,400001)
            P,_,_=Phi_derivs(U)
            gvals=U*P*He(n,z*U/np.sqrt(2))
            # cumulative from right: W(u)=int_u^inf g
            h=U[1]-U[0]
            cum=np.concatenate([[0],np.cumsum((gvals[1:]+gvals[:-1])*h/2)])
            W=cum[-1]-cum
            kp=kappa0p(U)
            sgn=np.sign(W); idx=np.where(sgn[1:]*sgn[:-1]<0)[0]
            # zeros tau (exclude u~0 where W(0)=f=0 numerically)
            idx=[i for i in idx if U[i]>0.02]
            bounds=[0]+[i+1 for i in idx]+[len(U)-1]
            Om=[np.trapezoid(np.abs(W[bounds[r]:bounds[r+1]+1])*kp[bounds[r]:bounds[r+1]+1],U[bounds[r]:bounds[r+1]+1]) for r in range(len(bounds)-1)]
            inc=all(Om[i]<Om[i+1] for i in range(len(Om)-1))
            res.append((d,m,len(idx),inc))
            if not inc or d in (9,20,30): print('d',d,'m',m,'#interior zeros of W',len(idx),'(expect d-2=%d)'%(d-2),'Omega increasing:',inc,'Omega ratios', [round(Om[i+1]/Om[i],3) for i in range(len(Om)-1)][:12])
    print('Thm 4.1(2) Omega increasing at m=d-1, k=0, d=4..30:', 'PASS' if all(r[3] for r in res) and all(r[2]==r[0]-2 for r in res) else 'FAIL')
