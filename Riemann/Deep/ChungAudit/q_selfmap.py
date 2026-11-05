# q_selfmap.py: Lemma A.4 Part 2a (p.55-56) claims T(Q) subset Q and a unique fixed point W_s=T(W_s) in Q.
# Test: is our saddle (root of printed F) in Q? Does T map sample points of Q into Q? (k=0 cases)
import sys, cmath, math, json, numpy as np
sys.path.insert(0,'.')
from saddle_scan import solve
pi=math.pi
def T(W,om,M,k,d):
    u=0.5*cmath.log(W/(2*pi))
    g=om*om/M; eta=om*om/(M*M); nu2=k*(k+1)/(om*om)
    V=1-eta*u*u-nu2/(u*u); Vp=-2*eta*u+2*nu2/u**3
    e2=cmath.exp(-2*u); e4=cmath.exp(-4*u)
    Q=1-15/(4*pi)*e2+15/(8*pi*pi)*e4; Qp=15/(2*pi)*e2-15/(2*pi*pi)*e4
    S=cmath.sqrt(V)
    if S.real<0: S=-S
    E=Qp/Q-Vp/(4*V)
    return 6.5+k/u+g*u+1j*om*S+E
def inQ(W,om,M,k,d):
    u=0.5*cmath.log(W/(2*pi)); g=om*om/M
    return (W.real>=6.40+0.65*k*u.real/abs(u)**2+g/2) and (W.real<=2*(d+k)+9.5) and (0<W.imag<=1.15*om) and (abs(W)>=max(14.49,0.70*om))
cases=[(4,5,1,16.5664741),(8,1,1,13.7020061),(22,0,19,58.0785),(30,0,27,75.3386)]
for d,k,m,om in cases:
    M=4*d+2*k-1
    u,res,Fu,S,W,V,g,eta,nu2=solve(om,M,k)
    # sample Q on a grid and test T(W) in Q
    g=om*om/M
    ReLo=6.40+g/2; ReHi=2*(d+k)+9.5
    tot=0; out=0; ex=None
    if ReHi>ReLo:
        for a in np.linspace(ReLo,ReHi,25):
            for b in np.linspace(0.01,1.15*om,25):
                Wq=complex(a,b)
                if not inQ(Wq,om,M,k,d): continue
                tot+=1
                Tw=T(Wq,om,M,k,d)
                if not inQ(Tw,om,M,k,d):
                    out+=1
                    if ex is None: ex=(round(Wq.real,2),round(Wq.imag,2),round(Tw.real,2),round(Tw.imag,2))
    print(('PASS ' if (inQ(W,om,M,k,d) and out==0) else 'FAIL ')+'(d,k,m)=(%d,%d,%d) omega=%.4f: saddle W_s=%.3f%+.3fi in Q: %s | Q-box Re in [%.2f,%.2f]; sampled %d points of Q, T(W) outside Q for %d (e.g. W->T(W): %s)'%(d,k,m,om,W.real,W.imag,inQ(W,om,M,k,d),ReLo,ReHi,tot,out,ex))
