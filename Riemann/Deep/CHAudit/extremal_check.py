
# extremal_check.py: numerical check of the band-limited approximants used by CH.
# (A) 2511.14736 Prop carlitfou: phi_lambda=1_[-1,1]*Phi_lambda, ||phihat - I_lambda||_1 = tanh(lambda/4)/lambda, decay kappa=1/(16 pi)
# (B) 2512.15709 Prop grahamvaaler: phi^{+-}_nu, one-sidedness F+>=I>=F-, L1 errors 1/(1-e^-nu)-1/nu and 1/nu-1/(e^nu-1).
import numpy as np
from numpy.polynomial.legendre import leggauss
x,w=leggauss(1500)
tp=(x+1)/2; wp=w/2
t=np.concatenate([-tp,tp]); wt=np.concatenate([wp,wp])
Y=150.0; h=0.005
y=np.round(np.arange(-Y,Y+h/2,h)/h)*h; i0=int(np.argmin(np.abs(y))); y[i0]=0.0
def fhat(vals):
    out=np.empty(len(y),dtype=complex)
    for k in range(0,len(y),400):
        yy=y[k:k+400]; out[k:k+400]=np.exp(-2j*np.pi*np.outer(yy,t))@(wt*vals)
    return out
def I(lam):
    return np.where(np.sign(lam)*y>0, np.exp(np.clip(-lam*y,-700,700)), 0.0)
def l1(F,lam):
    Il=I(lam); Ir=I(lam).copy()
    if lam>0: Ir[i0]=1.0
    else: Il[i0]=1.0
    return np.trapezoid(np.abs(F[:i0+1]-Il[:i0+1]),y[:i0+1])+np.trapezoid(np.abs(F[i0:]-Ir[i0:]),y[i0:])
for lam in [0.5,2.0,5.0,-1.0,-0.01]:
    Phi=(np.sign(lam)/4)*(1/np.tanh(np.pi*t/(2j)+lam/4)-np.tanh(lam/4))
    F=fhat(Phi); d=F-I(lam); big=np.abs(y)>=20
    C=np.max(y[big]**2*np.abs(d[big])); tail=2*C/Y
    print('CL lam=%6.2f: max|Im F|=%.1e  L1[-150,150]=%.6f (+tail<=%.1e)  claim tanh(l/4)/l=%.6f  sup_{|y|>=20} y^2|F-I|=%.5f (kappa=1/16pi=%.5f)'%(lam,np.max(np.abs(F.imag)),l1(F,lam),tail,np.tanh(lam/4)/lam,C,1/(16*np.pi)),flush=True)
for nu in [0.5,2.0,5.0]:
    for sgn,name in [(1,'+'),(-1,'-')]:
        wz=-2j*np.pi*t+nu
        Pc=0.5*(1/np.tanh(wz/2)+sgn*1)
        Ps=(1j/(2*np.pi))*((nu/2)/np.tanh(nu/2)-(wz/2)/np.tanh(wz/2)+sgn*np.pi*1j*t)
        F=fhat(Pc+np.sign(t)*Ps); d=F.real-I(nu); nz=np.abs(y)>h/2
        claim=(1/(1-np.exp(-nu))-1/nu) if sgn>0 else (1/nu-1/(np.exp(nu)-1))
        big=np.abs(y)>=20; tail=2*np.max(y[big]**2*np.abs(d[big]))/Y
        print('GV nu=%.1f %s: max|Im F|=%.1e  min(F-I)=%.2e max(F-I)=%.2e (y!=0)  F(0)=%.6f  L1=%.6f (+tail<=%.1e) claim=%.6f'%(nu,name,np.max(np.abs(F.imag)),np.min(d[nz]),np.max(d[nz]),F.real[i0],l1(F.real,nu),tail,claim),flush=True)
