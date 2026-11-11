
import numpy as np
from flint import acb, ctx
d=np.load('/Users/adriansager/playground/riemann/.scratch/chaudit/zp.npy'); g=d[0]
n=np.arange(1,len(g)+1)
F=(g/(2*np.pi))*np.log(g/(2*np.pi*np.e))+7/8
Qp=n-F; Qm=n-1-F
b=np.log(g)/5+2
print('lem:oldie check on zeros<=%.0f: max Q/(log t/5+2)=%.4f, min Q/(log t/5+2)=%.4f ; max|Q| for t<=280: %.3f'%(g[-1],max(Qp.max(),Qm.max()) and np.max(np.maximum(Qp,Qm)/b),np.min(np.minimum(Qp,Qm)/b),np.max(np.abs(np.concatenate([Qp[g<=280],Qm[g<=280]])))))
print('cor:brut check: max N(t)/((t/2pi)log(t/2pi)) at zeros = %.6f'%np.max(n/((g/(2*np.pi))*np.log(g/(2*np.pi)))))
ctx.prec=96
for T in [1e6, 1e6+1]:
    best=1e9; bs=None
    for k in range(2401):
        s=-1/64+(1+1/64)*k/2400; z=float(abs(acb(s,T).zeta()).mid())
        if z<best: best,bs=z,s
    print('chamlu T=%.0f: max 1/|zeta(s+iT)| on grid = %.9f at s=%.5f (paper for T=1e6+1: 0.894198297)'%(T,1/best,bs))
