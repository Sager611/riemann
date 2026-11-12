
# float verification of a primal_2gap solution on fine grids (A linearly interpolated, constant beyond U)
import numpy as np, sys, json
import importlib.util
spec=importlib.util.spec_from_file_location('pg','primal_2gap.py'); pg=importlib.util.module_from_spec(spec); spec.loader.exec_module(pg)
o=json.load(open(sys.argv[1])); theta=eval(sys.argv[2]) if len(sys.argv)>2 else 399993/700000
D=o['D']; s=o['s']; lam=o['lam']; U=o['U']; cvec=np.array(o['coef'])
f=lambda a: pg.basis_f(D,s,a)@cvec
fh=lambda u: pg.basis_fh(D,s,u)@cvec
a=np.arange(0,14,0.0005); fa=f(a); print('min f on [0,14] = %.3e (f(0)=%.4f, max f=%.4f) at a=%.4f'%(fa.min(),fa[0],fa.max(),a[np.argmin(fa)]))
x,wq=np.polynomial.legendre.leggauss(400); aq=(x+1)*theta/2; wq=wq*theta/2
J=f(np.array([0.0]))[0]+2*np.sum(f(aq)*wq*aq); print('J=%.6f (LP %.6f), fhat(0)=%.6f'%(J,o['J'],fh(np.array([0.0]))[0]))
lo = lam if not o['aug'] else U
u=np.arange(lo,14,0.0005); fu=fh(u); print('max fhat on [%.4f,14] = %.3e at u=%.4f'%(lo,fu.max(),u[np.argmax(fu)]))
if o['aug']:
    g=np.array(o['g']); A=np.array(o['A']); c=o['c']; beta=o['beta']
    Ai=lambda x: np.interp(x,g,A,right=A[-1])
    step=float(sys.argv[3]) if len(sys.argv)>3 else 0.002
    gg=np.arange(lam, U+0.3, step); G1=np.where(gg<U, fh(gg), 0.0)
    worst=-1e9; arg=None
    for i,g1 in enumerate(gg):
        s12=g1+gg; G12=np.zeros_like(gg); m=s12<U
        if m.any(): G12[m]=fh(s12[m])
        lhs=G1[i]+G12+Ai(g1)-Ai(gg)+beta*(1-g1)-c
        j=np.argmax(lhs)
        if lhs[j]>worst: worst=lhs[j]; arg=(g1,gg[j])
    print('2D: max LHS-c = %.3e at (g1,g2)=(%.4f,%.4f), step %.4f; c=%.5f beta=%.4f'%(worst,arg[0],arg[1],step,c,beta))
    ceff=c+max(worst,0)
    print('margin LP = %.6f ; margin with c_eff = %.6f (needs >0, plus f>=0 / fhat<=0 slack)'%(o['J']-1-2*c, J-fh(np.array([0.0]))[0]-2*ceff))
else:
    print('margin = %.6f'%(J-1))
