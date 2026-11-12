
import numpy as np, sys
base='/Users/adriansager/playground/riemann/.scratch/explore33/'
def load(f):
    d=np.loadtxt(base+f); return d[:,0], d[:,1]
out=[]
for f,name in [('curve_0.txt','primes'),('curve_1.txt','p+1'),('curve_2.txt','smooth'),('curve_3.txt','Poisson'),('curve_4.txt','iid eps=.5'),('curve_5.txt','shift e^-.25'),('curve_6.txt','planted .8+10i')]:
    x,M=load(f)
    sel=x>=1e4
    x=x[sel]; M=M[sel]; L=np.log(x)
    # growth exponent: max |M| over half-decades, regress
    edges=10**np.arange(4,8.01,0.5)
    mx=[]; xc=[]
    for a,b in zip(edges[:-1],edges[1:]):
        s=(x>=a)&(x<b)
        if s.sum(): mx.append(np.abs(M[s]).max()); xc.append(np.sqrt(a*b))
    mx=np.array(mx); xc=np.array(xc)
    slope=np.polyfit(np.log(xc[-4:]),np.log(mx[-4:]),1)[0]
    slope_all=np.polyfit(np.log(xc),np.log(mx),1)[0]
    line=f"{name:16s} max|M|/sqrt(x) on [1e4,1e8]={np.max(np.abs(M)/np.sqrt(x)):.3f}  growth exponent (half-decade maxima, last 2 decades)={slope:.3f} (all 4 decades {slope_all:.3f})"
    if name.startswith('planted'):
        s=x>=1e5
        y=M[s]/x[s]**0.8; g=10.0
        A=np.vstack([np.cos(g*L[s]),np.sin(g*L[s])]).T
        coef,res,_,_=np.linalg.lstsq(A,y,rcond=None)
        fit=A@coef; r2=1-np.sum((y-fit)**2)/np.sum((y-y.mean())**2)
        amp=np.hypot(*coef); ph=np.arctan2(-coef[1],coef[0])
        line+=f"\n   fit M/x^0.8 = a cos(10 log x) + b sin(10 log x): amp={amp:.4f} phase={ph:.3f} R^2={r2:.4f}"
        # scan beta and gamma to show optimum
        best=None
        for beta in np.arange(0.70,0.901,0.01):
            for gg in np.arange(9.5,10.51,0.05):
                y2=M[s]/x[s]**beta; A2=np.vstack([np.cos(gg*L[s]),np.sin(gg*L[s])]).T
                c2,_,_,_=np.linalg.lstsq(A2,y2,rcond=None); f2=A2@c2
                rr=1-np.sum((y2-f2)**2)/np.sum((y2-y2.mean())**2)
                # amplitude stability: fit separately on first and second half
                if best is None or rr>best[0]: best=(rr,beta,gg)
        line+=f"\n   best (beta,gamma) by R^2 over grid: beta={best[1]:.2f} gamma={best[2]:.2f} R^2={best[0]:.4f}"
        # amplitude drift test: fit on [1e5,1e6.5] and [1e6.5,1e8]
        for (a,b) in [(1e5,10**6.5),(10**6.5,1e8)]:
            s2=(x>=a)&(x<b); y2=M[s2]/x[s2]**0.8; A2=np.vstack([np.cos(10*L[s2]),np.sin(10*L[s2])]).T
            c2,_,_,_=np.linalg.lstsq(A2,y2,rcond=None); line+=f"\n   window [{a:.1e},{b:.1e}]: amp={np.hypot(*c2):.4f}"
    if name=='primes':
        s=x>=1e5; y=M[s]/x[s]**0.8; A=np.vstack([np.cos(10*L[s]),np.sin(10*L[s])]).T
        coef,_,_,_=np.linalg.lstsq(A,y,rcond=None); line+=f"\n   same fit for primes: amp={np.hypot(*coef):.5f}"
    out.append(line)
print("\n".join(out))
