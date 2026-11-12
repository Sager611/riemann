
import numpy as np, json
from fractions import Fraction as Fr
z=np.load('/Users/adriansager/playground/riemann/Riemann/Deep/ExploreHelio/results/zeta_zeros.npy'); z=np.sort(z[z>100])
dh=json.load(open('/Users/adriansager/playground/riemann/Riemann/Deep/ExploreHelio/results/dh_zeros.json'))
def unfold(g): return g/(2*np.pi)*np.log(g/(2*np.pi*np.e))
x=unfold(z); gaps=np.diff(x)
print('zeta: %d zeros in [%.0f, %.0f]; mean unfolded gap %.4f'%(len(z),z[0],z[-1],gaps.mean()))
for lam in [0.768,0.77,0.775,0.7851]: print('  D(%.4f) zeta = %.4f'%(lam,np.mean(gaps<=lam)))
# DH
try:
    if isinstance(dh,dict):
        keys=list(dh.keys()); print('dh keys',keys[:10])
        arr=None
        for k in keys:
            v=dh[k]
            if isinstance(v,list) and len(v)>50: arr=v; print('using',k,len(v)); break
    else: arr=dh
    a=np.array([complex(*e) if isinstance(e,(list,tuple)) else complex(e) if not isinstance(e,dict) else complex(e.get('re',e.get('beta',0.5)),e.get('im',e.get('gamma'))) for e in arr])
    on=a[np.abs(a.real-0.5)<1e-6].imag; on=np.sort(on[on>20])
    gd=np.diff(unfold(on)); print('DH on-line zeros used: %d, mean gap %.3f'%(len(on),gd.mean()))
    for lam in [0.768,0.7851]: print('  D(%.4f) DH on-line = %.4f'%(lam,np.mean(gd<=lam)))
except Exception as e: print('DH parse failed',repr(e)[:200])
# certificate f-hat on real zeros
d=json.load(open('/Users/adriansager/playground/riemann/.scratch/explore40/logs/cert768.json')); B=[float(Fr(b)) for b in d['B']]; B[0]+=1e-3; s=1.1; k=np.sqrt(2*np.pi); D=d['D']
def fh(u):
    v=k*s*np.asarray(u,float); H=[np.ones_like(v),2*v]
    for m in range(1,2*D): H.append(2*v*H[-1]-2*m*H[-2])
    return s*sum(B[n]*(-1)**n*H[2*n] for n in range(D+1))*np.exp(-v*v/2)
N=len(x); tot=0.0
for j in range(1,40):
    du=x[j:]-x[:-j]; tot+=2*fh(du).sum()
off=tot/N; print('real zeros: (1/N) sum_{j!=l} fhat = %.4f ; + fhat(0) = %.4f ; J_theta(f) = 0.3056 (input bound), hard-core bound fhat(0)+2c = %.4f'%(off, off+fh([0.0])[0], fh([0.0])[0]+2*(-0.347862)))
# Z_0 of real zeros: drop zeros with a neighbouring gap <= lam
lam=0.768; small=gaps<=lam; bad=np.zeros(N,bool); bad[:-1]|=small; bad[1:]|=small
x0=x[~bad]; g0=np.diff(x0); print('Z_0 (lam=0.768): N_0/N = %.3f, mean gap %.3f, min gap %.3f'%(len(x0)/N, g0.mean(), g0.min()))
e=(fh(g0[:-1])+fh(g0[:-1]+g0[1:])).mean(); print('  energy per point (1/N_0) sum_i [fhat(g_i)+fhat(g_i+g_{i+1})] = %.4f  (bound c_cert=-0.3479 is for density-1 hard-core sets; Z_0 here has density %.3f)'%(e, 1/g0.mean()))
