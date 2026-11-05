# disk_trend.py: Thm 5.3 disk |z|=eta|u_s|^2 at m=d-3, k=0 for larger d; also Thm 5.3 own lambda bound
import sys, json, mpmath as mp
sys.path.insert(0,'.')
from saddle_scan import solve, consts
mp.mp.dps=80
g=[mp.mpf(s) for s in json.load(open(sys.argv[1]))]
for d in [9,13,20,25,30,35,40,50,60,70,80,90]:
    k=0; N=d-1
    coeffs=[mp.binomial(N,j)*g[k+1+j] for j in range(N+1)]
    r=mp.polyroots(coeffs[::-1], maxsteps=800, extraprec=400)
    xs=sorted([mp.re(x) for x in r], reverse=True)
    M=4*d+2*k-1
    out=[]
    for m in [1, d-4, d-3]:
        z=mp.sqrt(-2*xs[m-1]); om=float(z*mp.sqrt(M)/2)
        u,res,Fu,S,W,V,gg,eta,nu2=solve(om,M,k)
        c=consts(u,Fu,S,W,V,gg,eta,nu2,om,M,k)
        out.append((m,round(om,3),round(c['sigma'],4),round(c['tau'],4),round(c['alpha'],4),'%.1e'%res))
    lam=(1-4.38167*(8*d-2)**(-2/3))**2
    print('PASS' if all(o[4]<=0.7945 for o in out) else 'FAIL', 'd',d,'(m,omega,sigma,tau,alpha=eta|us|^2,res)',out,' paper-chain lam_max=%.4f -> 1-(1-lam)(1-0.45lam)=%.4f'%(lam,1-(1-lam)*(1-0.45*lam)))
