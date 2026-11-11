
# zeros2.py t_lo t_hi out.npy : all zeros of Z(t) in (t_lo, t_hi] by Riemann-Siegel (C0+C1) on a 0.05 grid, plus a
# local-|Z|-minimum search for close pairs inside one grid cell; |zeta'| = |Z'| at each zero. Non-rigorous.
import sys, numpy as np, mpmath as mp, time, resource
from numpy.polynomial import chebyshev as C
t00=time.time(); tlo=float(sys.argv[1]); thi=float(sys.argv[2]); out=sys.argv[3]
mp.mp.dps=50
def Psi(p):
    p=mp.mpf(p); return mp.cos(2*mp.pi*(p*p-p-mp.mpf(1)/16))/mp.cos(2*mp.pi*p)
K=120; nodes=[(mp.cos(mp.pi*(k+mp.mpf(1)/2)/K)+1)/2 for k in range(K)]
cf=C.chebfit(np.array([float(2*x-1) for x in nodes]),np.array([float(Psi(x)) for x in nodes]),80)
d1=C.chebder(cf)*2.0; d3=C.chebder(cf,3)*8.0; d4=C.chebder(cf,4)*16.0
def th(t): return 0.5*t*np.log(t/(2*np.pi)) - 0.5*t - np.pi/8 + 1/(48*t) + 7/(5760*t**3)
def thp(t): return 0.5*np.log(t/(2*np.pi)) - 1/(48*t*t) - 7/(1920*t**4)
def Z(t):
    if len(t)>500: return np.concatenate([Z(t[k:k+500]) for k in range(0,len(t),500)])
    if len(t)==0: return np.zeros(0)
    a=np.sqrt(t/(2*np.pi)); N=np.floor(a).astype(np.int64); p=a-N
    n=np.arange(1,int(N.max())+1,dtype=np.float64); ln=np.log(n); isq=1/np.sqrt(n)
    ph=th(t)[:,None]-t[:,None]*ln[None,:]
    s=2*np.sum(np.where(n[None,:]<=N[:,None], isq[None,:]*np.cos(ph),0.0),axis=1)
    sg=np.where(N%2==1,1.0,-1.0); x=2*p-1
    return s+sg*a**-0.5*(C.chebval(x,cf)-C.chebval(x,d3)/(96*np.pi**2)/a)
def Zp(t):
    if len(t)>500: return np.concatenate([Zp(t[k:k+500]) for k in range(0,len(t),500)])
    if len(t)==0: return np.zeros(0)
    a=np.sqrt(t/(2*np.pi)); N=np.floor(a).astype(np.int64); p=a-N
    n=np.arange(1,int(N.max())+1,dtype=np.float64); ln=np.log(n); isq=1/np.sqrt(n)
    ph=th(t)[:,None]-t[:,None]*ln[None,:]
    s=-2*np.sum(np.where(n[None,:]<=N[:,None], isq[None,:]*np.sin(ph)*(thp(t)[:,None]-ln[None,:]),0.0),axis=1)
    sg=np.where(N%2==1,1.0,-1.0); x=2*p-1
    return s+sg*(a**-0.5*C.chebval(x,d1)-0.5*a**-1.5*C.chebval(x,cf)-a**-1.5*C.chebval(x,d4)/(96*np.pi**2))/(4*np.pi*a)
def bisect(a,b,fa):
    for it in range(17):
        m=0.5*(a+b); fm=Z(m); left=np.sign(fm)==np.sign(fa)
        a=np.where(left,m,a); fa=np.where(left,fm,fa); b=np.where(left,b,m)
    fb=Z(b); return a-fa*(b-a)/(fb-fa)
h=0.05; seg=2000.0; zs=[]; npairs=0; ncand=0
tstart=tlo
while tstart<thi:
    tend=min(tstart+seg,thi)
    grid=np.arange(tstart,tend+h/2,h); grid=grid[grid<=tend]
    if grid[-1]<tend: grid=np.append(grid,tend)
    v=Z(grid); sv=np.sign(v)
    sc=np.where(sv[:-1]!=sv[1:])[0]
    roots=[bisect(grid[sc],grid[sc+1],v[sc])]
    av=np.abs(v); i=np.arange(1,len(v)-1)
    cand=i[(av[i]<=av[i-1])&(av[i]<=av[i+1])&(sv[i-1]==sv[i])&(sv[i]==sv[i+1])]
    ncand+=len(cand)
    if len(cand):
        A=grid[cand-1].copy(); B=grid[cand+1].copy(); s=sv[cand]; gr=(np.sqrt(5)-1)/2
        for it in range(24):
            c=B-gr*(B-A); d=A+gr*(B-A); fc=s*Z(c); fd=s*Z(d)
            lft=fc<fd; B=np.where(lft,d,B); A=np.where(lft,A,c)
        tm=0.5*(A+B); fm=s*Z(tm); hit=fm<0
        if hit.any():
            npairs+=int(hit.sum()); a0=grid[cand[hit]-1]; b0=grid[cand[hit]+1]; tmh=tm[hit]
            roots.append(bisect(a0,tmh,Z(a0))); roots.append(bisect(tmh,b0,Z(tmh)))
    r=np.sort(np.concatenate(roots)); zs.append(np.vstack([r,np.abs(Zp(r))]))
    tstart=tend
    rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20
    if len(zs)%25==0: print('progress t=%.0f zeros %d pairs %d cand %d maxrss %.0f MB %.0fs'%(tend,sum(q.shape[1] for q in zs),npairs,ncand,rss,time.time()-t00),flush=True)
    if rss>250: print('ABORT memory'); raise SystemExit(1)
Zall=np.concatenate(zs,axis=1)
np.save(out,Zall)
print('done (%.4f,%.4f]: zeros %d, close pairs recovered %d, candidates %d, %.0fs'%(tlo,thi,Zall.shape[1],npairs,ncand,time.time()-t00))
