#!/usr/bin/env python3
"""alptoy item 4: Gz(X) = sum_{n<=X} Lambda(n)(X/n - 1) - X log X + (1+gamma) X  (Formalism/StatementOC.lean:11-13)
for all real X in [1, XMAX] via a segmented sieve. Gz is continuous, concave between prime powers:
on [n_j, n_{j+1}) Gz(X) = X*L_j - psi_j - X log X + (1+g)X, L_j = sum_{n<=n_j} Lambda(n)/n.
Extrema per interval: endpoints and X* = exp(L_j + g). Records per dyadic bin: max/min of (Gz - log 2pi)/sqrt X.
RH prediction (explicit formula): Gz(X) = log 2pi - sum_k X^-2k/(2k(2k+1)) + sum_rho X^rho/(rho(1-rho)),
so |Gz - log 2pi + r(X)| <= (2 + g - log 4pi) sqrt X = 0.0461914 sqrt X."""
import numpy as np, math, sys, json, time
XMAX=int(float(sys.argv[1])); OUT=sys.argv[2]
g=0.57721566490153286061; L2P=math.log(2*math.pi); C=2+g-math.log(4*math.pi)
r=int(math.isqrt(XMAX))+1
bs=np.ones(r+1,dtype=bool); bs[:2]=False
for i in range(2,int(r**.5)+1):
    if bs[i]: bs[i*i::i]=False
bp=np.nonzero(bs)[0]
pp=[]  # prime powers p^j, j>=2
for p in bp.tolist():
    v=p*p
    while v<=XMAX: pp.append((v,math.log(p))); v*=p
pp.sort(); ppn=np.array([v for v,_ in pp],dtype=np.int64); ppl=np.array([l for _,l in pp])
SEG=1<<22
psi_b=0.0; L_b=0.0   # bases (sums up to previous segment)
carry=None  # (n_last, psi_last, L_last) from previous segment
bins={}  # dyadic bin -> [max ratio, argmax, min ratio, argmin, maxabs Gz/sqrtX, count_pos, count]
def upd(X,G):
    # X,G arrays: record per dyadic bin
    ratio=(G-L2P)/np.sqrt(X); ab=np.abs(G)/np.sqrt(X)
    b=np.floor(np.log2(X)).astype(np.int64)
    for bb in np.unique(b):
        m=b==bb; rr=ratio[m]; xx=X[m]
        i1=int(np.argmax(rr)); i2=int(np.argmin(rr))
        e=bins.setdefault(int(bb),[-1e9,0,1e9,0,0.0,0,0])
        if rr[i1]>e[0]: e[0]=float(rr[i1]); e[1]=float(xx[i1])
        if rr[i2]<e[2]: e[2]=float(rr[i2]); e[3]=float(xx[i2])
        e[4]=max(e[4],float(ab[m].max())); e[5]+=int((rr>0).sum()); e[6]+=int(m.sum())
t0=time.time(); signchg=0; lastsign=None; samples={}
SAMP=[10**j for j in range(1,10) if 10**j<=XMAX]+[int(XMAX)]
lo=1
while lo<=XMAX:
    hi=min(lo+SEG,XMAX+1)
    seg=np.ones(hi-lo,dtype=bool)
    if lo<=1: seg[:2-lo]=False
    for p in bp.tolist():
        if p*p>=hi: break
        st=max(p*p,((lo+p-1)//p)*p)
        seg[st-lo::p]=False
    pr=np.nonzero(seg)[0].astype(np.int64)+lo
    lam=np.log(pr.astype(np.float64))
    m=(ppn>=lo)&(ppn<hi)
    n=np.concatenate([pr,ppn[m]]); la=np.concatenate([lam,ppl[m]])
    o=np.argsort(n,kind='stable'); n=n[o]; la=la[o]
    nf=n.astype(np.float64)
    cps=np.cumsum(la); cL=np.cumsum(la/nf)
    psi=psi_b+cps; L=L_b+cL
    # interval list: previous carry interval + intervals starting at each event
    starts=nf; Ls=L; ps=psi
    if carry is not None:
        starts=np.concatenate([[carry[0]],starts]); Ls=np.concatenate([[carry[2]],Ls]); ps=np.concatenate([[carry[1]],ps])
    ends=np.concatenate([starts[1:],[float(hi) if hi<=XMAX else float(XMAX)]])
    # endpoint values (left end and right-limit), interior max
    def Gz(X,Lj,pj): return X*Lj-pj-X*np.log(X)+(1+g)*X
    Ga=Gz(starts,Ls,ps); Gb=Gz(ends,Ls,ps)
    xs=np.exp(Ls+g); inn=(xs>starts)&(xs<ends)
    Xc=np.concatenate([starts,ends,xs[inn]]); Gc=np.concatenate([Ga,Gb,Gz(xs[inn],Ls[inn],ps[inn])])
    keep=Xc>=1.0
    upd(Xc[keep],Gc[keep])
    s=np.sign(Ga-L2P); s=s[s!=0]
    if s.size:
        if lastsign is not None and s[0]!=lastsign: signchg+=1
        signchg+=int((s[1:]!=s[:-1]).sum()); lastsign=s[-1]
    for X in SAMP:
        if starts[0]<=X<=ends[-1]:
            j=int(np.searchsorted(starts,X,side='right'))-1
            samples[X]=float(Gz(np.array([float(X)]),Ls[j:j+1],ps[j:j+1])[0])
    carry=(float(n[-1]),float(psi[-1]),float(L[-1])) if n.size else carry
    psi_b=math.fsum([psi_b,float(math.fsum(la.tolist()))]) if n.size else psi_b
    L_b=math.fsum([L_b,float(math.fsum((la/nf).tolist()))]) if n.size else L_b
    lo=hi
res=dict(XMAX=XMAX,C=C,secs=time.time()-t0,signchanges=signchg,bins=bins,samples=samples,psi=psi_b,L=L_b)
json.dump(res,open(OUT,'w'))
print('XMAX',XMAX,'secs %.1f'%res['secs'],'psi(XMAX)-XMAX=%.3f'%(psi_b-XMAX),'sign changes of Gz-log2pi:',signchg)
gmax=max(v[0] for k,v in bins.items() if k>=4); gmin=min(v[2] for k,v in bins.items() if k>=4)
print('X>=16: max (Gz-log2pi)/sqrtX = %.6f, min = %.6f, bound C = %.7f'%(gmax,gmin,C))
for k in sorted(bins):
    v=bins[k]
    if k%3==0 or k>=27: print('bin 2^%d: max %.5f @%.4g  min %.5f @%.4g  max|Gz|/sqrtX %.5f  frac>0 %.3f'%(k,v[0],v[1],v[2],v[3],v[4],v[5]/max(v[6],1)))
for X,G in sorted(samples.items()): print('Gz(%g) = %.6f   (Gz-log2pi)/sqrtX = %.6f'%(X,G,(G-L2P)/math.sqrt(X)))
