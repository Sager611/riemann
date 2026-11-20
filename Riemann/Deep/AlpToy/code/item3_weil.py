#!/usr/bin/env python3
"""alptoy item 3 (note items 4, 6): Weil quadratic form on indicators, normalised by the explicit formula.
W(f,g) := sum_rho fhat(rho) * ghat(1-rho), fhat(s) = int_0^inf f(x) x^(s-1) dx  (symmetric sum over zeros).
1_(0,a]^(s) = a^s/s  =>  W(1_(0,a],1_(0,b]) = b*H(a/b) (a>=b), H(y) := sum_rho y^rho/(rho(1-rho)),  H(y) = y H(1/y).
Explicit formula (proved from psi(t) = t - sum t^rho/rho - log 2pi - 1/2 log(1-t^-2)):
   Gz(y) = log 2pi - r(y) + H(y),  r(y) = sum_k y^(-2k)/(2k(2k+1)),  Gz from Formalism/StatementOC.lean:11-13.
Prime side: H_P(y) = Gz(y) - log 2pi + r(y), exact from prime powers <= y (y <= 1e7).
Zero side: H_Z(y) = 2 sqrt(y) sum_{j<=1e5} cos(g_j log y)/(1/4+g_j^2) (RH-form; valid since these zeros are on the line)."""
import numpy as np, math, sys, json, mpmath as mp
Z=np.loadtxt(sys.argv[1]); OUT=sys.argv[2]; NP=int(1e7)
g=0.57721566490153286061; L2P=math.log(2*math.pi); H1=2+g-math.log(4*math.pi)
T0=float(Z[-1]); nZ=len(Z)
# tail bound: sum_{gamma>T0} 1/gamma^2 = -N(T0)/T0^2 + 2 int_T0^inf N(t)/t^3 dt, N(T0)=1e5 (Odlyzko table complete);
# N(t) <= t/(2pi) log(t/(2pi e)) + 7/8 + 0.112 log t + 0.278 log log t + 2.510 + 0.2/t  (Trudgian 2014)
mp.mp.dps=30
Nup=lambda t: t/(2*mp.pi)*mp.log(t/(2*mp.pi*mp.e))+mp.mpf(7)/8+0.112*mp.log(t)+0.278*mp.log(mp.log(t))+2.510+0.2/t
tail=float(-nZ/mp.mpf(T0)**2+2*mp.quad(lambda t: Nup(t)/t**3,[T0,10*T0,mp.inf]))
tail2=2*tail   # bound for sum over +-gamma of 1/|rho|^2 beyond T0
print('zeros',nZ,'T0=%.3f'%T0,'tail bound sum_{|g|>T0} 1/|rho|^2 <= %.3e'%tail2)
w=1.0/(0.25+Z*Z)
print('H(1): closed form 2+g-log4pi = %.10f ; zeros 2*sum 1/(1/4+g^2) = %.10f ; diff %.3e (<= tail %.3e)'%(H1,2*w.sum(),H1-2*w.sum(),tail2))
def HZ(y):
    y=np.atleast_1d(np.asarray(y,dtype=float)); out=np.empty(len(y))
    for i in range(0,len(y),64):
        ly=np.log(y[i:i+64])[:,None]
        out[i:i+64]=2*np.sqrt(y[i:i+64])*(np.cos(ly*Z[None,:])*w[None,:]).sum(1)
    return out
# prime side tables
bs=np.ones(NP+1,dtype=bool); bs[:2]=False
for i in range(2,int(NP**.5)+1):
    if bs[i]: bs[i*i::i]=False
pr=np.nonzero(bs)[0]
ev=[pr]; lam=[np.log(pr)]
for j in range(2,25):
    q=pr[pr**j<=NP] if True else None
    q=pr[pr.astype(float)**j<=NP]
    if len(q)==0: break
    ev.append(q**j); lam.append(np.log(q))
n=np.concatenate(ev); la=np.concatenate(lam); o=np.argsort(n); n=n[o]; la=la[o]
cps=np.cumsum(la); cL=np.cumsum(la/n)
def r_(y):
    x=1.0/y
    return np.where(y>1.0+1e-12, 1-0.5*np.log1p(-x*x)-np.arctanh(np.minimum(x,1-1e-17))/x, 1-math.log(2))
def Gz(y):
    y=np.atleast_1d(np.asarray(y,dtype=float)); j=np.searchsorted(n,y,side='right')-1
    Ls=np.where(j>=0,cL[np.maximum(j,0)],0.0); ps=np.where(j>=0,cps[np.maximum(j,0)],0.0)
    return y*Ls-ps-y*np.log(y)+(1+g)*y
def HP(y): y=np.atleast_1d(np.asarray(y,dtype=float)); return Gz(y)-L2P+r_(y)
# 1) explicit-formula identity H_P = H_Z on a grid in [1,1e7]
ys=np.unique(np.concatenate([np.exp(np.linspace(0,math.log(NP),1500)),np.array([10.0**k for k in range(8)])]))
hp=HP(ys); hz=HZ(ys); d=np.abs(hp-hz)/np.sqrt(ys)
print('identity H_P=H_Z on %d points y in [1,1e7]: max |H_P-H_Z|/sqrt(y) = %.3e  (tail bound %.3e)'%(len(ys),d.max(),tail2))
res=dict(T0=T0,tail=tail2,H1=H1,H1_zeros=float(2*w.sum()),identity_max=float(d.max()))
# 2) W values and item 6: Q(c) = W(1_(0,1] + c 1_(0,1/X]) = H1 (1 + c^2/X) + 2 c H(X)/X, c = kappa X^eta
Xs=np.exp(np.linspace(math.log(2),math.log(NP),4000)); Hx=HP(Xs)
print('W(1_(0,1]) = H(1) = %.10f > 0 ; W(1_(0,1/X]) = H(1)/X ; cross = H(X)/X ; max |H(X)|/(H(1) sqrt X) on grid = %.4f'%(H1,np.max(np.abs(Hx)/np.sqrt(Xs))/H1))
item6={}
for eta in [0.1,0.25,0.4,0.5,0.6,0.75]:
    for kap in [0.1,1.0,10.0]:
        c=kap*Xs**eta; Q=H1*(1+c*c/Xs)+2*c*Hx/Xs
        item6[f'{eta},{kap}']=float((Q/H1).min())
print('item6 min_X Q/H(1) over X in [2,1e7]:',{k:round(v,4) for k,v in item6.items()})
# literal brief version W(1_(0,1]) + c X^eta W(1_(0,1/X]) = H1(1 + c X^(eta-1)) > 0 trivially
# worst c for each X (c* = -H(X)/H1): Q_min = H1 - H(X)^2/(X H1) >= 0  <=>  |H(X)| <= H1 sqrt X   (2x2 Gram determinant)
det2=(H1*H1/Xs-Hx*Hx/Xs**2)*Xs/H1**2
print('2x2 Gram det (normalised: 1 - H(X)^2/(H1^2 X)) min over grid = %.4f'%det2.min())
res['item6']=item6; res['det2_min']=float(det2.min())
# 3) k x k Gram matrices of 1_(0,x_i], x_i = X^(-i/(k-1)); certified via Weyl: lambda_min > k * entry error
def gram(xs):
    # normalised Gram matrix of f_i = 1_(0,x_i]/sqrt(x_i) (congruent to the plain Gram matrix: same signature)
    # N_ij = H(y)/sqrt(y), y = max/min; float error per entry <= dH(y)/sqrt(y), dH(y) = 1e-14*y*log(y+2)+1e-15 (~50 ulp of the terms)
    k=len(xs); G=np.empty((k,k)); E=np.empty((k,k))
    for i in range(k):
        for j in range(k):
            a,b=max(xs[i],xs[j]),min(xs[i],xs[j]); y=a/b
            G[i,j]=HP(y)[0]/math.sqrt(y); E[i,j]=(1e-14*y*math.log(y+2)+1e-15)/math.sqrt(y)
    return G,float(np.sqrt((E*E).sum()))
out=[]; rng=np.random.default_rng(1)
for k in [2,3,4,6,8]:
    for X in [1e2,1e4,1e6,1e7]:
        xs=X**(-np.arange(k)/(k-1)); G,err=gram(xs); ev_=np.linalg.eigvalsh(G)
        out.append(dict(k=k,X=X,lmin=float(ev_[0]),lmax=float(ev_[-1]),err=err,det=float(np.linalg.det(G))))
for trial in range(200):
    k=int(rng.integers(2,9)); xs=np.sort(np.exp(-rng.uniform(0,math.log(NP),k)))[::-1]
    if np.min(np.abs(np.diff(np.log(xs))))<1e-3: continue
    G,err=gram(xs); ev_=np.linalg.eigvalsh(G)
    out.append(dict(k=k,X=float(1/xs[-1]),lmin=float(ev_[0]),lmax=float(ev_[-1]),err=err,det=float(np.linalg.det(G)),rand=True))
neg=[o for o in out if o['lmin']<-o['err']]; unc=[o for o in out if abs(o['lmin'])<=o['err']]
print('Gram PSD tests (normalised, Weyl-certified lambda_min > ||E||_F): %d matrices (k<=8, points in [1e-7,1]); lambda_min < -err: %d ; uncertified: %d ; min lambda_min = %.3e ; min lambda_min/err = %.3e'%(len(out),len(neg),len(unc),min(o['lmin'] for o in out),min(o['lmin']/o['err'] for o in out)))
for o in out[:20:3]: print('  k=%d X=%.0e lmin=%.4e lmax=%.4e err=%.1e'%(o['k'],o['X'],o['lmin'],o['lmax'],o['err']))
res['gram']=out
# 4) cross-check item 4's sieve values at 1e8, 1e9 against the zero side (|diff| <= tail*sqrt X + ...)
for X,Gs in [(1e8,-109.729763),(1e9,-403.472697)]:
    hz=HZ([X])[0]; print('X=%.0e: sieve Gz=%.6f ; zeros log2pi - r + H_Z = %.6f ; diff %.4f (tail bound %.4f)'%(X,Gs,L2P-r_(np.array([X]))[0]+hz,Gs-(L2P-r_(np.array([X]))[0]+hz),tail2*math.sqrt(X)))
json.dump(res,open(OUT,'w'),default=float)
