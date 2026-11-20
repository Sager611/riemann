
import numpy as np, math, sys, time
t0=time.time()
x=6e12; L=math.log(x/(4*math.pi)); pr=L/4; us=L/2
th=0.8733458561067185768; y0=2*th-1
tB=16/L**2*math.log(math.cosh(L*y0/4))
print(f"L={L:.6f} pi*rho={pr:.6f} u*={us:.6f} y0={y0:.10f} tB={tB:.7f}")
NEX=10**7
s=np.ones(NEX+1,dtype=bool); s[:2]=False
for i in range(2,int(NEX**0.5)+1):
    if s[i]: s[i*i::i]=False
P=np.nonzero(s)[0].astype(np.int64); del s
ns=[];lam=[]
for k in range(1,25):
    Q=P[P<=int(NEX**(1.0/k))+1]
    if len(Q)==0: break
    pk=Q**k; m=(pk<=NEX)&(pk>1)
    if not m.any(): break
    ns.append(pk[m]); lam.append(np.log(Q[m].astype(float)))
n=np.concatenate(ns); lm=np.concatenate(lam); del ns,lam
u=np.log(n.astype(float)); mu=lm/np.sqrt(n.astype(float))
o=np.argsort(u); u=u[o]; mu=mu[o]; del n,lm,o
print("prime powers", len(u), "sum mu", mu.sum(), "time", round(time.time()-t0,1))
# check: sum_{n<=N} Lambda(n) n^{-1/2} ~ 2 sqrt(N)
print("check 2sqrt(N)=", 2*math.sqrt(NEX))
np.savez("/tmp/ceil3_lambda.npz", u=u, mu=mu)
