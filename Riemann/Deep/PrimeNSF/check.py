# Every n > 24 is a prime plus a non-squarefree integer (COROLLARIES; arXiv 2605.02426 Thm main3 without GRH).
# Run: python3 check.py (numpy only; about 11 s, < 200 MB). Output kept in check_output.txt.
import numpy as np, math, time
t0=time.time()
N=10**7
s=np.ones(N+1,dtype=bool); s[:2]=False
for i in range(2,int(N**0.5)+1):
    if s[i]: s[i*i::i]=False
P=np.nonzero(s)[0]
# (2) least primes in reduced classes mod q^2, q prime <= 200
worst=(0,0,0); ok=True
for q in [int(x) for x in P[P<=200]]:
    m=q*q; r=P%m
    u,idx=np.unique(r,return_index=True)
    red=[(a,P[i]) for a,i in zip(u,idx) if a%q!=0]
    if len(red)!=q*(q-1): ok=False; print('MISSING classes for q',q,len(red),q*(q-1))
    mx=max(p for a,p in red)
    if mx>worst[0]: worst=(int(mx),q,0)
print('q<=200 all reduced classes hit:',ok,' max least prime',worst[0],'at q',worst[1])
# (3) q in (200,3000]: 7.8978*log(0.72 q^2 log q^2) < theta(q-1)
th=0.0; pi=0; minmarg=1e9; 
Plist=[int(x) for x in P[P<=3000]]
cum={}
acc=0.0
for p in Plist:
    cum[p]=acc   # theta(p-1)
    acc+=math.log(p)
for q in Plist:
    if q<=200: continue
    m=q*q; lhs=7.8978*math.log(0.72*m*math.log(m)); rhs=cum[q]
    minmarg=min(minmarg, rhs-lhs)
print('q in (200,3000]: min theta(q-1) - 7.8978 log(0.72 q^2 log q^2) =', round(minmarg,3))
# (4) q > 3000: Rosser-Schoenfeld theta(x) > x(1-1/log x) for x >= 41; check at q=3001 and growth
q=3001; lhs=7.8978*math.log(0.72*q*q*math.log(q*q)); rhs=(q-1)*(1-1/math.log(q-1))
print('q=3001: RS lower bound', round(rhs,1), 'vs', round(lhs,1), '(lhs grows like 15.8 log q, rhs like q)')
# (5) every n in [25, N]: prime p < n with n - p divisible by some s in squares
sq=[4,9,25,49,121,169,289,361,529,841,961]
n=np.arange(N+1); good=np.zeros(N+1,dtype=bool)
for S in sq:
    L=np.full(S,np.iinfo(np.int64).max,dtype=np.int64)
    r=P%S
    u,idx=np.unique(r,return_index=True)
    L[u]=P[idx]
    good|= (n > L[n%S])
bad=np.nonzero(~good[25:])[0]+25
print('n in [25,1e7] without p:', len(bad), bad[:10])
print('time', round(time.time()-t0,1))
