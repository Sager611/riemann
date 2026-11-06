"""R408 check (orchestrator 21, 6 Nov 2026): imaginary quadratic fields with Polya index 2, unconditionally.
Exact class numbers h(-d) for every d <= 435,435 by counting reduced forms (a, b, c), |b| <= a <= c, b >= 0 if |b| = a or a = c.
Index 2 <=> h = 2 * #Cl[2] = 2^t, t = number of primes dividing d (genus theory; Akbary-Maarefparvar 2508.11125 l. 189-192, 381).
Every such class group is (Z/2)^(t-2) x Z/4, of exponent 4, so it lies in EKN's exponent-4 list (largest |D| = 435,435),
complete without GRH by Exp48 (Riemann/Deep/Exp48/REPORT.md). Output (2.4 s): fundamental 132365 index2 161 max 40755
one-class-per-genus 65 max 5460 (sanity: the classical 65); the 161 equal the paper's Table 'imaginiary quadratic t=2' exactly.
Run: python3 polya_index2.py
"""
import numpy as np, math, sys, json
N=435435
cnt=np.zeros(N+1,dtype=np.int32)
a=1
while 3*a*a<=N:
    for b in range(-(a-1),a+1):
        c0=a if b>=0 else a+1
        cmax=(N+b*b)//(4*a)
        if cmax<c0: continue
        s=4*a*c0-b*b; e=4*a*cmax-b*b
        cnt[s:e+1:4*a]+=1
    a+=1
# smallest prime factor sieve
spf=np.zeros(N+1,dtype=np.int32)
for p in range(2,N+1):
    if spf[p]==0:
        spf[p::p][spf[p::p]==0]=p
def fac(n):
    f=[]
    while n>1:
        p=int(spf[n]); k=0
        while n%p==0: n//=p; k+=1
        f.append((p,k))
    return f
def fund(d):
    if d%4==3:
        return all(k==1 for p,k in fac(d))
    if d%4==0:
        m=d//4
        if m%4 in (1,2): return all(k==1 for p,k in fac(m))
    return False
idx2=[]; one=[]; nf=0
for d in range(3,N+1):
    if not fund(d): continue
    nf+=1
    t=len(fac(d)); h=int(cnt[d])
    if h==2**t: idx2.append(d)
    if h==2**(t-1): one.append(d)
print('fundamental', nf, 'index2', len(idx2), 'max', max(idx2), 'one-class-per-genus', len(one), 'max', max(one))
print('h(-23)',cnt[23],'h(-163)',cnt[163],'h(-4)',cnt[4],'h(-3)',cnt[3],'h(-435435)',cnt[435435])
