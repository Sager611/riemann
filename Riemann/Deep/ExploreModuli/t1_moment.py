# T1: family-averaged trace-method moment for primitive characters mod prime q.
# Explicit formula, h(r) = R (sin(rR/2)/(rR/2))^2, g(u) = (1-|u|/R)_+ (h >= 0 on R and on iR).
# E_chi sum_rho h(gamma_rho) = arch(q,R) - (2/(q-2)) [ (q-1) S1 - S0 ],
#   S1 = sum_{n = 1 mod q} Lambda(n) n^{-1/2} g(log n),  S0 = sum_{q not | n} Lambda(n) n^{-1/2} g(log n).
import numpy as np, math, sys, time, json
import mpmath as mp
t0=time.time()
X=int(float(sys.argv[1])) if len(sys.argv)>1 else int(1.1e9)
Q=[101,1009,10007,100003]
AG={101:[0.5,1,1.25,1.5,1.75,2,2.5,3,3.5,4,4.5],1009:[0.5,1,1.25,1.5,1.75,2,2.5,3],10007:[0.5,1,1.25,1.5,1.75,2,2.25],100003:[0.5,1,1.25,1.5,1.75,1.8]}
th=[(q,A,q**A) for q in Q for A in AG[q] if q**A<=X]
xs=np.array(sorted(set(t[2] for t in th)))
A0=np.zeros(len(xs));A1=np.zeros(len(xs))
B0={q:np.zeros(len(xs)) for q in Q};B1={q:np.zeros(len(xs)) for q in Q}
r=int(math.isqrt(X))+2
sb=np.ones(r+1,dtype=bool);sb[:2]=False
for i in range(2,int(math.isqrt(r))+1):
    if sb[i]: sb[i*i::i]=False
base=np.nonzero(sb)[0]
seg=1<<24
for lo in range(0,X+1,seg):
    hi=min(lo+seg,X+1)
    s=np.ones(hi-lo,dtype=bool)
    if lo==0: s[:2]=False
    for p in base:
        p=int(p)
        if p*p>=hi: break
        st=max(p*p,((lo+p-1)//p)*p)
        s[st-lo::p]=False
    P=np.nonzero(s)[0].astype(np.int64)+lo
    if len(P)==0: continue
    lp=np.log(P.astype(np.float64)); w=lp/np.sqrt(P.astype(np.float64))
    cw=np.cumsum(w); cwl=np.cumsum(w*lp)
    idx=np.searchsorted(P,xs,side='right'); tk=idx>0
    A0[tk]+=cw[idx[tk]-1]; A1[tk]+=cwl[idx[tk]-1]
    for q in Q:
        m=(P%q)==1
        if not m.any(): continue
        Pq=P[m]; cq=np.cumsum(w[m]); cql=np.cumsum(w[m]*lp[m])
        iq=np.searchsorted(Pq,xs,side='right'); t2=iq>0
        B0[q][t2]+=cq[iq[t2]-1]; B1[q][t2]+=cql[iq[t2]-1]
# prime powers p^k, k>=2
qpow={q:np.zeros(len(xs)) for q in Q}; qpowl={q:np.zeros(len(xs)) for q in Q}
for p in base:
    p=int(p); lpp=math.log(p); n=p*p; k=2
    while n<=X:
        wv=lpp/math.sqrt(n); ln=k*lpp
        msk=xs>=n
        A0[msk]+=wv; A1[msk]+=wv*ln
        for q in Q:
            if n%q==1: B0[q][msk]+=wv; B1[q][msk]+=wv*ln
        n*=p; k+=1
# powers of q itself (excluded from S0)
for q in Q:
    n=q;k=1
    while n<=X:
        wv=math.log(q)/math.sqrt(n); ln=k*math.log(q); msk=xs>=n
        qpow[q][msk]+=wv; qpowl[q][msk]+=wv*ln; n*=q;k+=1
tsieve=time.time()-t0
# arch term: log(q/pi) + E_parity (1/pi) int_0^inf 2 (sin v/v)^2 Re psi((1/2+a+2iv/R)/2) dv
def archpsi(R,a,K=1500):
    tot=mp.mpf(0)
    for k in range(K):
        Ik=mp.si(2*(k+1)*mp.pi)-mp.si(2*k*mp.pi)
        v=(k+0.5)*mp.pi
        tot+=Ik*mp.re(mp.digamma((mp.mpf(1)/2+a+2j*v/R)/2))
    # tail: (sin v/v)^2 ~ 1/(2v^2) on average, Re psi(z) ~ log|z|
    V=K*mp.pi
    tail=mp.quad(lambda v: mp.re(mp.log((mp.mpf(1)/2+a+2j*v/R)/2))/(2*v*v),[V,mp.inf])
    return float(2*(tot+tail)/mp.pi)
rows=[]
for (q,A,x) in th:
    j=int(np.searchsorted(xs,x)); R=A*math.log(q)
    S0=(A0[j]-A1[j]/R)-(qpow[q][j]-qpowl[q][j]/R)
    S1=B0[q][j]-B1[q][j]/R
    D=(q-1)*S1-S0
    ne=(q-3)/2; no=(q-1)/2
    arch=math.log(q/math.pi)+(ne*archpsi(R,0)+no*archpsi(R,1))/(ne+no)
    Eprime=2*D/(q-2)
    Z=arch-Eprime
    Zup=arch+2*S0/(q-2)   # unconditional: drops S1 >= 0
    M=(4/R)*(math.exp(R/2)-1)-2
    rows.append(dict(q=q,A=A,R=R,S0=S0,S1=S1,D=D,arch=arch,Z=Z,Zup=Zup,M=M))
    print(f"q={q:6d} A={A:4.2f} R={R:6.2f} S0={S0:12.2f} M={M:12.2f} (q-1)S1={(q-1)*S1:12.2f} D={D:10.2f} D/sqrt(q)={D/math.sqrt(q):7.3f} arch={arch:6.3f} Z={Z:7.3f} Zup={Zup:10.3f} Zup/arch={Zup/arch:8.3f}",flush=True)
json.dump(rows,open('t1_rows.json','w'),indent=1)
print('sieve s',round(tsieve,1),'total s',round(time.time()-t0,1))
