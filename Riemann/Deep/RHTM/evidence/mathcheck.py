import math, sys, time
t0=time.time()
def q_of(n):
    q=0
    while 4**(q+1)<=3*n: q+=1
    return q
def thr(n):
    q=q_of(n); return (q+1)<<q
def beta_exact(b,B):
    # beta = min d with B <= M_d, M_d=b+...+b^d
    if B==0: return 0
    e=max(1,(B.bit_length()-1)*1.0/math.log2(b)); d=max(0,int(e)-2)
    Md=(b**(d+1)-b)//(b-1)
    while B>Md:
        d+=1; Md=Md*b+b
    return d
def beta_loop(b,B):
    c=0
    while B>0: B=(B-1)//b; c+=1
    return c
# 1) exact check for 2<=n<=2656 (Prop small range) + spot check beta_exact==beta_loop
N1=2656
A=1; isp=bytearray([1])*(N1+1); isp[0]=isp[1]=0
for i in range(2,int(N1**0.5)+1):
    if isp[i]: isp[i*i::i]=bytearray(len(isp[i*i::i]))
worst=(0.0,0,0,0)
for n in range(2,N1+1):
    A*=n
    if isp[n]: A*=n
    B=A-1
    b=beta_exact(n,B)
    if n<=200: assert b==beta_loop(n,B),n
    D=max(n-1-b,0); T=thr(n)
    assert D<T,(n,D,T)
    r=D/T
    if r>worst[0]: worst=(r,n,D,T)
print('exact range 2..%d: criterion holds; max Delta/thr = %.4f at n=%s (Delta=%s thr=%s); %.1fs'%(N1,worst[0],worst[1],worst[2],worst[3],time.time()-t0),flush=True)
# 2) clock lemma exhaustively to 10^6
for n in range(2,10**6+1):
    z=n-1;q=0
    while z>0: z=(z-1)//4; q+=1
    if q!=q_of(n): print('CLOCK FAIL',n); break
else: print('clock lemma ok to 1e6',flush=True)
# 3) Prop-small analytic chain: n/log n < thr(n) for 3<=n<=2656 (float, margin printed)
mm=min(thr(n)-n/math.log(n) for n in range(3,2657)); print('min thr - n/log n on [3,2656] = %.3f'%mm)
# 4) float upper bound Delta_n <= max(n-1-floor(L),0), L=log A_n/log n, up to N2
N2=int(sys.argv[1]) if len(sys.argv)>1 else 2*10**6
sieve=bytearray([1])*(N2+1); sieve[0]=sieve[1]=0
for i in range(2,int(N2**0.5)+1):
    if sieve[i]: sieve[i*i::i]=bytearray(len(sieve[i*i::i]))
theta=0.0; comp=0.0; worst=(0.0,0,0,0); close=0; maxD=(0,0)
for n in range(2,N2+1):
    if sieve[n]:
        y=math.log(n)-comp; t=theta+y; comp=(t-theta)-y; theta=t
    L=(math.lgamma(n+1)+theta)/math.log(n)
    fl=math.floor(L)
    if L-fl<1e-7 or fl+1-L<1e-7: close+=1
    Dub=max(n-1-fl,0)
    T=thr(n)
    if Dub>=T: print('FLOAT BOUND FAILS at',n,Dub,T); break
    if Dub/T>worst[0]: worst=(Dub/T,n,Dub,T)
    if Dub>maxD[0]: maxD=(Dub,n)
print('float range to %d: Delta upper bound < thr everywhere; max ratio %.4f at n=%d (Dub=%d thr=%d); max Dub=%s; near-integer L cases=%d; %.1fs'%(N2,worst[0],worst[1],worst[2],worst[3],maxD,close,time.time()-t0))
