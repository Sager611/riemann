
import math
def primes(N):
    s=bytearray([1])*(N+1); s[0]=s[1]=0
    for i in range(2,int(N**.5)+1):
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
    return [i for i in range(N+1) if s[i]]
P=primes(2*10**6)
for name,K in [('shift eps=.25',math.exp(0.25)),('iid U[-.5,.5]',math.sinh(0.5)/0.5)]:
    lg=0.0
    for p in P: lg+=math.log(1-K/p)-K*math.log(1-1/p)
    G=math.exp(lg)
    print(name,'K=%.6f G=%.5f 1/Gamma(-K)=%.5f'%(K,G,1/math.gamma(-K)))
    for x in [1e6,1e7,1e8]:
        L=math.log(x); print('   x=%.0e lead=%.1f'%(x,G/math.gamma(-K)*x*L**(-1-K)))
# equal weights top at 1e8
Q=60792694; n=5761455; print('equal-weights top |D|-2n =',Q-2*n, (Q-2*n)/Q)

