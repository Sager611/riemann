
from math import sqrt, prod
from itertools import combinations
def primes(x): return [p for p in range(2,x+1) if all(p%d for d in range(2,int(p**.5)+1))]
for x in [10,30,60,100]:
    P=primes(x); n=len(P)
    # enumerate down-set
    D=[]
    def dfs(i,m,S):
        D.append(tuple(S))
        for j in range(i,n):
            if m*P[j]>x: break
            dfs(j+1,m*P[j],S+[j])
    dfs(0,1,[])
    M=sum((-1)**len(S) for S in D)
    q=[1/(p+1) for p in P]
    tf=0.0
    for S in D:
        s=set(S); w=1.0
        for i in range(n):
            if i in s: w*=q[i]*sqrt(P[i])
            else: w*=(1-q[i])*(-1/sqrt(P[i]))
        tf+=w
    amp=prod(-(p+1)/sqrt(p) for p in P)
    Ef=sum(prod(q[i] if i in set(S) else 1-q[i] for i in range(n)) for S in D)
    print(x,n,'M=',M,'amp*tf=',round(amp*tf,10),'E_pi f=',round(Ef,4),'|tf|=',abs(tf))

