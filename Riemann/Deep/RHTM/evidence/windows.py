import sys, time
from collections import deque
def read_canon(p):
    t={}
    for L in open(p):
        if L.startswith('#') or not L.strip(): continue
        q,s,w,d,n=L.split(); t[int(q),int(s)]=(int(w),1 if d=='R' else -1,-1 if n=='H' else int(n))
    return t
t=read_canon(sys.argv[1]); nst=len(t)//2
for r in map(int,sys.argv[2:]):
    t0=time.time(); mask=(1<<(2*r+1))-1
    seen=[set() for _ in range(nst)]; seen[0].add(0); dq=deque([(0,0)]); cnt=1
    while dq:
        q,wd=dq.popleft(); s=(wd>>r)&1
        w,d,n=t[q,s]
        if n<0: continue
        ch=(wd&~(1<<r))|(w<<r)
        for b in (0,1):
            nw=(((ch<<1)&mask)|b) if d==1 else ((ch>>1)|(b<<(2*r)))
            if nw not in seen[n]: seen[n].add(nw); dq.append((n,nw)); cnt+=1
    reads={(q,(wd>>r)&1) for q in range(nst) for wd in seen[q]}
    excl=sorted((q,s) for q in range(nst) for s in (0,1) if (q,s) not in reads)
    print('radius',r,'pairs',cnt,'excluded',excl,'%.1fs'%(time.time()-t0),flush=True)
