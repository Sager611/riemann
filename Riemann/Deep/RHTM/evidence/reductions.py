"""Independent re-check of the 122->121->120 reductions.
1) Homomorphism raw-named -> 120 via the names map, listing every (state,read) pair where it fails.
2) Own 5-cell window over-approximation on the RAW named table (not the canonical one).
"""
import json
from collections import deque
def read_named(p):
    t={}
    for L in open(p):
        if L.startswith('#') or not L.strip(): continue
        q,s,w,d,n=L.split(); t[q,int(s)]=(int(w),d,n)
    return t
def read_canon(p):
    t={}
    for L in open(p):
        if L.startswith('#') or not L.strip(): continue
        q,s,w,d,n=L.split(); t[int(q),int(s)]=(int(w),d,n if n=='H' else int(n))
    return t
raw=read_named('../repro/machines/RH_121_raw_named.tm')
t120=read_canon('../repro/machines/RH_120_paper.tm')
pi=json.load(open('../repro/machines/RH_120.names.json'))
states={q for q,s in raw}
print('raw states',len(states),'raw transitions',len(raw))
fails=[]
for q in sorted(states):
    for s in (0,1):
        if (q,s) not in raw:
            fails.append((q,s,'missing in raw ->120 gives',t120[pi[q],s])); continue
        w,d,n=raw[q,s]
        tw,td,tn=t120[pi[q],s]
        img='H' if n=='halt' else pi[n]
        if (w,d,img)!=(tw,td,tn): fails.append((q,s,(w,d,n,img),(tw,td,tn)))
print('homomorphism failures:',fails)
# window invariant on raw table, radius 2
r=2; mask=(1<<(2*r+1))-1
start='0a.boot1.A'
seen=set([(start,0)]); dq=deque([(start,0)])
while dq:
    q,wd=dq.popleft()
    s=(wd>>r)&1
    if (q,s) not in raw: continue   # missing transition: machine would halt (completion) - record
    w,d,n=raw[q,s]
    if n=='halt': continue
    ch=(wd&~(1<<r))|(w<<r)
    for b in (0,1):
        nw=(((ch<<1)&mask)|b) if d=='R' else ((ch>>1)|(b<<(2*r)))
        if (n,nw) not in seen: seen.add((n,nw)); dq.append((n,nw))
reads={(q,(wd>>r)&1) for q,wd in seen}
print('admitted pairs',len(seen))
for q,s in [p[:2] for p in fails]:
    print('pair',(q,s),'reachable under window invariant:',(q,s) in reads)
print('prep_1 read 1 reachable under window invariant:',('1b.reg.prep_1',1) in reads)
unreach=sorted((q,s) for q in states for s in (0,1) if (q,s) not in reads)
print('pairs excluded by invariant:',unreach)
