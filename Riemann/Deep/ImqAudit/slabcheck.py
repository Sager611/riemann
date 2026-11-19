import re
s=open('Riemann/Deep/ImqClass2/mjx_slabs.log').read()
R=re.findall(r'(cl\d+mod\d+)\.(\d+): entries (\d+), first (\d+), last (\d+), h<=2000: (\d+), h<=10000: (\d+), min h (\d+) at d (\d+)',s)
S=2**28; seen={}; bad=[]; tot=0
for nm,k,n,f,l,a,b,mh,dm in R:
    k=int(k); f=int(f); l=int(l); tot+=int(n)
    if not (k*S<=f<=l<(k+1)*S): bad.append((nm,k,f,l))
    seen.setdefault(k,set()).add(nm)
ks=sorted(seen)
print('lines',len(R),'slabs',len(ks),'k range',ks[0],ks[-1],'all 4 classes each:',all(len(seen[k])==4 for k in ks),'contiguous:',ks==list(range(ks[0],ks[-1]+1)))
print('files with first/last outside slab:',len(bad),bad[:5])
print('total entries',tot)
# min h over slabs 105..110 (overlap with sieve)
mins=[(int(mh),nm,k) for nm,k,n,f,l,a,b,mh,dm in R if int(k)>=104]
print('min h in slabs>=104:',sorted(mins)[:3])
