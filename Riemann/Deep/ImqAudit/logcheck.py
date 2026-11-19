import re,glob,os
L=[]
for f in sorted(glob.glob('Riemann/Deep/ImqClass2/sieve_logs/p_*.log')):
    s=open(f).read()
    hd=re.search(r'isieve2 N=(\d+) LO=(\d+) HI=(\d+) pairs (\w+) PT=(\d+)',s)
    dn=re.findall(r'DONE N=(\d+) \[(\d+),(\d+)\]: lanes (\d+) = certF (\d+) \+ cert1 (\d+) \+ t3 (\d+) \(t3pass (\d+), nonfund (\d+), t4 (\d+)\); tested (\d+), tier3 rechecked (\d+), bad (\d+)',s)
    assert hd and len(dn)==1, f
    mm='LANE COUNT MISMATCH' in s or 'FAIL' in s
    N,lo,hi,pairs,PT=hd.groups(); d=[int(x) for x in dn[0]]
    L.append((d[1],d[2],os.path.basename(f),int(N),pairs,int(PT),d,mm))
L.sort()
def cnt(lo,hi):  # number of d in [lo,hi] with d%16 in S
    S=(3,7,11,15,4,8); c=0
    for r in S:
        # count x in [lo,hi], x%16==r
        f=lambda n: (n - r)//16 + 1 if n>=r else 0  # count of x in [0,n] with x%16==r
        c+=f(hi)-f(lo-1)
    return c
tot=0; ok=True; t4=0; prev=None
for lo,hi,name,N,pairs,PT,d,mm in L:
    e=[]
    if N!=10000 or d[0]!=10000: e.append('N')
    if pairs!='on': e.append('pairs')
    if PT!=250000: e.append('PT')
    if d[1]!=lo or d[2]!=hi: e.append('hdr')
    if d[3]!=cnt(lo,hi): e.append('lanes!=count')
    if d[3]!=d[4]+d[5]+d[6]: e.append('sum')
    if d[6]!=d[7]+d[8]+d[9]: e.append('t3sum')
    if d[12]!=0: e.append('bad')
    if mm: e.append('mismatch/FAIL in log')
    if prev is not None and lo!=prev+1: e.append('GAP/OVERLAP prev_hi=%d'%prev)
    prev=hi; tot+=d[3]; t4+=d[9]
    if e: ok=False; print(name,lo,hi,e)
print('chunks',len(L),'first LO',L[0][0],'last HI',L[-1][1],'total lanes',tot,'expected',cnt(28000000000,762700000000),'t4',t4,'ALLOK' if ok else 'PROBLEMS')
