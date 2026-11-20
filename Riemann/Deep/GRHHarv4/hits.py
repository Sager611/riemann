import sys, re, glob, os
# usage: python3 hits.py <id> [maxlines] [width]: GRH/ERH/zero-free hits with file:line, truncated
i=sys.argv[1]; mx=int(sys.argv[2]) if len(sys.argv)>2 else 12; w=int(sys.argv[3]) if len(sys.argv)>3 else 150
pat=re.compile(r'GRH|ERH\b|[Gg]enerali[sz]ed Riemann|[Ee]xtended Riemann|zero-free|zero free|quasi-?R|Linnik|Chebotarev')
out=[]
for f in sorted(glob.glob(i+'/**/*.tex',recursive=True)):
    for k,l in enumerate(open(f,errors='replace'),1):
        if l.lstrip().startswith('%'): continue
        if pat.search(l):
            m=pat.search(l); s=max(0,m.start()-w//3)
            out.append(f"{os.path.basename(f)}:{k}: "+l.strip()[s:s+w])
print(i, len(out), 'hits')
for o in out[:mx]: print('  '+o)
