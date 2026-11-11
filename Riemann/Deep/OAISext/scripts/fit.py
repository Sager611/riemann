import os,json,datetime,time
ROOT='/Users/adriansager/playground/riemann'
d=json.load(open(ROOT+'/.scratch/oaisext/closure.json'))
lib=os.path.join(ROOT,'.lake/build/lib/lean')
def lines(m):
    try: return open(os.path.join(ROOT,*m.split('.'))+'.lean',encoding='utf-8').read().count('\n')
    except: return 0
b=[]
for dp,dn,fn in os.walk(os.path.join(lib,'Riemann/Deep/QRH')):
    for f in fn:
        if f.endswith('.olean'):
            p=os.path.join(dp,f); m=os.path.relpath(p,lib)[:-6].replace('/','.')
            b.append((os.path.getmtime(p),m))
b.sort()
def ts(s): return datetime.datetime.strptime(s,'%Y-%m-%dT%H:%M:%SZ').replace(tzinfo=datetime.timezone.utc).timestamp()
B0=ts('2026-11-11T22:58:50Z')
wb=[x for x in b if x[0]>=B0]
rows=[]
prev=None
for t,m in wb:
    if prev is not None: rows.append((t-prev,lines(m),m.split('.')[-1]))
    prev=t
import statistics
xs=[r[1] for r in rows]; ys=[r[0] for r in rows]
n=len(rows); mx=sum(xs)/n; my=sum(ys)/n
bb=sum((x-mx)*(y-my) for x,y in zip(xs,ys))/sum((x-mx)**2 for x in xs); aa=my-bb*mx
print(f'n={n} fit duration = {aa:.1f} s + {bb*1000:.1f} ms/line; median {statistics.median(ys):.1f} s; mean {my:.1f} s')
for r in sorted(rows)[:5]+sorted(rows)[-5:]: print(f'  {r[0]:6.1f}s {r[1]:5d} lines {r[2]}')
# same for window A
A0,A1=ts('2026-11-11T22:49:15Z'),ts('2026-11-11T22:58:33Z')
wa=[x for x in b if A0<=x[0]<A1]; ra=[]; prev=None
for t,m in wa:
    if prev is not None: ra.append((t-prev,lines(m)))
    prev=t
xa=[r[1] for r in ra]; ya=[r[0] for r in ra]; na=len(ra); mxa=sum(xa)/na; mya=sum(ya)/na
ba=sum((x-mxa)*(y-mya) for x,y in zip(xa,ya))/sum((x-mxa)**2 for x in xa); a_=mya-ba*mxa
print(f'A n={na} fit = {a_:.1f} s + {ba*1000:.1f} ms/line; median {statistics.median(ya):.1f} s')
up=set(d['up']); down=set(d['down'])
built=set(m for _,m in b)
rem_up=[m for m in up if m not in built]
fork_lines=[lines(m) for m in down]
ru=[lines(m) for m in rem_up]
for name,(A,B) in {'bg fit':(aa,bb),'utility fit':(a_,ba)}.items():
    tu=sum(A+B*l for l in ru); tf=sum(A+B*l for l in fork_lines)
    print(f'{name}: remaining upstream {len(ru)} modules -> {tu/3600:.1f} h; fork 241 -> {tf/3600:.1f} h; total {(tu+tf)/3600:.1f} h')
