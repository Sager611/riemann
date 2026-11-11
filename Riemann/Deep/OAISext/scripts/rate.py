import os,json,re,time,datetime,sys
ROOT='/Users/adriansager/playground/riemann'
d=json.load(open(ROOT+'/.scratch/oaisext/closure.json'))
up=set(d['up']); down=set(d['down'])
def src(m): return os.path.join(ROOT,*m.split('.'))+'.lean'
def lines(m):
    try: return open(src(m),encoding='utf-8').read().count('\n')
    except: return 0
lib=os.path.join(ROOT,'.lake/build/lib/lean')
built={}
for base in ['Riemann/Deep/QRH','Riemann/Deep/OAISext']:
    for dp,dn,fn in os.walk(os.path.join(lib,base)):
        for f in fn:
            if f.endswith('.olean'):
                rel=os.path.relpath(os.path.join(dp,f),lib)[:-6].replace('/','.')
                built[rel]=os.path.getmtime(os.path.join(dp,f))
def ts(s): return datetime.datetime.strptime(s,'%Y-%m-%dT%H:%M:%SZ').replace(tzinfo=datetime.timezone.utc).timestamp()
A0,A1=ts(sys.argv[1]),ts(sys.argv[2]); B0=ts(sys.argv[3]); now=time.time()
wa=[m for m,t in built.items() if A0<=t<A1]; wb=[m for m,t in built.items() if t>=B0]
la=sum(lines(m) for m in wa); lb=sum(lines(m) for m in wb)
da=A1-A0; db=now-B0
fork=set('Riemann.Deep.OAISext.'+m[len('Riemann.Deep.QRH.'):] for m in down)|{'Riemann.Deep.OAISext.BetaEff','Riemann.Deep.OAISext.Wrapper','Riemann.Deep.OAISext.Balanced'}
need=up|fork
rem=[m for m in need if m not in built]
rl=sum(lines(m) for m in rem)
print(f'built total {len(built)}; window A (utility) {len(wa)} modules {la} lines in {da/60:.1f} min: {da/max(1,len(wa)):.1f} s/module, {la/da*60:.0f} lines/min')
print(f'window B (background) {len(wb)} modules {lb} lines in {db/60:.1f} min: {db/max(1,len(wb)):.1f} s/module, {lb/max(db,1)*60:.0f} lines/min')
avgA=la/max(1,len(wa)); avgB=lb/max(1,len(wb))
print(f'avg lines/module: A {avgA:.0f}, B {avgB:.0f}; remaining needed (up + fork) {len(rem)} modules, {rl} lines (avg {rl/max(1,len(rem)):.0f})')
if wb:
    print(f'projection B: by modules {len(rem)*db/len(wb)/3600:.1f} h; by lines {rl/(lb/db)/3600:.1f} h')
if wa:
    print(f'projection A: by modules {len(rem)*da/len(wa)/3600:.1f} h; by lines {rl/(la/da)/3600:.1f} h')
