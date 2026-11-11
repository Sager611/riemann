# fit_m4.py — table and fits for logs/m4_ladder.log, logs/m4_random.log, logs/m4_prin.log, logs/m4_lengths.log
import re, math, sys
import numpy as np
D='logs/'
def kv(line):
    return {k: v for k, v in re.findall(r'(\w+)=([-+0-9.eE]+)', line)}
runs=[kv(l) for l in open(D+'m4_ladder.log') if l.startswith('run ')]
rm={}
try:
    for l in open(D+'m4_random.log'):
        if l.startswith('rm '): d=kv(l); rm[float(d['X'])]=d
except FileNotFoundError: pass
prin={}
try:
    for l in open(D+'m4_prin.log'):
        if l.startswith('prin '): d=kv(l); prin[float(d['U'])]=d
except FileNotFoundError: pass
print('| U | X | rows | M4_6 | M4_6/rows | E4_6(X) | M4_6/(rows E4_6) | M4_3 | M4_3/rows | M4_3/(rows E4_3) | P4_6/M4_6 | P4_3/M4_3 | time s |')
print('|---|---|---|---|---|---|---|---|---|---|---|---|---|')
U=[];m6=[];m3=[];R=[]
for d in runs:
    u=float(d['U']); x=float(d['X']); rows=float(d['rows']); a=float(d['M4_6']); b=float(d['M4_3'])
    U.append(u); m6.append(a); m3.append(b); R.append(rows)
    e=rm[min(rm,key=lambda k:abs(k-x))] if rm and min(abs(k-x) for k in rm)<0.01*x else None; p=prin.get(u)
    e6=float(e['E4_6']) if e else float('nan'); e3=float(e['E4_3']) if e else float('nan')
    p6=float(p['P4_6'])/a if p else float('nan'); p3=float(p['P4_3'])/b if p else float('nan')
    print(f"| {u:.0e} | {x:.0f} | {rows:.0f} | {a:.4e} | {a/rows:.5f} | {e6:.5f} | {a/rows/e6:.4f} | {b:.4e} | {b/rows:.5f} | {b/rows/e3:.4f} | {p6:.3g} | {p3:.3g} | {float(d['time']):.1f} |")
U=np.array(U); m6=np.array(m6); m3=np.array(m3); R=np.array(R)
Xs=np.sqrt(U); E6=np.array([float(rm[min(rm,key=lambda k:abs(k-x))]['E4_6']) for x in Xs]); E3=np.array([float(rm[min(rm,key=lambda k:abs(k-x))]['E4_3']) for x in Xs])
for name,y,E in (('sextic',m6,E6),('cubic',m3,E3)):
    s=U>=1e6; c=np.polyfit(np.log(U[s]),np.log(y[s]/(R[s]*E[s])),1); print(f'{name}: M4/(rows*E4) ~ U^{c[0]:+.4f} over U>=1e6 (random-model ratio; MC error of E4 ~0.5%)')
for name,y in (('sextic',m6),('cubic',m3)):
    loc=[math.log(y[i+1]/y[i])/math.log(U[i+1]/U[i]) for i in range(len(U)-1)]
    print(name,'local exponents:',' '.join(f'{v:.4f}' for v in loc))
    for lo in (1e5,1e6):
        s=U>=lo; A=np.vstack([np.log(U[s]),np.ones(s.sum())]).T; c,res,_,_=np.linalg.lstsq(A,np.log(y[s]),rcond=None)
        rss=float(((A@c-np.log(y[s]))**2).sum())
        # alternative: y = rows*(a+b*log U)
        B=np.vstack([R[s],R[s]*np.log(U[s])]).T; c2,_,_,_=np.linalg.lstsq(B,y[s],rcond=None); rss2=float(((np.log(B@c2)-np.log(y[s]))**2).sum())
        print(f'  fit U>={lo:.0e}: M4 ~ U^{c[0]:.4f} (RSS log {rss:.2e}); M4 = rows*({c2[0]:.4f}+{c2[1]:.5f} log U) (RSS log {rss2:.2e})')
if prin:
    us=sorted(prin); 
    for key,lab in (('P4_6','sixth-power rows'),('P4_3','cube rows')):
        y=np.array([float(prin[u][key]) for u in us]); x=np.log(np.array(us))
        s=np.array(us)>=1e5; c=np.polyfit(x[s],np.log(y[s]),1)
        print(f'excluded {lab}: contribution ~ U^{c[0]:.4f} (paper/heuristic: U^(1+1/{6 if key=="P4_6" else 3}))')
try:
    for l in open(D+'m4_lengths.log'):
        if l.startswith('run '): d=kv(l); print('length run: U=%s m=%s rows=%s M4_6/rows=%.5f M4_3/rows=%.5f'%(d['U'],d['m'],d['rows'],float(d['M4_6'])/float(d['rows']),float(d['M4_3'])/float(d['rows'])))
except FileNotFoundError: pass

# Q(zeta_12)
try:
    zr=[kv(l) for l in open(D+'m4z12.log') if l.startswith('z12 ')]
    zm={float(d['U']):d for d in (kv(l) for l in open(D+'m4z12.log') if l.startswith('z12rm '))}
    print('\nQ(zeta_12):')
    print('| U | X | rows | M4_6 | M4_6/rows | E4_6(X) | M4_6/(rows E4_6) | M4_3 | M4_3/rows | M4_3/(rows E4_3) | time s |')
    print('|---|---|---|---|---|---|---|---|---|---|---|')
    UU=[];y6=[];y3=[];RR=[];e6=[];e3=[]
    for d in zr:
        u=float(d['U']); rows=float(d['rows']); a=float(d['M4_6']); b=float(d['M4_3']); e=zm.get(u)
        E6=float(e['E4_6']) if e else float('nan'); E3=float(e['E4_3']) if e else float('nan')
        UU.append(u); y6.append(a); y3.append(b); RR.append(rows); e6.append(E6); e3.append(E3)
        print(f"| {u:.0e} | {float(d['X']):.0f} | {rows:.0f} | {a:.4e} | {a/rows:.5f} | {E6:.5f} | {a/rows/E6:.4f} | {b:.4e} | {b/rows:.5f} | {b/rows/E3:.4f} | {float(d['time']):.1f} |")
    UU=np.array(UU); y6=np.array(y6); y3=np.array(y3); RR=np.array(RR); e6=np.array(e6); e3=np.array(e3)
    for name,y,E in (('sextic',y6,e6),('cubic',y3,e3)):
        s=UU>=1e5; A=np.vstack([np.log(UU[s]),np.ones(s.sum())]).T; c,_,_,_=np.linalg.lstsq(A,np.log(y[s]),rcond=None); rss=float(((A@c-np.log(y[s]))**2).sum())
        ok=s&np.isfinite(E); msg=''
        if ok.sum()>=2: c2=np.polyfit(np.log(UU[ok]),np.log(y[ok]/(RR[ok]*E[ok])),1); msg=f'; M4/(rows*E4) ~ U^{c2[0]:+.4f}'
        print(f'z12 {name}: fit U>=1e5: M4 ~ U^{c[0]:.4f} (RSS log {rss:.2e}){msg}')
except FileNotFoundError: pass
