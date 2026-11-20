#!/usr/bin/env python3
"""alptoy item 2: compact table for t (formal exponential) and sigma, both exponents gam(p)=p/(p^2-1) (Lean InertBd, Fnum)
and rho(p)=1/(p^2-1) (Lean PhiCBd), Formalism/IncomingBounds/IncomingBound.lean:8-10."""
import json, math, sys
rows=json.load(open(sys.argv[1]))+json.load(open(sys.argv[2]))
def lean_form(a,p,w):
    s=0.0; i=1
    while a/p**i>w: s+=a/p**i-w; i+=2
    return s
print('%-10s %-5s %2s %5s %8s %8s %8s | %7s %7s %7s | %7s %7s | %5s %5s'%('K','f','p','n','slope','gam','rho','defG','excG','C_ln','defR','excR','w*','D^p'))
summ=[]
for r in rows:
    if r['kind']!='inert': continue
    p=r['p']; V={int(a):v for a,v in r['V'].items() if v is not None}
    pts=sorted(V.items()); g=p/(p*p-1); rh=1/(p*p-1)
    n=len(pts); sx=sum(a for a,_ in pts); sy=sum(v for _,v in pts); sxx=sum(a*a for a,_ in pts); sxy=sum(a*v for a,v in pts)
    slope=(n*sxy-sx*sy)/(n*sxx-sx*sx)
    defG=max(g*a-v for a,v in pts); excG=max(v-g*a for a,v in pts)
    defR=max(rh*a-v for a,v in pts); excR=max(v-rh*a for a,v in pts)
    C=max((g*a-v)/math.log(a) for a,v in pts if a>=2)
    lo,hi=0.0,60.0
    for _ in range(50):
        mid=(lo+hi)/2
        if all(v>=lean_form(a,p,mid)-1e-12 for a,v in pts): hi=mid
        else: lo=mid
    s=dict(K=r['K'],f=r['f'],p=p,n=n,slope=slope,gam=g,rho=rh,defG=defG,excG=excG,C_ln=C,defR=defR,excR=excR,wstar=hi,Dp=r['Dp_zero'],amax=max(V))
    summ.append(s)
    print('%-10s %-5s %2d %5d %8.5f %8.5f %8.5f | %7.3f %7.3f %7.3f | %7.3f %7.3f | %5.2f %5s'%(r['K'],r['f'],p,n,slope,g,rh,defG,excG,C,defR,excR,hi,r['Dp_zero']))
T=[s for s in summ if s['f']=='t']; S=[s for s in summ if s['f']=='sigma']
print('t (formal exponential, %d cases): max (gam a - v) = %.3f ; max (v - gam a) = %.3f ; max C_ln = %.3f ; max w* = %.3f ; |slope/gam - 1| <= %.2e'%(
   len(T),max(s['defG'] for s in T),max(s['excG'] for s in T),max(s['C_ln'] for s in T),max(s['wstar'] for s in T),max(abs(s['slope']/s['gam']-1) for s in T)))
print('sigma (%d cases): max (rho a - v) = %.3f ; max (v - rho a) = %.3f ; max (gam a - v) = %.1f ; |slope/rho - 1| <= %.3f'%(
   len(S),max(s['defR'] for s in S),max(s['excR'] for s in S),max(s['defG'] for s in S),max(abs(s['slope']/s['rho']-1) for s in S)))
json.dump(summ,open(sys.argv[3],'w'))
