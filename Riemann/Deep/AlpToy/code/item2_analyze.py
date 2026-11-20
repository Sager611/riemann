#!/usr/bin/env python3
"""alptoy item 2 analysis: compare v_p(D^a f(0)) with a/p + a/p^3 + ... = a p/(p^2-1) and with the Lean form
sum_{i odd} max(a/p^i - w, 0) (Formalism/SizeWindowBrackets/FirstPrimeNeSecondPrime.lean:75-77 askP; odd i = 'not InPhi (i * frob p)')."""
import json, math, sys
rows=json.load(open(sys.argv[1]))
STRUCT={('Q(i)','t'):lambda a:a%4==1,('Q(i)','sigma'):lambda a:a%4==1,('Q(sqrt-3)','t'):lambda a:a%3==1,
        ('Q(sqrt-3)','sigma'):lambda a:a%6==1,('Q(sqrt-7)','t'):lambda a:True,('Q(sqrt-7)','sigma'):lambda a:a%2==1,
        ('Q(sqrt-23)','t'):lambda a:True,('Q(sqrt-23)','sigma'):lambda a:a%2==1}
def vfact(a,p):
    s=0; q=p
    while q<=a: s+=a//q; q*=p
    return s
def lean_form(a,p,w):
    s=0.0; i=1
    while a/p**i>w: s+=a/p**i-w; i+=2
    return s
out=[]
print('%-10s %-5s %2s %-5s %5s %8s %8s %8s %8s %8s %6s %6s'%('K','f','p','kind','n','slope','p/(p2-1)','maxDef','C_ln','C_logp','w*','D^p0'))
for r in rows:
    p=r['p']; K=r['K']; f=r['f']; V={int(a):v for a,v in r['V'].items()}; st=STRUCT[(K,f)]
    pts=[(a,v) for a,v in V.items() if v is not None]
    unexp=[a for a,v in V.items() if v is None and st(a)]          # zero mod p^prec where nonzero expected
    nonstruct=[a for a,v in pts if not st(a)]                        # nonzero where structural zero expected
    if r['kind']=='split':
        mv=max(v for a,v in pts); out.append(dict(K=K,f=f,p=p,kind='split',maxv=mv,n=len(pts),unexp=len(unexp)))
        print('%-10s %-5s %2d %-5s %5d   max v_p = %d (no growth)'%(K,f,p,'split',len(pts),mv)); continue
    c=p/(p*p-1)
    n=len(pts); sx=sum(a for a,_ in pts); sy=sum(v for _,v in pts); sxx=sum(a*a for a,_ in pts); sxy=sum(a*v for a,v in pts)
    slope=(n*sxy-sx*sy)/(n*sxx-sx*sx)
    defs=[(c*a-v,a) for a,v in pts]
    maxdef=max(defs)
    C_ln=max((c*a-v)/math.log(a) for a,v in pts if a>=2)
    C_lp=max((c*a-v)/(math.log(a)/math.log(p)) for a,v in pts if a>=2)
    minexc=min(v-c*a for a,v in pts)
    # minimal w with v >= lean_form(a,p,w) for all a
    lo,hi=0.0,50.0
    for _ in range(60):
        mid=(lo+hi)/2
        if all(v>=lean_form(a,p,mid)-1e-12 for a,v in pts): hi=mid
        else: lo=mid
    w1=all(v>=lean_form(a,p,1.0)-1e-12 for a,v in pts)
    fails1=[a for a,v in pts if v<lean_form(a,p,1.0)-1e-12]
    full_fail=sum(1 for a,v in pts if v<vfact(a,p))
    worst_full=max(vfact(a,p)-v for a,v in pts)
    out.append(dict(K=K,f=f,p=p,kind='inert',n=n,slope=slope,c=c,maxdef=maxdef[0],argmaxdef=maxdef[1],C_ln=C_ln,C_logp=C_lp,
                    min_excess=minexc,wstar=hi,w1_holds=w1,w1_fails=fails1[:10],n_w1_fails=len(fails1),full_fact_fail=full_fail,worst_full=worst_full,
                    unexp=len(unexp),nonstruct=len(nonstruct),Dp_zero=r['Dp_zero'],amax=max(V)))
    print('%-10s %-5s %2d %-5s %5d %8.5f %8.5f %8.3f %8.3f %8.3f %6.3f %6s  w=1:%s(%d fails, first %s) v<v_p(a!) at %d a (worst %d) unexp0=%d'%(
        K,f,p,'inert',n,slope,c,maxdef[0],C_ln,C_lp,hi,r['Dp_zero'],w1,len(fails1),fails1[:3],full_fail,worst_full,len(unexp)))
json.dump(out,open(sys.argv[2],'w'))
ii=[o for o in out if o['kind']=='inert']
print('INERT: all D^p t == 0 mod p:',all(o['Dp_zero'] for o in ii if o['f']=='t'),'; max C_ln = %.3f ; max C_logp = %.3f ; max w* = %.3f ; min excess v - ap/(p^2-1) = %.3f'%(
    max(o['C_ln'] for o in ii),max(o['C_logp'] for o in ii),max(o['wstar'] for o in ii),min(o['min_excess'] for o in ii)))
