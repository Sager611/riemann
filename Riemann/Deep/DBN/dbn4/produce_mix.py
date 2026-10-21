
"""Producer: classical-strip comparison with the extended single catalogue (evidence/) plus per-cell
two-shift mixed probes (same provider/height), each certified exactly on the fly (certify_mix.py).
Stepping/rounding rules identical to Planat's produce_classical_portfolio.py."""
import json, sys, time, math, gzip, hashlib
from pathlib import Path
from fractions import Fraction as F
import numpy as np
sys.path.insert(0,'.')
from certify_mix import certify_mix
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def ceilq(x,d): return F(-((-x.numerator*d)//x.denominator),d)
def exponential_upper(x):
    assert 0<x<1
    term=total=F(1)
    for k in range(1,19): term*=-x/k; total+=term
    return total
def logarithm_upper(z):
    assert 0<z<1
    return sum(((-1)**(k+1)*z**k/k for k in range(1,12)),F(0))
def cap_ok(cap,r,v):
    R=r*r;A=cap;B=cap*(2*R-2)+2;C=cap*(1+R)**2-2*(1+R)
    vals=[C,A*v*v+B*v+C];vx=-B/(2*A)
    if 0<vx<v: vals.append(A*vx*vx+B*vx+C)
    return min(vals)>=0
SH=['9/20','19/40','49/100','1/2','51/100','21/40','107/200','113/200','23/40','3/5']
names=['SOURCE_OLD_ENHANCED.json','SOURCE_RETUNED.json']
ev=Path('evidence'); S=[json.loads((ev/n).read_text()) for n in names]
# float machinery for the mixture search
zz=np.linspace(1e-7,1-1e-9,8001); DD=zz/(1-zz)
def Jf(r,v,B):
    R=r*r; D=DD; C=D+R+1+B; A=C*C+4*D*R-4*B
    D4=A*A-16*D*R*C*C; G=A*(C-2*B)-8*D*R*C
    return 2*G/D4
def Tf(v,B):
    Qs=(DD+v-B)**2+4*DD*B; return 4*(DD+v-B)/Qs
def s_r(w,r): R=r*r; return 2*(1+R-w)/((1+R-w)**2+4*R*w)
def cap_rational(r,rho):
    """smallest 1e-9-lattice M with s_r<=M on [0,rho^2], verified exactly"""
    w=np.linspace(0,float(rho)**2,4001); M=float(np.max(s_r(w,float(F(r)))))
    Mr=F(math.ceil(M*10**9),10**9)
    while not cap_ok(Mr,F(r),F(rho)**2): Mr+=F(1,10**9)
    return Mr
THETAS=[F(k,40) for k in range(1,40)]
def best_mixture(i,q,single_best):
    best=None
    for p in (0,1):
        h=F(S[p]['plan']['rows'][i]['height']); hf=float(h)
        rho=F(math.ceil(math.sqrt(float(q))/hf*(1+1e-7)*10**7),10**7)
        if rho>=F(48,100): continue
        v=float(rho)**2; B=v*(1-1e-13)
        T=Tf(v,B); Js=[Jf(float(F(r)),v,B) for r in SH]
        ells=[F(S[p]['rows'][i]['r_'+r.replace('/','_')]) for r in SH]
        caps=[cap_rational(r,rho) for r in SH]
        for a in range(10):
            for b in range(a+1,10):
                fa=float(ells[a]-caps[a]/h); fb=float(ells[b]-caps[b]/h)
                for th in THETAS:
                    t=float(th); lt=min(2.0,float(np.nanmin(T/((1-t)*Js[a]+t*Js[b]))))
                    c=2*lt/hf*((1-t)*fa+t*fb)
                    if best is None or c>best[0]: best=(c,p,a,b,th,lt,rho,h,caps[a],caps[b],ells[a],ells[b])
    if best is None or best[0]<=float(single_best)*(1+1e-7): return None
    c,p,a,b,th,lt,rho,h,Ma,Mb,la,lb=best
    lam=min(F(2),F(math.floor((lt-3e-6)*10**7),10**7))
    for attempt in range(8):
        l1=(1-th)*lam; l2=th*lam
        cert=certify_mix(SH[a],SH[b],rho,l1,l2)
        if cert['geometry']: break
        lam-=F(1,10**5)*2**attempt
    else: return None
    rate=2/h*(l1*(la-Ma/h)+l2*(lb-Mb/h))
    return dict(provider=p,r1=SH[a],r2=SH[b],rho=rho,lam1=l1,lam2=l2,cap1=Ma,cap2=Mb,floor1=la,floor2=lb,height=h,rate=rate,cert=cert)
def single_candidates(i):
    out=[]
    for pi in (0,1):
        h=F(S[pi]['plan']['rows'][i]['height'])
        for ki,k in enumerate(S[pi]['kernels']):
            r=F(k['r']);rho=F(k['rho']);lam=F(k['lambda']);cap=F(k['self_cap'])
            ell=F(S[pi]['rows'][i]['r_'+k['r'].replace('/','_')])
            rate=2*lam/h*(ell-cap/h)
            if rate>0 and abs(h*r)<2: out.append((rate,pi,ki,h*h*rho*rho,h,ell,cap))
    return sorted(out,key=lambda x:(-x[0],x[1],x[2]))
def produce(parts=32,mix_cells=345,label='classical_dbn4_mix'):
    dt=F(1,4000*parts);q0=F(1000001,1000000)**2;q=q0;qe=F(1,10**8)
    records=[];mixes=[];exp_cache={};t0=time.time()
    for i in range(608):
        cand=single_candidates(i)
        mix=None
        if i<mix_cells:
            sb=next((c for c in cand if q<=c[3]),None)
            mix=best_mixture(i,q,sb[0] if sb else 0)
            if mix:
                mix['id']=len(mixes); mix['cell']=i; mixes.append(mix)
        for j in range(parts):
            el=next((c for c in cand if q<=c[3]),None)
            use_mix = mix is not None and q<=mix['height']**2*mix['rho']**2 and (el is None or mix['rate']>el[0])
            if use_mix:
                rate=mix['rate']; key=('m',mix['id'])
                rec={'cell':i,'sub':j,'type':'mix','provider':mix['provider'],'mix':mix['id'],'q_in':str(q),'height':str(mix['height']),'rate':str(rate)}
            else:
                rate,pi,ki,thr,h,ell,cap=el; key=('s',i,pi,ki)
                rec={'cell':i,'sub':j,'type':'single','provider':pi,'kernel':ki,'q_in':str(q),'height':str(h),'floor':str(ell),'self_cap':str(cap),'rate':str(rate)}
            if not records: assert q0-(2+rate*q0)*F(1,10**8)>1
            now=F(i,4000)+j*dt
            ck=(key,i)
            if ck not in exp_cache: exp_cache[ck]=exponential_upper(rate*dt)
            nxt=ceilq((q+2/rate)*exp_cache[ck]-2/rate,2**96)
            records.append(rec)
            if nxt<=qe:
                z=rate*(q-qe)/(2+rate*qe);lag=logarithm_upper(z)/rate
                assert 0<lag<=dt
                rec['terminal_lag_upper']=str(lag)
                upper=ceilq(now+lag+qe/2,10**15); assert upper<F(41,250)
                mixout=[{k:(str(v) if isinstance(v,F) else v) for k,v in m.items()} for m in mixes]
                hist={'schema':'DBN4_MIXED_PORTFOLIO_V1','label':label,'parts':parts,
                      'source_files':names,'source_hashes':[sha(ev/n) for n in names],
                      'initial_q':str(q0),'initial_strip_q':'1','initial_delay':'1/100000000','terminal_q':str(qe),
                      'rounding_grid_bits':96,'upper_exact':str(upper),'upper_decimal':format(float(upper),'.15f'),
                      'steps':len(records),'mixed_kernels':mixout,'records':records}
                Path('histories').mkdir(exist_ok=True)
                with gzip.open('histories/'+label+'.json.gz','wt') as f: json.dump(hist,f,separators=(',',':'),sort_keys=True)
                print('PRODUCED',label,hist['upper_decimal'],len(records),'mixed kernels',len(mixes),'mixed steps',sum(1 for r in records if r['type']=='mix'),'%.0fs'%(time.time()-t0),flush=True)
                return hist
            assert qe<nxt<q
            rec['q_out']=str(nxt);q=nxt
        if i%50==0: print('cell',i,'q=%.6f'%float(q),'mixes',len(mixes),'%.0fs'%(time.time()-t0),flush=True)
if __name__=='__main__':
    produce()
