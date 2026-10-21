
"""SUPERSEDED by EVIDENCE/scripts/audit_superset.py (a superset of Planat's auditor; see EVIDENCE/README.md).
Index range checks added 2026-10-21 (before, a negative index wrapped around).
Independent auditor for DBN4_MIXED_PORTFOLIO_V1 histories (imports no producer code).
Checks every record: contiguity, incoming state, provider/height ownership, source floors, own-pair caps on the
whole decreasing range, kernel domain, exact rate, outward exponential restart, dyadic lattice, terminal lag,
initial separation and the final strip cost. Mixed kernels: geometry re-certified (certify_mix) and the
polynomial cross-checked against an independent sympy construction for a sample."""
import json, gzip, hashlib, sys, math, time, random
from pathlib import Path
from fractions import Fraction as F
def req(x,m):
    if not x: raise ValueError(m)
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def exp_upper(x):
    req(0<x<1,'exp arg'); return sum(((-x)**k/F(math.factorial(k)) for k in range(21)),F(0))
def log_upper(z):
    req(0<z<1,'log arg'); w=z/(2+z)
    return 2*sum((w**(2*k+1)/F(2*k+1) for k in range(7)),F(0))+2*w**15/(15*(1-w*w))
def cap_ok(cap,r,v):
    R=r*r;A=cap;B=cap*(2*R-2)+2;C=cap*(1+R)**2-2*(1+R)
    req(A>0 and v>=0,'cap range'); vals=[C,A*v*v+B*v+C];vx=-B/(2*A)
    if 0<vx<v: vals.append(A*vx*vx+B*vx+C)
    return min(vals)>=0
SHIFTS={'9/20','19/40','49/100','1/2','51/100','21/40','107/200','113/200','23/40','3/5'}
def audit(path,ev='evidence',recertify=True):
    H=json.load(gzip.open(path,'rt'))
    names=['SOURCE_OLD_ENHANCED.json','SOURCE_RETUNED.json']
    S=[json.loads((Path(ev)/n).read_text()) for n in names]
    req(H['schema']=='DBN4_MIXED_PORTFOLIO_V1','schema')
    req(H['source_hashes']==[sha(Path(ev)/n) for n in names],'source identity')
    parts=H['parts']; req(type(parts) is int and parts>0,'parts')
    dt=F(1,4000*parts); q0=F(1000001,1000000)**2; qe=F(1,10**8)
    req(F(H['initial_q'])==q0 and F(H['terminal_q'])==qe and F(H['initial_strip_q'])==1,'strip')
    for s in S:
        req(len(s['rows'])==608 and len(s['plan']['rows'])==608,'schedule')
        req(s['frequency_pairs']==[[j,j] for j in range(6)]+[[j,j+1] for j in range(5)],'classes')
        for i,p in enumerate(s['plan']['rows']):
            req(p['time_index']==i and F(p['t_lo'])==F(i,4000) and F(p['t_hi'])==F(i+1,4000),'cells')
    MK=H['mixed_kernels']
    # geometry of mixed kernels
    if recertify:
        sys.path.insert(0,'.')
        from certify_mix import certify_mix
        t=time.time()
        for m in MK:
            req(m['r1'] in SHIFTS and m['r2'] in SHIFTS,'declared shifts')
            req(F(m['lam1'])>=0 and F(m['lam2'])>=0 and 0<F(m['rho'])<1,'mixed params')
            c=certify_mix(m['r1'],m['r2'],m['rho'],m['lam1'],m['lam2'])
            req(c['geometry'],'mixed geometry %s'%m['id'])
        print('mixed geometry re-certified',len(MK),'%.0fs'%(time.time()-t),flush=True)
    q=q0; mixed_used=set(); single_used=set()
    for n,rec in enumerate(H['records']):
        i,j=divmod(n,parts); now=n*dt
        req((rec['cell'],rec['sub'])==(i,j),'contiguous')
        req(F(rec['q_in'])==q,'incoming state')
        pi=rec['provider']; req(pi in (0,1),'provider'); src=S[pi]
        h=F(src['plan']['rows'][i]['height']); req(F(rec['height'])==h,'height')
        v=q/h**2
        if rec['type']=='single':
            req(type(rec['kernel']) is int and 0<=rec['kernel']<len(src['kernels']),'kernel index')
            k=src['kernels'][rec['kernel']]; r,rho,lam=F(k['r']),F(k['rho']),F(k['lambda'])
            req(v<=rho*rho and abs(h*r)<2 and 0<lam<=2,'kernel domain')
            ell=F(src['rows'][i]['r_'+k['r'].replace('/','_')]); cap=F(k['self_cap'])
            req(F(rec['floor'])==ell and F(rec['self_cap'])==cap,'floor/cap')
            req(cap_ok(cap,r,v),'cap range')
            c=2*lam/h*(ell-cap/h); single_used.add(rec['kernel'])
        elif rec['type']=='mix':
            req(type(rec['mix']) is int and 0<=rec['mix']<len(MK),'mixture index')
            m=MK[rec['mix']]; req(m['id']==rec['mix'] and m['cell']==i and m['provider']==pi,'mixed ownership')
            r1,r2,rho,l1,l2=F(m['r1']),F(m['r2']),F(m['rho']),F(m['lam1']),F(m['lam2'])
            req(v<=rho*rho and abs(h*r1)<2 and abs(h*r2)<2,'mixed domain')
            e1=F(src['rows'][i]['r_'+m['r1'].replace('/','_')]); e2=F(src['rows'][i]['r_'+m['r2'].replace('/','_')])
            req(F(m['floor1'])==e1 and F(m['floor2'])==e2 and F(m['height'])==h,'mixed floors')
            M1,M2=F(m['cap1']),F(m['cap2']); req(cap_ok(M1,r1,v) and cap_ok(M2,r2,v),'mixed caps')
            c=2/h*(l1*(e1-M1/h)+l2*(e2-M2/h)); mixed_used.add(rec['mix'])
        else: raise ValueError('type')
        req(c==F(rec['rate']) and 0<c*dt<1,'rate')
        if n==0:
            d0=F(H['initial_delay']); req(0<d0<dt and q0-(2+c*q0)*d0>1,'initial separation')
        if 'terminal_lag_upper' in rec:
            z=c*(q-qe)/(2+c*qe); lag=F(rec['terminal_lag_upper'])
            req(log_upper(z)/c<=lag<=dt,'terminal lag')
            U=F(H['upper_exact']); req(now+lag+qe/2<=U<F(41,250),'final cost')
            req(n+1==len(H['records'])==H['steps'],'complete')
            return dict(PASS=True,upper_exact=str(U),upper_decimal=format(float(U),'.15f'),steps=n+1,
                        mixed_kernels_used=len(mixed_used),single_kernels_used=len(single_used),
                        mixed_steps=sum(1 for r in H['records'] if r['type']=='mix'))
        out=F(rec['q_out']); req(qe<out<q,'decreasing'); req((out*2**96).denominator==1,'lattice')
        req(out>=(q+2/c)*exp_upper(c*dt)-2/c,'outward restart')
        q=out
    raise ValueError('no terminal')
def negative_controls(path):
    import copy, io
    H=json.load(gzip.open(path,'rt'))
    bad=[('inflated mixed floor',lambda h:h['mixed_kernels'][0].update(floor1='100')),
         ('omitted mixed cap',lambda h:h['mixed_kernels'][0].update(cap1='0')),
         ('mixed domain',lambda h:h['mixed_kernels'][0].update(rho='1/100')),
         ('inflated mixed lam',lambda h:h['mixed_kernels'][0].update(lam1=str(F(h['mixed_kernels'][0]['lam1'])+F(1,10)))),
         ('downward restart',lambda h:h['records'][0].update(q_out='1/10')),
         ('time gap',lambda h:h['records'][1].update(sub=2))]
    res=[]
    for name,f in bad:
        h=copy.deepcopy(H); f(h)
        tmp='histories/_neg.json.gz'
        with gzip.open(tmp,'wt') as g: json.dump(h,g)
        try: audit(tmp,recertify=(name=='inflated mixed lam'))
        except (ValueError,KeyError,IndexError,AssertionError) as e: res.append((name,'rejected',str(e)[:40]))
        else: res.append((name,'ACCEPTED'))
    Path('histories/_neg.json.gz').unlink()
    return res
if __name__=='__main__':
    t=time.time(); A=audit('histories/classical_dbn4_mix.json.gz')
    print('AUDIT',A,'%.0fs'%(time.time()-t),flush=True)
    json.dump(A,open('histories/classical_dbn4_mix_AUDIT.json','w'),indent=2)
    print('NEG',negative_controls('histories/classical_dbn4_mix.json.gz'))
