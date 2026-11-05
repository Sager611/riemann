
# test_nsolve.py -- random promised norm equations; checks every claim along the way
import random, sys, time
from ffpoly import *
from nsolve import NormSolver, Frac
import numpy as np

def prim_root_unity(q,p):
    for g in range(2,p):
        z=pow(g,(p-1)//q,p)
        if z!=1: return z
def primary_gen(q,p):
    n=p-1; s=0; t=n
    while t%q==0: t//=q; s+=1
    for g in range(2,p):
        w=pow(g,t,p)
        if pow(w,q**(s-1),p)!=1: return w
def rand_sqfree(d,p):
    while True:
        f=[random.randrange(p) for _ in range(d)]+[1]
        if pgcd(f,pderiv(f,p),p)==[1]: return f

def norm_poly(F, coeffs, q, p):
    # norm of sum_j coeffs[j] Y^j (polynomial coefficients) as polynomial
    ns = NormSolver.__new__(NormSolver)
    ns.p, ns.q, ns.F = p, q, F
    h = [Frac(c, [1], p) for c in coeffs]
    G = ns.norm_of(h)
    assert G.d == [1]
    return G.n

def run(seed, cases, Mextra=0, verbose=False):
    random.seed(seed)
    stats = {'ok':0,'fail':0,'dimS':[], 'groups':0, 'a0_in_B':0}
    for (p,q,N,dh,neg) in cases:
        zeta=prim_root_unity(q,p); omega=primary_gen(q,p)
        F=rand_sqfree(N,p)
        # random h0 = (sum a_j Y^j) / c(X)
        while True:
            if dh<0:
                coeffs=[[random.randrange(1,p)],[random.randrange(0,p)]]+[[] for _ in range(q-2)]
            else:
                coeffs=[norm([random.randrange(p) for _ in range(random.randrange(0,dh+2))]) for _ in range(q)]
            if any(coeffs): break
        if random.random()<0.3:
            coeffs=[[]]*(q-1)+[coeffs[0] or [1]]  # monomial times Y^{q-1}
        A=norm_poly(F,coeffs,q,p)
        cden=[1]
        if dh>=0 and random.random()<0.5:
            cden=norm([random.randrange(p) for _ in range(random.randrange(1,3))]+[1])
        B=ppow(cden,q,p)
        k=random.randrange(1,p)
        A=pscal(A,pow(k,q,p),p)
        if neg:
            # multiply by a non-q-th power constant: generally not a norm (test of failure behaviour)
            A=pscal(A,omega,p)
        try:
            t0=time.time()
            ns=NormSolver(p,q,zeta,omega,F,A,B)
            if Mextra: ns.M += Mextra
            h=ns.run()
            Nh=ns.norm_of(h)
            G=Frac(ns.A,ns.B,p)
            good = (Nh == G)
            info=ns.info
            stats['dimS'].append(info['dimS'])
            stats['a0_in_B'] += int(info['a0_in_B'])
            ok = good and info['idempotent'] and info['trd_e']==Frac([1],[1],p) and info['a0_in_B']
            if verbose or not ok:
                print('case',(p,q,N,dh,neg),'D',ns.D,'M',ns.M,'groups',info['groups'],'unk',info['unknowns'],'eq',info['equations'],'dimS',info['dimS'],'dimB',info['dimB'],'dimN',info['dimN'],'dimU',info['dimU'],'a0inB',info['a0_in_B'],'idem',info['idempotent'],'trd1',info['trd_e']==Frac([1],[1],p),'N(h)=G',good,'%.1fs'%(time.time()-t0),flush=True)
            stats['ok' if ok else 'fail'] += 1
        except Exception as ex:
            stats['fail'] += 1
            print('EXC case',(p,q,N,dh,neg),type(ex).__name__,str(ex)[:120],flush=True)
    return stats

if __name__=='__main__':
    seed=int(sys.argv[1]); which=sys.argv[2]
    if which=='smoke':
        cases=[(10009,3,3,1,False)]*3+[(10111,5,5,0,False)]*1
        print(run(seed,cases,verbose=True))
    elif which=='q3':
        cases=[(10009,3,3,random.Random(i).randrange(0,3),False) for i in range(40)]+[(10009,3,6,1,False) for i in range(10)]
        print(run(seed,cases))
    elif which=='small':
        cases=[(pp,3,3,1,False) for pp in [7,13,19,31,37,43] for _ in range(5)]
        print(run(seed,cases))
    elif which=='q5':
        cases=[(10111,5,5,random.Random(i).randrange(0,2),False) for i in range(8)]
        print(run(seed,cases))
    elif which=='q7':
        cases=[(10093,7,7,0,False) for i in range(2)]
        print(run(seed,cases,verbose=True))
    elif which=='neg':
        cases=[(10009,3,3,1,True) for i in range(10)]
        print(run(seed,cases,verbose=True))

    elif which=='batch':
        import random as _r
        cases=[(10009,3,3,0,False) for i in range(16)]+[(10009,3,3,1,False) for i in range(4)]+[(1009,3,3,0,False) for i in range(4)]+[(331,3,3,0,False) for i in range(4)]
        print('q3 batch',run(seed,cases,verbose=True),flush=True)
        cases=[(10009,3,3,0,True) for i in range(6)]
        print('neg batch',run(seed+1,cases,verbose=True),flush=True)
        cases=[(10009,3,3,0,False) for i in range(4)]
        s1=run(seed+2,cases); s2=run(seed+2,cases,Mextra=6)
        print('M stability dimS',s1['dimS'],'vs M+6',s2['dimS'],flush=True)
    elif which=='q5small':
        cases=[(10111,5,5,-1,False) for i in range(2)]
        print('q5 batch',run(seed,cases,verbose=True),flush=True)
    elif which=='mstab':
        cases=[(10009,3,3,1,False) for i in range(6)]
        s1=run(seed,cases); s2=run(seed,cases,Mextra=6)
        print('dimS M_eff',s1['dimS'],'dimS M_eff+6',s2['dimS'])

