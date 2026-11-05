
# p3_lifting.py -- Lemma local-lifting (uniform lifting classes) for small primes
import sys
sys.path.insert(0,'.')
from p3_local import *

def build_hl(h,prefer,ps,Lmax):
    J=12; z={}; m={}
    for p in ps:
        zp=unit_root_padic(h,p,J+2,prefer)
        exact = prefer is not None and any(zz%p and peval(h,zz)==0 and zz%(p**(J+2))==zp for zz in prefer)
        z[p]=zp; m[p]=mult_at(h,[zz for zz in prefer if zz%p and peval(h,zz)==0][0] if exact else zp,p,J,exact)
    def lam(n):
        r=1
        for p in ps:
            if n%p==0: r*=p**(vp(n,p)*m[p])
        return r
    def rl(l):
        if l==1: return 0
        r=0; mod=1
        for p in ps:
            if l%p: continue
            pe=p**vp(l,p); t=((z[p]-r)*pow(mod,-1,pe))%pe; r=r+mod*t; mod*=pe
        r%=l
        return r-l if r>0 else 0
    out={}
    for l in range(1,Lmax+1):
        r=rl(l); num=compose_affine(h,r,l); L=lam(l)
        out[l]=(r,[c//L for c in num])
    return out

def run(name,h,prefer,Lmax=30):
    ps=primes_upto(Lmax)
    fam=build_hl(h,prefer,ps,Lmax)
    k=len(h)-1
    res={}
    bad=0
    for p in [2,3,5,7]:
        Sp=0; best={}
        for l,(r,hl) in fam.items():
            dh=deriv(hl)
            cands=[x for x in range(1,k*p+1) if (r+l*x)%p]
            s=min(vp(peval(dh,x),p) for x in cands)
            x0=[x for x in cands if vp(peval(dh,x),p)==s][0]
            best[l]=(x0,s); Sp=max(Sp,s)
        Ep=2*Sp+1
        for l,(r,hl) in fam.items():
            x0,s=best[l]
            mod=p**(Ep-s)
            # class constancy mod p^Ep and unit condition on the class (check representatives)
            base=peval(hl,x0)%(p**Ep)
            for t in range(0,p**2):
                y=x0+t*mod
                if (r+l*y)%p==0: bad+=1
                if peval(hl,y)%(p**Ep)!=base: bad+=1
            # uniform pushforward for j = Ep+1, Ep+2
            for j in (Ep+1,Ep+2):
                cnt={}
                n=p**(j-Ep+s)    # residues of y in the class mod p^j
                for t in range(n):
                    y=x0+t*mod
                    v=peval(hl,y)%(p**j)
                    cnt[v]=cnt.get(v,0)+1
                target=[v for v in range(p**j) if v%(p**Ep)==base]
                if set(cnt)!=set(target) or len(set(cnt.values()))!=1: bad+=1
        res[p]=(Sp,Ep)
    return res,bad

if __name__=='__main__':
    h8=pmul(pmul([-2,0,1],[-17,0,1]),[-34,0,1])
    for name,h,pref in [('(x-1)^2',[1,-2,1],[1]),('x^2-1',[-1,0,1],[1,-1]),('x^3-x',[0,-1,0,1],[1,-1]),('(x-1)^2(x+1)',[1,-1,-1,1],[1,-1]),('(x-1)^3',[-1,3,-3,1],[1]),('h8',h8,None)]:
        print(name,'(S_p,E_p) per p, violations:',run(name,h,pref),flush=True)

