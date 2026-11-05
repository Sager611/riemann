
# p4_checks.py -- finite checks for "Weighted dilation graphs, smooth shifted primes and totient fibers"
import math, itertools

def determinant_bijection(E0=12, M0=60, H0=60):
    # forward set: (e,n,m,mp,h,hp) positive, m!=mp, m n = 2 e h + 1, mp n = 2 e hp + 1
    fwd=set()
    for e in range(1,E0+1):
        for h in range(1,H0+1):
            v=2*e*h+1
            for m in range(1,min(M0,v)+1):
                if v%m: continue
                n=v//m
                # mp with mp n = 2 e hp + 1
                for hp in range(1,H0+1):
                    w=2*e*hp+1
                    if w%n: continue
                    mp=w//n
                    if mp!=m and mp<=M0: fwd.add((e,n,m,mp,h,hp))
    # map forward -> (k, m, mp, h, hp) via identities, check them
    img=set(); bad=0
    for (e,n,m,mp,h,hp) in fwd:
        if (mp-m)%(2*e): bad+=1; continue
        k=(mp-m)//(2*e)
        if hp-h!=k*n: bad+=1
        if m*hp-mp*h!=k: bad+=1
        img.add((k,m,mp,h,hp))
    # converse: all (k,m,mp,h,hp) with mp*h=z, m*hp=z+k, mp=m mod 2|k|, e=(mp-m)/(2k)>0, gcd(m,k)=1
    conv=set(); conv_nocop=0
    for m in range(1,M0+1):
        for mp in range(1,M0+1):
            if mp==m: continue
            for h in range(1,H0+1):
                for hp in range(1,H0+1):
                    k=m*hp-mp*h
                    if k==0 or (mp-m)%(2*abs(k)): continue
                    e=(mp-m)//(2*k)
                    if e<=0 or e>E0: continue
                    if math.gcd(m,abs(k))!=1:
                        # paper: coprimality necessary for converse; check whether it would fail
                        if (2*e*h+1)%m: conv_nocop+=1
                        else: conv.add((k,m,mp,h,hp))
                        continue
                    if (2*e*h+1)%m: bad+=1; continue
                    n=(2*e*h+1)//m
                    if hp-h!=k*n or mp*n!=2*e*hp+1: bad+=1
                    conv.add((k,m,mp,h,hp))
    # compare (restricted to ranges)
    missing=img-conv; extra=conv-img
    return len(fwd),len(img),len(conv),len(missing),len(extra),bad,conv_nocop

def block_sieve(z=13,h=2,N=200000):
    ps=[p for p in range(2,z+1) if all(p%q for q in range(2,int(p**0.5)+1))]
    blocks={}
    for p in ps:
        j=0
        while not (z**(2.0**(-j-1))<p<=z**(2.0**(-j))): j+=1
        blocks.setdefault(j,[]).append(p)
    J=sorted(blocks)
    def Uj_val(j,t): # sum_{a<=h_j} (-1)^a C(t,a)
        hj=(j+1)*h
        return sum((-1)**a*math.comb(t,a) for a in range(0,min(t,hj)+1))
    def Ej_val(j,t):
        hj=(j+1)*h
        return math.comb(t,hj+1)
    bad=0; maxlevel=0
    for n in range(1,N+1):
        ts={j:sum(1 for p in blocks[j] if n%p==0) for j in J}
        U=1
        for j in J: U*=Uj_val(j,ts[j])
        corr=0
        for j in J:
            prod=Ej_val(j,ts[j])
            for k2 in J:
                if k2!=j: prod*=Uj_val(k2,ts[k2])
            corr+=prod
        V=U-corr
        I=1 if all(n%p for p in ps) else 0
        if not (V<=I<=U and U>=0): bad+=1
    # coefficient support/level: monomials of U are products over blocks of <= h_j primes
    lvl=0.0
    for j in J:
        hj=(j+1)*h
        largest=sorted(blocks[j])[-hj:] if len(blocks[j])>=hj else blocks[j]
        lvl+=sum(math.log(p) for p in largest)
    return N,bad,round(lvl/math.log(z),3),4*h+2

if __name__=='__main__':
    print('determinant bijection (fwd, img, conv, missing, extra, bad, noncoprime-failures):',determinant_bijection(),flush=True)
    for z,h in [(13,2),(30,2),(30,4)]:
        print('block sieve z=%d h=%d (N, violations, level/log z, bound 4h+2):'%(z,h),block_sieve(z,h,60000),flush=True)

