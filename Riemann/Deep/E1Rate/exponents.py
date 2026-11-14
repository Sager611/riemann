# e1rate exponents: Conrey-structure (Vaughan + DI 1984 Lemma 1 + Weil) bound for SL's dual sum in windows H=T^theta.
# All quantities are exponents of T; error terms relative to the main term H. See DERIVATION.md section 3.
import sys
def gDI(a,b,n,v):  # DI Lemma 1 (Baluyot L3.1.6) relative gain: {..}/(ABNV)^{1/2}
    m=max(n,a)
    return -(a+b+n+v)/2+max((v+b)/2, m/4+max(v+b+m+max(v,2*a), n+2*a+2*b)/4)
def bounds(th,pv,kq,n):
    # pv: modulus variable exponent, kq: decomposed variable exponent, n: dual length exponent
    x0=(n+1-pv-kq)/2
    etr=0.75*n+0.25*(pv+kq)-0.25
    esep=etr+max(0.0,x0)/2
    return etr,esep,x0
def eII(th,pv,a,b,n):
    etr,esep,x0=bounds(th,pv,a+b,n); return esep+gDI(a,b,n,pv)
def eI(th,pv,kq,a,n):
    etr,esep,x0=bounds(th,pv,kq,n); return etr+max(-pv, pv/2+max(0.0,x0)/2-(kq-a))
def Bdec(th,pv,kq,n,us,K=16):
    etr=bounds(th,pv,kq,n)[0]; best=9.0
    for u in us:
        if kq<=u: worst=etr
        else:
            worst=-9.0
            for i in range(K+1):   # piece 2: a<=U^2 arbitrary, smooth rest
                a=min(2*u,kq)*i/K
                worst=max(worst,min(etr,eI(th,pv,kq,a,n),eII(th,pv,a,kq-a,n),eII(th,pv,kq-a,a,n)))
            if kq-u>u:
                for i in range(K+1):   # piece 3: a,b>U
                    a=u+(kq-2*u)*i/K
                    worst=max(worst,min(etr,eII(th,pv,a,kq-a,n),eII(th,pv,kq-a,a,n)))
        best=min(best,worst)
    return best
def E(th,nu,G=16,withK=True):
    us=[nu*j/40 for j in range(0,21)]; worst=-9.0; arg=None
    for i in range(G+1):
        for j in range(G+1):
            P=nu*i/G; Q=nu*j/G; n0=P+Q+1-2*th
            if n0<0: continue
            for n in (n0, 0.75*n0, 0.5*n0, 0.25*n0, 0.0):
                etr=bounds(th,P,Q,n)[0]
                b=etr if not withK else min(etr,Bdec(th,P,Q,n,us),Bdec(th,Q,P,n,us))
                if b>worst: worst,arg=b,(round(P,4),round(Q,4),round(n,4))
    return worst,arg
def nuK(th,withK=True,lo=0.0,hi=0.75):
    for _ in range(17):
        mid=(lo+hi)/2
        if E(th,mid,withK=withK)[0]<0: lo=mid
        else: hi=mid
    return lo
if __name__=="__main__":
    for th in [float(x) for x in sys.argv[1:]]:
        a=nuK(th,withK=False); b=nuK(th)
        print(f"theta={th:.4f} nu_E1(grid)={a:.5f} (3th-1)/4={(3*th-1)/4:.5f} nu_K(grid)={b:.5f} (4/7)(2th-1)={4*(2*th-1)/7:.5f} worst@nuK-0.002={E(th,b-0.002)}",flush=True)
