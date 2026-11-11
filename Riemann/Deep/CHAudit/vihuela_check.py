
# vihuela_check.py: direct numerical check of Prop tritura (adago) and the assembled bound of Prop vihuela at T=1e7 and T=3e12+1+pi/3.
import mpmath as mp
mp.mp.dps=30
def F(u): return 1/mp.pi-(1-u)*mp.cot(mp.pi*(1-u))
C1=mp.nsum(lambda n: mp.zeta(2*n)*(1/(2*n)-2/(2*n+1)+1/(2*n+2)),[1,mp.inf])
C2=mp.nsum(lambda n: mp.zeta(2*n)*(1/(2*n)**2-2/(2*n+1)**2+1/(2*n+2)**2),[1,mp.inf])
t0=mp.mpf(20000)
for T in [mp.mpf(10)**7, mp.mpf(10)**9, 3*mp.mpf(10)**12+1+mp.pi/3]:
    eps=t0/T
    # (2pi/T) int_{t0}^T sqrt(F(t/T)^2+(1-t/T)^2) log(t/2pi) dt = 2pi int_eps^1 sqrt(F(u)^2+(1-u)^2) log(uT/2pi) du
    f=lambda u: mp.sqrt(F(u)**2+(1-u)**2)*mp.log(u*T/(2*mp.pi))
    pts=[eps]+[mp.mpf(10)**k for k in range(int(mp.floor(mp.log10(eps)))+1,0)]+[mp.mpf(1)]
    LHS=2*mp.pi*mp.quad(f,pts)
    y=mp.log(T/(2*mp.pi)); y0=mp.log(t0/(2*mp.pi))
    RHS=y**2-y0**2-2*C1*mp.log(mp.e*T/(2*mp.pi))+2*C2
    # assembled vihuela bound for sigma in {0,1}: adar main + err + salmon low part
    for sig in [0,1]:
        a=abs(sig-mp.mpf(1)/2)
        err1=2*(mp.mpf(2)/5*mp.log(t0)+mp.mpf(21)/5)+(mp.log(mp.e*t0/(2*mp.pi))/mp.pi+4/(5*t0)*(mp.log(t0)+mp.mpf(41)/4))*a
        err2=(mp.mpf('2.78')*a+1)*y-mp.mpf(2)/5+mp.pi**2*t0/T*(mp.mpf(2)/5*mp.log(t0)+4)
        P1=2*C1*y+2*(C1-C2)
        adar=(y**2-y0**2-P1)/(2*mp.pi)+err1/t0+err2/T
        low=mp.mpf('10.3193172')+(1+(mp.mpf('2.78')*a+1)/T)*(t0/T)*y0   # 2 sum 1/gamma (our zeros) + c t0/T log(t0/2pi)
        tot=adar+low; claim=y**2/(2*mp.pi)-mp.mpf('1.01')/(6*mp.pi)*y
        print('T=%s sigma=%d: tritura LHS=%s RHS=%s (ok=%s) | assembled bound=%s  vihuela claim=%s  margin=%s'%(mp.nstr(T,8),sig,mp.nstr(LHS,10),mp.nstr(RHS,10),LHS<=RHS,mp.nstr(tot,10),mp.nstr(claim,10),mp.nstr(claim-tot,5)))
