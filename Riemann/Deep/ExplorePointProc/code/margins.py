
from mpmath import mp, ci, si, euler, log, cos, pi, mpf, floor
mp.dps=20
def Cin(x): return euler+log(x)-ci(x)
def A(L): return Cin(2*pi*L)/pi**2
def varGUE(L):
    c=2*pi*L
    return (Cin(c)+1-cos(c)+c*(pi/2-si(c)))/pi**2
def tail_conj(L): return varGUE(L)-A(L)   # int_{|a|>1} W(a) da with F=1
print(' L     A(L)    VarGUE   th(1-th)  margin=th(1-th)-A   conj tail   ratio margin/conjtail')
for L in [0.5,0.9,1.0,1.1,1.2,1.25,1.3,1.4,1.5,1.6,1.7,1.75,1.8,1.9,2.0,2.25,2.5,2.75,3.0,3.5,4.5]:
    L=mpf(L); th=L-floor(L); y=th*(1-th)
    a=A(L); v=varGUE(L); t=tail_conj(L)
    print('%5.2f  %.5f  %.5f  %.5f  %+.5f   %.5f   %.3f'%(L,a,v,y,y-a,t,(y-a)/t))
