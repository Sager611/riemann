# usage: python3 Riemann/Deep/Harvest6/k_est.py (stdlib only; < 1 s)
# rough float estimate of Gracia Llorente (Zenodo 24044336) eq. (20) majorant of E(x); NOT a certificate
import math
H0=3e12; h=math.log(H0); tau=3e-12
tab=[(20/32,21/32,9.387,8.290,147.0),(21/32,22/32,10.70,7.997,141.9),(22/32,23/32,12.11,7.704,136.7),(23/32,24/32,13.66,7.411,131.5),
     (24/32,25/32,15.35,7.118,126.3),(25/32,26/32,17.21,6.826,121.1),(26/32,27/32,19.31,6.533,115.9),(27/32,28/32,21.77,6.240,110.7)]
def Bs(s):
    for a,b,B1,B2,B3 in tab:
        if a<=s<=b+1e-15: return B1,B2,B3
    raise ValueError(s)
def W(s):
    al=3*(1-s)/(2-s); p=(7-5*s)/(2-s); q=2-al
    I3=28**(p-3)*H0**(-q)*(h**3/q+3*h*h/q**2+6*h/q**3+6/q**4)
    I2=H0**-2*(h*h/2+h/2+0.25); I1=H0**-2*(h/2+0.25)
    B1,B2,B3=Bs(s)
    return min(tau/2, 4*(B1*I3+B2*I2+B3*I1))
def E_bound(L, smax):
    x=math.exp(L); J=2048
    ss=[5/8+j/(4*J) for j in range(J+1)]
    Wv=[W(s) if s<=smax else 0.0 for s in ss]
    Wm=[]; m=float('inf')
    for w in Wv: m=min(m,w); Wm.append(m)
    c=[tau-2*Wm[0]]+[2*(Wm[j-1]-Wm[j]) for j in range(1,J)]+[2*Wm[J-1]]
    TJ=0.5*sum(cj*(math.exp((s-1)*L)+math.exp(-s*L)) for cj,s in zip(c,ss))
    f=1/L+1/L**2
    G0=13/25*(2.794e-12)**2
    main=f*2.134*math.exp(-L/2)
    return main, f*(1+1e-12)*TJ, f*2/(x*L), G0/L
for L,smax,lab in [(106,7/8,'GL check: L=106, s<=7/8'),(108.57,7/8,'L=108.57, s<=7/8'),(108.57,0.87338196661320752,'L=108.57, s<=B2')]:
    a,b,c2,d=E_bound(L,smax); tot=a+b+c2+d
    print(f'{lab}: main {a:.3e} zeros>H0 {b:.3e} G {d:.3e} total {tot:.3e} = 2^{math.log2(tot):.2f} -> K_max={math.floor(-math.log2(tot))-1}')
