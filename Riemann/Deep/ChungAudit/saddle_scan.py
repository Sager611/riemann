# saddle_scan.py: solve Chung's saddle equation (App. A.2/A.4) at ACTUAL interior roots and test the
# uniform claims of Thm 5.3, Lemma 6.2b, Lemma 6.4, Lemma A.4. Double precision is enough (O(1) quantities).
# Run: python3 saddle_scan.py ROOTSJSON
import sys, json, cmath, math
pi=math.pi
def model(u, om, M, k):
    g=om*om/M; eta=om*om/(M*M); nu2=k*(k+1)/(om*om)
    W=2*pi*cmath.exp(2*u)
    V=1-eta*u*u-nu2/(u*u); Vp=-2*eta*u+2*nu2/u**3; Vpp=-2*eta-6*nu2/u**4
    e2=cmath.exp(-2*u); e4=cmath.exp(-4*u)
    Q=1-15/(4*pi)*e2+15/(8*pi*pi)*e4; Qp=15/(2*pi)*e2-15/(2*pi*pi)*e4; Qpp=-15/pi*e2+30/(pi*pi)*e4
    return g,eta,nu2,W,V,Vp,Vpp,Q,Qp,Qpp
def F_and_Fu(u, om, M, k, Sref):
    g,eta,nu2,W,V,Vp,Vpp,Q,Qp,Qpp=model(u,om,M,k)
    S=cmath.sqrt(V)
    if Sref is not None and abs(S-Sref)>abs(-S-Sref): S=-S
    E=Qp/Q-Vp/(4*V)
    F=W-6.5-k/u-g*u-1j*om*S-E
    LQp=(Qpp*Q-Qp*Qp)/(Q*Q)
    Fu=2*W+k/(u*u)-g-1j*om*Vp/(2*S)-LQp+(Vpp*V-Vp*Vp)/(4*V*V)
    return F,Fu,S,W,V,g,eta,nu2
def newton(u,om,M,k,Sref,it=60):
    for _ in range(it):
        F,Fu,S,*_=F_and_Fu(u,om,M,k,Sref)
        du=-F/Fu
        # damp
        if abs(du)>0.2: du=du*0.2/abs(du)
        u=u+du; Sref=S
        if abs(du)<1e-14: break
    F,Fu,S,W,V,g,eta,nu2=F_and_Fu(u,om,M,k,Sref)
    return u,abs(F),Fu,S,W,V,g,eta,nu2
def solve(om_t,M_t,k_t,steps=200):
    # homotopy from anchor (om=12.272, M=35, k=0) to target, continuous in (om, M, k)
    u=complex(0.470776,0.486693); Sref=None
    om0,M0,k0=12.272,35.0,0.0
    for i in range(1,steps+1):
        t=i/steps
        om=om0+(om_t-om0)*t; M=M0+(M_t-M0)*t; k=k0+(k_t-k0)*t
        u,res,Fu,S,W,V,g,eta,nu2=newton(u,om,M,k,Sref)
        Sref=S
    return u,res,Fu,S,W,V,g,eta,nu2
def consts(u,Fu,S,W,V,g,eta,nu2,om,M,k):
    s,t=u.real,u.imag
    Z=Fu/2
    Aeff=W.real-g/2
    a1=W/2
    Phia=cmath.phase(W-4.5-6/(W-3))
    thu=math.atan2(t,s)
    dTh0=Phia-thu
    kT=abs(Z)/abs(W)
    alpha=eta*abs(u)**2; beta=nu2/abs(u)**2
    Ka0=abs(1/u-2)
    return dict(sigma=s,tau=t,Reff=abs(W),Zabs=abs(Z),Acore=Z.real,Aeff=Aeff,kT=kT,dTheta0=dTh0,Phia=Phia,
                alpha=alpha,beta=beta,ReS=S.real,absV=abs(V),Ka0=Ka0,ReW=W.real,ImW=W.imag)
if __name__=='__main__':
    # anchor check
    u,res,Fu,S,W,V,g,eta,nu2=solve(12.272,35,0,steps=1)
    c=consts(u,Fu,S,W,V,g,eta,nu2,12.272,35,0)
    print('ANCHOR', {k_:round(v,6) for k_,v in c.items()}, 'res',res)
    R=json.load(open(sys.argv[1]))
    claims={'sigma>=0.470776':lambda c:c['sigma']>=0.470776-1e-9,'tau>=0.420':lambda c:c['tau']>=0.420,
            'tau<pi/4':lambda c:c['tau']<pi/4,'Acore>=6.10':lambda c:c['Acore']>=6.10,'Aeff>=6.40':lambda c:c['Aeff']>=6.40,
            'kT>=0.72':lambda c:c['kT']>=0.72,'dTheta0>=0.4388':lambda c:c['dTheta0']>=0.4388,'Phia<=1.2604':lambda c:c['Phia']<=1.2604,
            'alpha<=0.7945':lambda c:c['alpha']<=0.7945,'beta<=1.10':lambda c:c['beta']<=1.10,'|V|>=0.22@us':lambda c:c['absV']>=0.22,
            'Zabs>=13.83':lambda c:c['Zabs']>=13.83,'ReS>=0.420':lambda c:c['ReS']>=0.420}
    viol={k_:[] for k_ in claims}; worst={k_:None for k_ in claims}
    keyf={'sigma>=0.470776':'sigma','tau>=0.420':'tau','tau<pi/4':'tau','Acore>=6.10':'Acore','Aeff>=6.40':'Aeff','kT>=0.72':'kT',
          'dTheta0>=0.4388':'dTheta0','Phia<=1.2604':'Phia','alpha<=0.7945':'alpha','beta<=1.10':'beta','|V|>=0.22@us':'absV','Zabs>=13.83':'Zabs','ReS>=0.420':'ReS'}
    bad=[]; out={}
    for key,val in R.items():
        d,k=map(int,key.split(','))
        M=4*d+2*k-1
        for m in range(1,d-2):
            om=float(val['omega'][m-1])
            u,res,Fu,S,W,V,g,eta,nu2=solve(om,M,k)
            if res>1e-8: bad.append((d,k,m,res)); continue
            c=consts(u,Fu,S,W,V,g,eta,nu2,om,M,k)
            out[f"{d},{k},{m}"]=c
            for name,f in claims.items():
                if not f(c): viol[name].append((d,k,m,round(c[keyf[name]],5)))
    print('unsolved',len(bad),bad[:5])
    for name in claims:
        v=viol[name]
        print(name,'violations:',len(v),'examples:',v[:6], 'PASS' if not v else 'FAIL')
    json.dump(out,open(sys.argv[2],'w'))
