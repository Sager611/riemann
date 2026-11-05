# T2b: Davenport-Heilbronn off-line zero rho = 0.808517182457 + 85.699348485378 i (alpha = 0.308517182457) in both kernels.
import numpy as np, math
a=0.308517182457; x0=85.699348485378
def h_line(z,R):
    w=R*z/2; return R*(np.sin(w)/w)**2
def ReH_strip(x,y,R,du=1e-3):
    u=np.arange(-R,R+du/2,du); k=np.clip(1-np.abs(u)/R,0,None); G=k/np.cosh(u/2)
    return float(np.sum(G*np.cosh(y*u)*np.cos(x*u))*du)
print("R   line: 2Re h(x0 - i a)/h(0)   line: sup_{|x-x0|<1} |2Re h(x - i a)|/h(0)   strip: ReH(x0 - i a)/ReH(x0) [cap 3.2235]")
for R in (5,10,20,40,80):
    v=2*np.real(h_line(np.array([x0-1j*a]),R))[0]/R
    xs=np.linspace(x0-1,x0+1,4001); s=np.abs(2*np.real(h_line(xs-1j*a,R))).max()/R
    rs=ReH_strip(x0,a,R)/ReH_strip(x0,0.0,R)
    print(f"{R:3d}  {v: .4e}   {s:.4e}   {rs:.4f}")
