import numpy as np, math
def H(xs,a,R,du=2e-3):
    u=np.arange(-R,R+du/2,du); k=np.clip(1-np.abs(u)/R,0,None); G=k/np.cosh(u/2); W=G*np.cosh(a*u)*du
    return np.array([np.sum(W*np.cos(x*u)) for x in xs])
xs=np.linspace(0,3,301)
for R in (5,10,40,160):
    h=H(xs,0,R); half=xs[np.argmax(h<h[0]/2)]
    # Markov lower bound for strip kernels: Z ~ k(0) log q, pair weight <= 2 C(a) H(0)
    print(f"R={R:4d}: H(0)={h[0]:.4f} (2pi={2*math.pi:.4f})  half-width in x = {half:.3f}  (sech(pi x) half-width {math.acosh(2)/math.pi:.3f})")
for q in (1e4,1e6,1e12,1e30):
    for a in (0.1,0.25):
        C=math.tan(math.pi/4+math.pi*a/2)
        print(f"q={q:.0e} a={a}: Markov bound >= log q/(4 pi C) = {math.log(q)/(4*math.pi*C):.3f}")
