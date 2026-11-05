# T2: positivity vs amplification.  Zero rho = 1/2 + alpha + i x enters the explicit formula as H(x - i alpha).
# (a) line-positive kernel (Selberg/AM type): g = (1-|u|/R)_+, h(z) = R (sin(Rz/2)/(Rz/2))^2 >= 0 on R u iR.
# (b) full-strip-positive kernel (Poitou/Odlyzko type): G(u) = k(u)/cosh(u/2), k = (1-|u|/R)_+ positive definite.
#     Re H >= 0 on |Im z| <= 1/2; Harnack (strip -> disc via tanh(pi z/2)): Re H(x - i a) <= tan(pi/4 + pi a/2) Re H(x).
import numpy as np, math
def h_line(z,R):
    w=R*z/2
    return R*(np.sin(w)/w)**2
def H_strip(xs,a,R,du=2e-3):
    u=np.arange(-R,R+du/2,du); k=np.clip(1-np.abs(u)/R,0,None); G=k/np.cosh(u/2)
    # H(x - i a) = int G(u) e^{i(x - i a)u} du = int G(u) e^{a u} e^{i x u} du ; Re part = int G cosh(a u) cos(x u)
    W=G*np.cosh(a*u)*du
    return np.array([np.sum(W*np.cos(x*u)) for x in xs])
print("(a) line-positive kernel: amplification h(i a)/h(0) and worst paired weight min_x 2Re h(x - i a) / h(0)")
for R in (5,10,20,40):
    xs=np.linspace(1e-6,12,24001)
    for a in (0.1,0.25,0.3085):
        amp=float(np.real(h_line(np.array([1j*a]),R))[0])/R
        pair=2*np.real(h_line(xs-1j*a,R))
        print(f"  R={R:3d} a={a:.4f}: h(ia)/h(0)={amp:12.4g}  min_x pair/h(0)={pair.min()/R:12.4g} at x={xs[pair.argmin()]:.3f}")
print("(b) full-strip kernel: min Re H over strip grid, amplification Re H(-ia)/Re H(0), max_x ReH(x-ia)/ReH(x), Harnack cap tan(pi/4+pi a/2)")
for R in (5,10,20,40):
    xs=np.linspace(0,12,601)
    H0=H_strip(xs,0.0,R)
    mins=min(H_strip(xs,y,R).min() for y in (0.1,0.25,0.4,0.49,0.5))
    for a in (0.1,0.25,0.3085,0.45):
        Ha=H_strip(xs,a,R)
        ratio=np.max(Ha[H0>1e-9*H0.max()]/H0[H0>1e-9*H0.max()])
        cap=math.tan(math.pi/4+math.pi*a/2)
        print(f"  R={R:3d} a={a:.4f}: ReH(-ia)/ReH(0)={Ha[0]/H0[0]:8.4f}  max_x ratio={ratio:8.4f}  cap={cap:7.4f}   (min ReH on strip grid {mins:.3e})")
