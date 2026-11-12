
import numpy as np, scipy.linalg as sl, sys, time
np.set_printoptions(precision=6, suppress=True)
t0=time.time()
n=int(sys.argv[1]) if len(sys.argv)>1 else 600
x=(np.arange(n)+0.5)/n; w=1.0/n
def Aop(phi, theta=1.0):
    # A on [0,theta]: (Af)(x)= i phi f + (2i/pi) int_0^x sin(pi phi u)/u f(x-u) du ; trapezoid-type Volterra
    xs=theta*x; ws=theta*w
    U=xs[:,None]-xs[None,:]
    K=np.where(U>0, np.sin(np.pi*phi*np.where(U>0,U,1.0))/np.where(U>0,U,1.0), 0.0)
    K=K+np.diag(np.full(n,0.5*np.pi*phi))*0.0
    # diagonal cell: integral over u in [0, ws/2] of kernel ~ pi phi * ws/2
    M=(2j/np.pi)*(K*ws + np.diag(np.full(n,np.pi*phi*ws/2)))
    return 1j*phi*np.eye(n)+M, xs, ws
def ip(f,g,ws): return np.sum(f*np.conj(g))*ws
for phi in [0.3,0.5,0.6,0.75,0.9]:
    a=np.pi*phi
    A,xs,ws=Aop(phi)
    C=sl.expm(np.pi*A)
    # closed form C f = e^{ia}(f + 2ia int_0^x e^{ia(x-y)} f(y) dy)
    U=xs[:,None]-xs[None,:]
    Cc=np.exp(1j*a)*(np.eye(n)+2j*a*np.where(U>0,np.exp(1j*a*U),0.0)*ws+np.diag(np.full(n,2j*a*ws/2)))
    g=np.exp(-1j*a*xs)-np.exp(1j*a)*np.exp(1j*a*xs)
    gn=g/np.sqrt(ip(g,g,ws).real)
    # moments
    mom=[ip(np.linalg.matrix_power(A,j)@gn,gn,ws) for j in range(8)]
    # density p0,p1 from resolvent of C (numerical expm)
    X=np.linspace(-12,12,4801); dx=X[1]-X[0]
    I=np.eye(n)
    p0=np.array([np.exp(np.pi*s)*ip(np.linalg.solve(np.exp(np.pi*s)*I-C,gn),gn,ws).imag for s in X[::8]])
    p1=np.array([np.exp(np.pi*s)*ip(np.linalg.solve(-np.exp(np.pi*s)*I-C,gn),gn,ws).imag for s in X[::8]])
    Xs=X[::8]; d=Xs[1]-Xs[0]
    mnu=[np.sum(Xs**j*p0)*d+np.sum((Xs+1j)**j*p1)*d for j in range(8)]
    # p0(0) closed form
    p00=np.sin(a)*np.sin(2*a)/(a*ip(g,g,ws).real)
    i0=np.argmin(abs(Xs))
    # Im C min eigenvalue
    ImC=(C-C.conj().T)/(2j); ev=np.linalg.eigvalsh((ImC+ImC.conj().T)/2)
    print(f"phi={phi}: |expm-closed|={np.abs(C-Cc).max():.2e}  p0(0) num={p0[i0]:.6f} closed={p00:.6f}  min p0={p0.min():.4f} min p1={p1.min():.4f}  mass0={np.sum(p0)*d:.5f} mass1={np.sum(p1)*d:.5f}  lam_min(ImC)={ev.min():.5f}")
    print("   moments <A^j g,g>:", np.round(np.array(mom[:6]),5))
    print("   moments of nu    :", np.round(np.array(mnu[:6]),5))
print("time",time.time()-t0)
