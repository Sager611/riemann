# T2: EP normal form of de Bruijn-Newman collisions (Lehmer pair 7005) + moved-pair adversarial
import numpy as np, sys, math
d = np.loadtxt('zeros_6670_6750.txt')
n, g = d[:,0].astype(int), d[:,1]
z0 = 2.0*g                      # zeros of H_0(z) = xi(1/2+iz/2)/8
ip = int(np.where(n==6709)[0][0])   # Lehmer pair = (6709, 6710)
a_, b_ = z0[0]-0.5*(z0[1]-z0[0]), z0[-1]+0.5*(z0[-1]-z0[-2])
D = lambda x: math.log(x/(4*math.pi))/(4*math.pi)   # zero density in z
Dl, Dr, L = D(a_), D(b_), 400.0
def Eext(x):  # static field 2*int D/(x-y) dy from [a-L,a] and [b,b+L]
    return 2*(Dl*np.log(np.abs((x-(a_-L))/(x-a_))) + Dr*np.log(np.abs((x-b_)/(x-(b_+L)))))
def Eext_d(x): # derivative of Eext
    return 2*(Dl*(1/(x-(a_-L)) - 1/(x-a_)) + Dr*(1/(x-b_) - 1/(x-(b_+L))))
def F(z):
    dz = z[:,None]-z[None,:]
    np.fill_diagonal(dz, np.inf)
    return 2*np.sum(1.0/dz, axis=1) + Eext(z)
def S_of(z, i, j):  # bath curvature at the pair midpoint, excluding the pair
    xm = 0.5*(z[i]+z[j]); m = np.ones(len(z),bool); m[[i,j]] = False
    return np.sum(1.0/(xm-z[m])**2) - 0.5*Eext_d(xm)
def rk4(z, h):
    k1=F(z); k2=F(z+0.5*h*k1); k3=F(z+0.5*h*k2); k4=F(z+h*k3)
    return z + h*(k1+2*k2+2*k3+k4)/6
# ---- backward flow to the EP
z = z0.copy(); t = 0.0; rec = []
d0 = z[ip+1]-z[ip]; S0 = S_of(z, ip, ip+1)
while True:
    dl = z[ip+1]-z[ip]
    rec.append((t, dl*dl))
    if dl < 2e-4: break
    h = -min(2e-6, 0.02*dl*dl/8)
    z = rk4(z, h); t += h
gaps = np.diff(z); 
tc = t - (z[ip+1]-z[ip])**2/8
rec = np.array(rec)
# slope of delta^2 vs t over the last stretch and at t=0
sl_end = np.polyfit(rec[-200:,0], rec[-200:,1], 1)[0]
sl_0 = (rec[5,1]-rec[0,1])/(rec[5,0]-rec[0,0])
tc_pb = -(1/(4*S0))*math.log(1/(1-S0*d0*d0/2))
print('Lehmer pair z = %.6f, %.6f; gap d0 = %.6f; mean gap %.4f; S0 = %.4f' % (z0[ip], z0[ip+1], d0, 4*math.pi/math.log(z0[ip]/(4*math.pi)), S0))
print('EP (full ODE, %d zeros): t_c = %.6e ; -d0^2/8 = %.6e ; pair+bath closed form %.6e' % (len(z0), tc, -d0*d0/8, tc_pb))
print('slope d(delta^2)/dt: at t=0 %.4f (pred 8-4*S0*d0^2 = %.4f), near EP %.5f (EP normal form: 8)' % (sl_0, 8-4*S0*d0*d0, sl_end))
print('min other gap at t_c: %.4f (others still real, so this is the first EP in the window)' % np.min(np.delete(gaps, ip)))
# ---- EP times of all adjacent pairs (pair+bath), real or complex
out = []
for i in range(5, len(z0)-6):
    dd = z0[i+1]-z0[i]; S = S_of(z0, i, i+1); x = S*dd*dd/2
    if x < 1: out.append((-(1/(4*S))*math.log(1/(1-x)), 0.0, n[i]))
    else: out.append(((1/(4*S))*math.log(x-1), math.pi/(4*S), n[i]))
out.sort(key=lambda r: -r[0] if r[1]==0 else 1e9)
real = [r for r in out if r[1]==0]; cplx = [r for r in out if r[1]>0]
print('pairs %d: real EP %d, complex EP %d; nearest real EPs:' % (len(out), len(real), len(cplx)), ' '.join('%d:%.2e' % (r[2], r[0]) for r in real[:4]))
if cplx: print('complex EPs (Re t, |Im t|) e.g.', ' '.join('%d:(%.3f,%.3f)' % (r[2], r[0], r[1]) for r in cplx[:3]))
# ---- adversarial: move the Lehmer pair off the line, forward flow to re-entry (EP at t>0)
for eta in (0.05, 0.2, 0.6):
    zc = z0.astype(complex).copy(); xm = 0.5*(z0[ip]+z0[ip+1])
    zc[ip], zc[ip+1] = xm+1j*eta, xm-1j*eta
    Sb = S_of(z0, ip, ip+1)
    t = 0.0
    while abs(zc[ip].imag) > 2e-4:
        y = abs(zc[ip].imag); h = min(2e-5, 0.02*y*y/2)
        zc = rk4(zc, h); t += h
    te = t + zc[ip].imag**2/2
    print('moved pair eta=%.2f: forward EP at t_e = %.6f ; de Bruijn eta^2/2 = %.6f ; bath formula %.6f' % (eta, te, eta*eta/2, (1/(4*Sb))*math.log(1+2*Sb*eta*eta)))
