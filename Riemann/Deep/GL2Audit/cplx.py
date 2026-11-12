# gl2audit: does jointly complex (P, Q) beat the real ceiling? (sesquilinear extension; NOT in Bernard's stated class)
# local test at the real optimum (Hessian in the imaginary directions) + line search + direct double-quadrature check
import sys
args = sys.argv[1:]; sys.argv = ['myceil', args[0], args[1], args[2]]
from myceil import *
R = mp.mpf(args[3]); KE = int(args[4]) if len(args) > 4 else K
F, lam, bb = Malt(R)
print('real: R', mp.nstr(R, 9), 'kappa', mp.nstr(1-mp.log(F)/R, 13), 'lam', mp.nstr(lam, 9))
E = mp.exp(2*R)
def Qr(y):
    t = 1-2*y; v = bb[0]; d = mp.mpf(0)
    for k in range(1, K+1):
        T, dT = cheb(2*k-1, t); v += bb[k]*(T-1); d += -2*bb[k]*dT
    return v, d
def chi(k, y):
    T, dT = cheb(2*k, 1-2*y); return T-1, -2*dT
# Gram data of the even (imaginary) basis and the coupling functional l_k = int e^{2Ry}(chi_k Qr' - Qr chi_k')
G0i = mp.zeros(KE, KE); G1i = mp.zeros(KE, KE); l = mp.zeros(KE, 1); W0r = mp.mpf(0); W1r = mp.mpf(0)
for y, wt in zip(YS, WS):
    e = wt*mp.exp(2*R*y); q, dq = Qr(y); W0r += e*q*q; W1r += e*(R*q+dq)**2
    cv = [chi(k, y) for k in range(1, KE+1)]
    for i in range(KE):
        l[i] += e*(cv[i][0]*dq-q*cv[i][1])
        for j in range(KE):
            G0i[i, j] += e*cv[i][0]*cv[j][0]; G1i[i, j] += e*(R*cv[i][0]+cv[i][1])*(R*cv[j][0]+cv[j][1])
Xr = mp.sqrt(W0r*W1r); z = XI*Xr/W0r
fX = mp.coth(z)-z/mp.sinh(z)**2; fW = (XI*Xr**2/W0r**2)/mp.sinh(z)**2
S = (fX/(2*Xr))*(W0r*G1i+W1r*G0i)+fW*G0i
Sl = mp.lu_solve(S, l); test = (fX/(2*Xr))*sum(l[i]*Sl[i] for i in range(KE))
print('KE', KE, 'test (>1 means complex helps locally):', mp.nstr(test, 12))
# second-order coefficient along the optimal direction e = S^-1 l:  q(e) = e^T S e - (fX/2Xr)(l.e)^2
qe = sum(Sl[i]*l[i] for i in range(KE)) - (fX/(2*Xr))*sum(l[i]*Sl[i] for i in range(KE))**2
print('q(e*) =', mp.nstr(qe, 6), '(negative = descent)')
def Fc(eps):  # exact F_c along d = eps*e*
    A = W0r+eps**2*sum(Sl[i]*Sl[j]*G0i[i, j] for i in range(KE) for j in range(KE))
    C = W1r+eps**2*sum(Sl[i]*Sl[j]*G1i[i, j] for i in range(KE) for j in range(KE))
    I = eps*sum(l[i]*Sl[i] for i in range(KE)); X = mp.sqrt(A*C-I*I)
    return (1+E*(bb[0]-2*sum(bb[1:]))**2)/2+X/mp.tanh(XI*X/A), A, C, I
for eps in ['1e-6', '1e-4', '1e-3', '1e-2', '0.1', '1']:
    f = Fc(mp.mpf(eps))[0]; print('eps', eps, 'kappa_c', mp.nstr(1-mp.log(f)/R, 13))
