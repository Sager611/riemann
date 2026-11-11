
import numpy as np, mpmath as mp
d=np.load('/Users/adriansager/playground/riemann/.scratch/chaudit/zp.npy'); g=d[0]; zp=d[1]
print('zeros',len(g),'max gamma',g[-1])
# (1) Platt sums: delta*sum_{0<gamma<T} |coth(delta rho)|/|zeta'(rho)|, delta=pi/(2T)
def platt(T):
    dl=np.pi/(2*T); m=g<T; z=dl*(0.5+1j*g[m])
    return dl*np.sum(np.abs(1/np.tanh(z))/zp[m])
for T in [1e3,1e4,1e5,2e5,5e5,1e6+1]:
    print('PLATT T=%-10g  one-sided sum = %.12f'%(T,platt(T)))
print('paper (lem:platt109) T=1e6+1: 2.66161277991001')
# also plain sums
for T in [1e3,1e4,1e5,1e6]:
    m=g<=T; print('T=%g  sum 1/|rho zeta\'|=%.8f   sum 1/|zeta\'|=%.4f  N=%d'%(T,np.sum(1/(np.abs(0.5+1j*g[m])*zp[m])),np.sum(1/zp[m]),m.sum()))
# (2) vihuela sums
def vih(T,sigma,xi):
    m=g<=T; rho=0.5+1j*g[m]
    u=(rho-sigma)/(1j*T)          # theta_{T,sigma}=1-u
    om=(1-u)/np.tan(np.pi*u)      # -theta cot(pi theta) = (1-u)cot(pi u)
    y=np.pi*(1-sigma)/T
    c=(1/np.pi)*(y/np.tanh(y) if y!=0 else 1.0)
    om=om+c
    th1=1-(rho-1)/(1j*T)
    return (2*np.pi/T)*np.sum(np.abs(om+xi*th1*1j))
def B(T): L=np.log(T/(2*np.pi)); return L*L/(2*np.pi)-1.01*L/(6*np.pi)
def CT(T): L=np.log(T/(2*np.pi)); return L*L/(2*np.pi)-L/(6*np.pi)
for T in [1e3,1e4,1e5,1e6,1.13e6]:
    vals=[vih(T,s,x) for s in (0,1) for x in (1,-1)]
    print('VIH T=%-8g sums(s0x+,s0x-,s1x+,s1x-)='%T+' '.join('%.5f'%v for v in vals)+'  bound(1.01)=%.5f  C_T=%.5f  max-bound=%.5f'%(B(T),CT(T),max(vals)-B(T)))
m=g<=2e4; print('2*sum_{gamma<=2e4} 1/gamma = %.7f (paper 10.319317...)'%(2*np.sum(1/g[m])), ' (1/2pi)log^2(t0/2pi)-that = %.6f (paper 0.03435)'%(np.log(2e4/(2*np.pi))**2/(2*np.pi)-2*np.sum(1/g[m])))
# (3) constants C1, C2 (Lemma konstanz): C_k = sum_n zeta(2n)(1/(2n)^k - 2/(2n+1)^k + 1/(2n+2)^k)
mp.mp.dps=30
def Ck(k): return mp.nsum(lambda n: mp.zeta(2*n)*(1/(2*n)**k-2/(2*n+1)**k+1/(2*n+2)**k),[1,mp.inf])
C1=Ck(1); C2=Ck(2)
print('C1=%s C2=%s  2C1=%s  2(C1-C2)=%s'%(mp.nstr(C1,12),mp.nstr(C2,12),mp.nstr(2*C1,12),mp.nstr(2*(C1-C2),12)))
print('1.01/3=%.6f'%(1.01/3))
# C_T at the paper's T
T=3e12+1+np.pi/3; print('C_T(3e12+1+pi/3)= %.6f'%CT(T))
