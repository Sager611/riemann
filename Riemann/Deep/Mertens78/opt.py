# opt.py -- float exploration of parameters (not rigorous; rigorous check is cert.py)
import math, sys
from mpmath import zeta as mz
B = 0.875
def z(x): return float(mz(x))
def consts(s1, eta, masses, T1=1000.0):
    kap = (0.125 + eta)/2
    A = math.log(3*z(1+eta)) + kap*math.log(8/15)
    Lam = sum(l for _,_,l in masses); a = (Lam-1)*kap
    K = (Lam-1)*A; R = s1
    for sig,y,l in masses:
        c = abs(y)
        if l > 0: K += l*math.log(z(sig)/z(2*sig)); R += l*((1+kap)*(c+1+sig)+2*c+2)
        else: K += -l*math.log(z(sig)); R += -l*((3+2*kap)*c+sig)
    return a, K + R/T1
def best(s1, eta, masses, T1=1000.0):
    a, logC = consts(s1, eta, masses, T1)
    if a <= 0 or a >= 1: return None
    th = (s1 + a)/(1 + a)
    logD = logC + math.log(1/a + 1/(1-a)) + a*math.log(2) - math.log(math.pi)   # T4 limit = D c_h^{-a}
    logch = (math.log(4*a) + logD)/(1+a); logCC = logch - math.log(4) + math.log(1 + 1/a)
    return dict(a=a, logC=logC, theta=th, logch=logch, C=math.exp(min(logCC, 700)), log10C=logCC/math.log(10))
def X_of(theta_target, r):
    # |M(x)| <= C x^th <= x^target for x >= C^{1/(target-th)}
    if theta_target <= r['theta']: return float('inf')
    return r['log10C']/(theta_target - r['theta'])   # log10 X
def ramare_cross(r):
    # smallest log10 x with C x^th <= 0.013 x / log x
    lo = 1.0
    for L10 in [k*0.5 for k in range(2, 20000)]:
        L = L10*math.log(10)
        if r['log10C']*math.log(10) + r['theta']*L <= math.log(0.013) + L - math.log(L): return L10
    return None
if __name__ == "__main__":
    res = []
    for s1 in [0.90,0.905,0.91,0.915,0.92,0.925,0.93,0.94,0.95]:
        for s0 in [1.002,1.005,1.01,1.02,1.03,1.05,1.1]:
            for eta in [0.002,0.005,0.01,0.02,0.03,0.05,0.1]:
                rho = (s0-B)/(s1-B)
                r = best(s1, eta, [(s0,0.0,rho)])
                if r: res.append((s1,s0,eta,r))
    for tgt in [0.95,0.94,0.935,0.93,0.925]:
        bb = min(res, key=lambda q: X_of(tgt, q[3]))
        s1,s0,eta,r = bb
        print(f"target {tgt}: log10X={X_of(tgt,r):.1f} s1={s1} s0={s0} eta={eta} theta*={r['theta']:.5f} log10C={r['log10C']:.2f} cross={ramare_cross(r)}")
    bb = min(res, key=lambda q: q[3]['theta'])
    print("min theta*", bb[0],bb[1],bb[2], bb[3]['theta'], bb[3]['log10C'])
    # best Ramare crossover
    bc = min(res, key=lambda q: ramare_cross(q[3]) or 1e9)
    print("best crossover", bc[0],bc[1],bc[2], bc[3]['theta'], bc[3]['log10C'], ramare_cross(bc[3]))

