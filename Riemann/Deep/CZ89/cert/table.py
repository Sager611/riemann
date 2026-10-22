
# Rigorous 1D tables for q = K_*^2 with python-flint (Arb ball arithmetic).
import sys, time, numpy as np
from flint import arb, ctx
ctx.prec = 96
NT = int(sys.argv[1]) if len(sys.argv) > 1 else 13   # table step h_t = 2^-NT
HT = 2.0**-NT
XMAXN = 9 * 2**NT          # grid covers [0, 9]
PI = arb.pi()
B = arb(1) / arb(2).sqrt()                 # b = 1/sqrt 2
C0 = B / (2 * B.sin())                     # b/(2 sin b)
AK = [arb(2)/25, arb(1)/25, -arb(1)/200]   # a_1,a_2,a_3
NS = 12
FACT = [arb(1)]
for i in range(1, 2*NS+3): FACT.append(FACT[-1]*i)
CN = [arb((-1)**n) / FACT[2*n+1] for n in range(NS+1)]
TAIL = arb(0, 1e-18)       # |tail| <= 2/(2NS+1)! < 1e-18 for |z|<=1
def sder(z):
    zm = abs(float(z.mid())); zr = float(z.rad())
    if zm > 0.5 + zr:
        sn = z.sin(); cs = z.cos()
        s0 = sn / z
        s1 = (cs - s0) / z
        s2 = -s0 - 2*s1/z
        return s0, s1, s2
    assert zm + zr <= 1.0
    z2 = z*z
    s0 = arb(0); s1 = arb(0); s2 = arb(0)
    for n in range(NS, -1, -1):
        s0 = s0*z2 + CN[n]
    for n in range(NS, 0, -1):
        s1 = s1*z2 + CN[n]*(2*n)
        s2 = s2*z2 + CN[n]*(2*n)*(2*n-1)
    s1 = s1*z
    return s0 + TAIL, s1 + TAIL, s2 + TAIL
def Kders(x):
    t = PI*x
    k0 = arb(0); k1 = arb(0); k2 = arb(0)
    for z, c in ((t - B, C0), (t + B, C0)):
        a0, a1, a2 = sder(z); k0 += c*a0; k1 += c*a1; k2 += c*a2
    for k in (1, 2, 3):
        c = AK[k-1] / 2
        for z in (PI*(x - k), PI*(x + k)):
            a0, a1, a2 = sder(z); k0 += c*a0; k1 += c*a1; k2 += c*a2
    return k0, PI*k1, PI*PI*k2
def lohi(v):
    m = float(v.mid()); r = float(v.rad())
    e = r + 1e-15*(1.0 + abs(m))
    return m - e, m + e
if __name__ == "__main__":
    t0 = time.time()
    N = XMAXN
    qlo = np.empty(N+1); qhi = np.empty(N+1); q1lo = np.empty(N+1); q1hi = np.empty(N+1)
    cq2 = np.empty(N); cqb = np.empty(N)
    maxrad = 0.0
    for m in range(N+1):
        x = arb(m) * arb(HT)          # exact dyadic
        k0, k1, k2 = Kders(x)
        q = k0*k0; q1 = 2*k0*k1
        maxrad = max(maxrad, float(q.rad()), float(q1.rad()))
        qlo[m], qhi[m] = lohi(q); q1lo[m], q1hi[m] = lohi(q1)
        if m < N:
            X = arb(arb(2*m+1) * arb(HT/2), HT/2)   # ball = cell [m h, (m+1) h]
            k0, k1, k2 = Kders(X)
            q2 = 2*(k1*k1 + k0*k2)
            lo, hi = lohi(q2)
            cq2[m] = max(abs(lo), abs(hi))
            Q = k0*k0
            cqb[m] = max(0.0, lohi(Q)[0])        # direct ball lower bound of q on the cell
        if m % 8192 == 0:
            print("m", m, "t", round(time.time()-t0,1), flush=True)
    qlo = np.maximum(qlo, 0.0)                   # q is a square
    cqmin = np.maximum(np.maximum(np.minimum(qlo[:-1], qlo[1:]) - cq2*HT*HT/8.0, cqb), 0.0)
    np.savez("tables.npz" if NT == 13 else "tables_%d.npz" % NT, qlo=qlo, qhi=qhi, q1lo=q1lo, q1hi=q1hi, cq2=cq2, cqmin=cqmin, HT=HT)
    print("done", N, "maxrad", maxrad, "max cq2", cq2.max(), "time", round(time.time()-t0,1))

