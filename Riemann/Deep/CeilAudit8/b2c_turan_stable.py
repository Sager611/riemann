# ceilaudit8 B2c: stable L_1(f^(m)) (Turan-type = Fleischer on |f^(m)|) for the zeta off-line pair model.
# f = P_near(x) * pair(x) * g_far(x); Taylor in u = h/gbar: T = prod((x-r)/gbar + u), pair = (w^2+d^2) + 2w u + u^2 (gbar units),
# g_far/g_far(x) = sum_i B_i(D^u)/i! u^i with D^u_k = gbar^k (log g_far)^(k) (far zeros + mirrors + smooth tail; no cancellation).
import numpy as np, math, sys
Z = np.load('Riemann/Deep/ExploreHelio/results/zeta_zeros.npy'); T1 = Z[-1]; M = 9; RN = 9.0
def far_D(x, far_real, far_cpl, gbar):
    S = np.zeros((M+1, x.size))
    for j in range(0, far_real.size, 1000):   # sized: 1201 x 1000 doubles = 9.6 MB per temporary
        r = 1.0/(x[:,None]-far_real[None, j:j+1000]); p = r.copy()
        for k in range(1, M+1): S[k] += p.sum(1); p *= r
    for t in far_cpl:
        r = 1.0/(x-t); p = r.copy()
        for k in range(1, M+1): S[k] += p.real; p = p*r
    g = np.geomspace(T1, T1*1e4, 4001); w = np.gradient(g)*np.log(g/(2*math.pi))/(2*math.pi)
    for j in range(0, x.size, 200):
        xa = x[j:j+200, None]
        for k in range(1, M+1): S[k][j:j+200] += (((xa-g[None,:])**(-k) + (xa+g[None,:])**(-k))*w[None,:]).sum(1)
    D = np.zeros((M+1, x.size))
    for k in range(1, M+1): D[k] = ((-1)**(k-1))*math.factorial(k-1)*S[k]*gbar**k
    B = [np.ones(x.size)]
    for n in range(0, M): B.append(sum(math.comb(n, i)*B[n-i]*D[i+1] for i in range(0, n+1)))
    return np.array([B[i]/math.factorial(i) for i in range(M+1)])
def conv(a, b):
    out = np.zeros((M+1, a.shape[1]))
    for i in range(min(b.shape[0], M+1)): out[i:] += b[i]*a[:M+1-i]
    return out
def L1m(x, real, pair, gbar, x0):
    near = real[np.abs(real-x0) <= RN]; far = real[np.abs(real-x0) > RN]
    T = np.zeros((M+1, x.size)); T[0] = 1.0
    for r in near:
        a = (x-r)/gbar; T[1:] = a*T[1:] + T[:-1]; T[0] = a*T[0]
    cpl = []
    if pair is not None:
        g0, d = pair; w = (x-g0)/gbar; dd = d/gbar
        P = np.zeros((3, x.size)); P[0] = w*w+dd*dd; P[1] = 2*w; P[2] = 1.0; T = conv(T, P); cpl = [-g0+1j*d, -g0-1j*d]
    F = conv(T, far_D(x, far, cpl, gbar))
    out = []
    for m in range(0, M-1):
        a = math.factorial(m+1)**2*F[m+1]**2; b = math.factorial(m)*math.factorial(m+2)*F[m]*F[m+2]
        out.append(((a-b)/(a+np.abs(b))).min())
    return out
def run(ia, ib, dgs, tag):
    ga, gb = Z[ia], Z[ib]; g0 = (ga+gb)/2; gbar = 2*math.pi/math.log(g0/(2*math.pi))
    keep = np.ones(Z.size, bool); keep[[ia, ib]] = False; zz = Z[keep]; real = np.concatenate([zz, -zz])
    x = np.linspace(g0-6, g0+6, 1201)+1.3e-6
    print(f"== {tag} g0={g0:.4f} gbar={gbar:.4f}; min over |x-g0|<=6 of normalised L_1(f^(m)), m=0..{M-2}")
    q = L1m(x, np.concatenate([Z, -Z]), None, gbar, g0)
    print("control (true zeros): " + " ".join(f"{v:+.3f}" for v in q)); sys.stdout.flush()
    for dg in dgs:
        q = L1m(x, real, (g0, dg*gbar), gbar, g0)
        print(f"delta/gbar={dg:4.2f}: " + " ".join(f"{v:+.3f}" for v in q) + "  fail m=" + str([m for m, v in enumerate(q) if v < -1e-9])); sys.stdout.flush()
run(6708, 6709, [0.5, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0], 'Lehmer')
run(6606, 6607, [1.0, 1.25, 1.5, 2.0, 2.5, 3.0], 'generic')
