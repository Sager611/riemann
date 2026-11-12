# gl2audit: independent float ceiling (own minimiser; mpmath quadrature, Chebyshev basis; no gl2ceil code)
import sys, os, mpmath as mp
mp.mp.dps = int(os.environ.get('DPS', '30'))
XI = mp.mpf(sys.argv[1])/mp.mpf(sys.argv[2]) if len(sys.argv) > 2 else mp.mpf(25)/256
K = int(sys.argv[3]) if len(sys.argv) > 3 else 16
N = 90
# Gauss-Legendre nodes on [0,1], Newton-refined in mpmath
def gl(n):
    xs, ws = [], []
    for i in range(1, n+1):
        x = mp.cos(mp.pi*(i-mp.mpf(1)/4)/(n+mp.mpf(1)/2))
        for _ in range(60):
            p0, p1 = mp.mpf(1), x
            for k in range(2, n+1):
                p0, p1 = p1, ((2*k-1)*x*p1-(k-1)*p0)/k
            dp = n*(x*p1-p0)/(x*x-1)
            dx = p1/dp; x -= dx
            if abs(dx) < mp.mpf(10)**(2-mp.mp.dps): break
        xs.append((1-x)/2); ws.append(1/((1-x*x)*dp*dp))  # weight 2/((1-x^2)P'^2) /2
    return xs, ws
YS, WS = gl(N)
def cheb(n, t):  # T_n(t), T_n'(t)
    T0, T1 = mp.mpf(1), t; d0, d1 = mp.mpf(0), mp.mpf(1)
    if n == 0: return T0, d0
    for k in range(2, n+1):
        T0, T1 = T1, 2*t*T1-T0; d0, d1 = d1, 2*T0+2*t*d1-d0
    return T1, d1
# basis: phi_0 = 1, phi_k = T_{2k-1}(1-2y) - 1 (real class); values and y-derivatives at nodes
BV = [[(mp.mpf(1), mp.mpf(0))] + [(lambda T: (T[0]-1, -2*T[1]))(cheb(2*k-1, 1-2*y)) for k in range(1, K+1)] for y in YS]
def mats(R):
    n = K+1; G0 = mp.zeros(n, n); G1 = mp.zeros(n, n)
    for y, wt, bv in zip(YS, WS, BV):
        e = wt*mp.exp(2*R*y); f = [v for v, d in bv]; g = [R*v+d for v, d in bv]
        for i in range(n):
            for j in range(i, n):
                G0[i, j] += e*f[i]*f[j]; G1[i, j] += e*g[i]*g[j]
    for i in range(n):
        for j in range(i): G0[i, j] = G0[j, i]; G1[i, j] = G1[j, i]
    return G0, G1
def CB(lam):  # by quadrature of P = sinh(lam x)/sinh lam (independent of closed forms)
    if lam == 0: return mp.mpf(1), mp.mpf(1)/3
    s = mp.sinh(lam); C = sum(wt*(lam*mp.cosh(lam*y)/s)**2 for y, wt in zip(YS, WS)); B = sum(wt*(mp.sinh(lam*y)/s)**2 for y, wt in zip(YS, WS))
    return C, B
def solveQ(R, G0, G1, al, ga):
    n = K+1; E = mp.exp(2*R)
    H = mp.matrix(K, K); g = mp.matrix(K, 1)
    for i in range(1, n):
        g[i-1] = al*G0[i, 0]+ga*G1[i, 0]-E
        for j in range(1, n): H[i-1, j-1] = al*G0[i, j]+ga*G1[i, j]+2*E
    b = -mp.lu_solve(H, g)
    c0 = (1+E)/2+al*G0[0, 0]+ga*G1[0, 0]
    J = c0+sum(g[i]*b[i] for i in range(K))
    bb = [mp.mpf(1)]+[b[i] for i in range(K)]
    W0 = sum(bb[i]*bb[j]*G0[i, j] for i in range(n) for j in range(n)); W1 = sum(bb[i]*bb[j]*G1[i, j] for i in range(n) for j in range(n))
    return J, W0, W1, bb
def Mlam(R, lam, G=None):
    G0, G1 = G or mats(R); C, B = CB(lam); return solveQ(R, G0, G1, C/XI, XI*B)[0]
def Malt(R, lam0=mp.mpf('0.77'), it=60):  # alternating P/Q descent; returns min_Q min_P c and lambda
    G0, G1 = mats(R); lam = lam0
    for _ in range(it):
        C, B = CB(lam); J, W0, W1, bb = solveQ(R, G0, G1, C/XI, XI*B)
        nl = XI*mp.sqrt(W1/W0)
        if abs(nl-lam) < mp.mpf(10)**-14: lam = nl; break
        lam = nl
    E = mp.exp(2*R); w1sq = E*(1-2*sum(bb[1:]))**2
    F = (1+w1sq)/2+mp.sqrt(W0*W1)/mp.tanh(XI*mp.sqrt(W1/W0))
    return F, lam, bb
def kap(R): F, lam, bb = Malt(R); return 1-mp.log(F)/R, lam
if __name__ == '__main__':
    R0 = mp.mpf(sys.argv[4]) if len(sys.argv) > 4 else mp.mpf('9.867773')
    k, lam = kap(R0); print('K', K, 'R', R0, 'kappa', mp.nstr(k, 13), 'lam', mp.nstr(lam, 9))
    a, b = R0-mp.mpf('0.6'), R0+mp.mpf('0.6'); gr = (mp.sqrt(5)-1)/2
    c, d = b-gr*(b-a), a+gr*(b-a); fc, fd = kap(c)[0], kap(d)[0]
    for _ in range(int(sys.argv[5]) if len(sys.argv) > 5 else 22):
        if fc > fd: b, d, fd = d, c, fc; c = b-gr*(b-a); fc = kap(c)[0]
        else: a, c, fc = c, d, fd; d = a+gr*(b-a); fd = kap(d)[0]
    Rm = (a+b)/2; k, lam = kap(Rm); print('max over R: R', mp.nstr(Rm, 9), 'kappa', mp.nstr(k, 13), 'lam', mp.nstr(lam, 9), flush=True)
