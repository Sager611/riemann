# FodgeAudit bookkeeping.py: Prop. C.1 block by block at many t >= T1 (exact integers N_k, L_k; mpmath dps 30):
# (i) per-block claims of the proof: u_k in [1/r, u+], N_k > 1e4, L_k <= N_k/75, x_k >= y0, mu'_k <= mu'+, alpha_k <= abar(x_k), beta_k <= bbar,
#     E_k <= (e1 x_k + e0) t^(-1/3), u/sqrt(1+u) <= kappa log(1+u), u <= kappa'' log(1+u);
# (ii) total = head + sum_k 2 N_k^(-1/2) sqrt((L_k/x_k + c)(alpha_k L_k + beta_k x_k^2)) + RS(t)  <=  b1 t^(1/6) log t + b2 t^(1/6) + S;
# (iii) Lemma C.4 on exact block sums: max_{L'} |sum_{n<L'} e(f(n))|^2 / RHS for random t in [T1, 3 T1] (first, middle, last block).
import sys, math, random
import mpmath as mp
mp.mp.dps = 30
c, r, nu = mp.mpf('2.44'), 80, mp.mpf('0.95'); T1 = mp.mpf('2.9e12')
b1, b2, S = mp.mpf('0.3080532'), mp.mpf('2.405897'), mp.mpf('506.6')
# constants of the proof (point values; consts.py encloses them)
up = 1/mp.mpf(r) + 1/(nu*T1**(mp.mpf(1)/3)); N0p = nu*T1**(mp.mpf(1)/3) - mp.mpf(1)/2
phi = lambda v: 2*(1+v)/(1+mp.sqrt(1+v))
mupp = (1 + 1/(2*N0p))**1.5*phi(up*(N0p+mp.mpf(1)/2)/N0p); A = 8*mupp/15; ainf = 1/c + A*mp.sqrt(c)
abar = lambda y: 1/c + A*mp.sqrt(c + 1/y)
y0 = mp.cbrt(mp.pi)*nu*(1 + mp.mpf(1)/r)
K, C0, E1 = mp.mpf('4.524'), mp.mpf('1.0004'), mp.mpf('2.24')
bbar = 8*K/(3*mp.sqrt(c)) + C0/y0 + E1*(c + 1/y0)/(3*y0**3)
e1 = bbar*mp.cbrt(mp.pi)*(1+up)*r + c*bbar*mp.cbrt(mp.pi)**2*(1+up)**2*r*r/mp.cbrt(T1)
e0 = c*abar(y0)*mp.cbrt(mp.pi)*(1+up)*r
kap = (2+up)/(2*mp.sqrt(1+up)); kpp = up/mp.log(1+up)
def blocks(t):
    n1 = int(mp.floor(mp.sqrt(t/(2*mp.pi)))); N = int(mp.ceil(nu*mp.cbrt(t))); out = []
    while N <= n1:
        L = -(-N//r); out.append((N, L)); N += L
    return n1, out
def run_t(t):
    n1, bl = blocks(t); N0 = bl[0][0]; worst = {}
    def w(k, v):
        if v > worst.get(k, -1e99): worst[k] = v
    head = 2*mp.fsum(mp.mpf(n)**-0.5 for n in range(1, N0)) if N0 < 2*10**5 else 4*mp.sqrt(N0-1) - 2
    w('head/(4 sqrt(nu) t^1/6 - 2)', head/(4*mp.sqrt(nu)*t**(mp.mpf(1)/6) - 2))
    tot = head + 1.48*t**-0.25 + 0.127*t**-0.75; sumint = 0
    for (N, L) in bl:
        u = mp.mpf(L)/N; x = mp.cbrt(mp.pi/t)*(N+L); W = x**3; Np = N - mp.mpf(1)/2
        mu1 = 2*mp.mpf(N+L)**1.5/L*(1/mp.sqrt(Np) - 1/mp.sqrt(Np+L))
        al = 1/c + 8*mu1/15*mp.sqrt(c + 1/x); be = 8*K/(3*mp.sqrt(c)) + C0/x + E1*(c + 1/x)/(3*W)
        Ek = be*x*x/L + c*al*x/L + c*be*x**3/L**2
        w('1/(u r)', 1/(u*r)); w('u/u+', u/up); w('1e4/N', mp.mpf(10**4)/N); w('75 L/N', 75*u); w('y0/x', y0/x)
        w('mu/mu+', mu1/mupp); w('alpha/abar', al/abar(x)); w('beta/bbar', be/bbar); w('E/((e1x+e0)t^-1/3)', Ek/((e1*x + e0)*t**(-mp.mpf(1)/3)))
        w('u/sqrt(1+u)/(kap log(1+u))', u/mp.sqrt(1+u)/(kap*mp.log(1+u))); w('u/(kpp log(1+u))', u/(kpp*mp.log(1+u)))
        tot += 2/mp.sqrt(N)*mp.sqrt((L/x + c)*(al*L + be*x*x))
    Nend = bl[-1][0] + bl[-1][1]; w('Nend/((1+u+) n1)', mp.mpf(Nend)/((1+up)*n1))
    bound = b1*t**(mp.mpf(1)/6)*mp.log(t) + b2*t**(mp.mpf(1)/6) + S
    return len(bl), tot, bound, worst
if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'blocks'
    if mode == 'blocks':
        agg = {}
        ts = [T1, T1*1.0000001, mp.mpf('3e12'), mp.mpf('1e13'), mp.mpf('1e14'), mp.mpf('1e15'), mp.mpf('1e16'), mp.mpf('1e18'), mp.mpf('1e20'),
              mp.mpf('1e22'), mp.mpf('1e25'), mp.mpf('1e30'), mp.mpf('1e35'), mp.mpf('1e40')]
        random.seed(7); ts += [T1*mp.mpf(10)**(random.random()*8) for _ in range(30)]
        for t in ts:
            nb, tot, bound, worst = run_t(t)
            for k, v in worst.items(): agg[k] = max(agg.get(k, -1e99), v)
            print(f't={mp.nstr(t,8)} blocks={nb} total={mp.nstr(tot,10)} bound={mp.nstr(bound,10)} ratio={mp.nstr(tot/bound,6)}', flush=True)
        print('max over all t of each per-block ratio (each must be <= 1):')
        for k, v in agg.items(): print(f'  {k}: {mp.nstr(v,8)}')

    if mode == 'c4':
        import numpy as np
        random.seed(11); mx = (0, None); cnt = 0
        for it in range(int(sys.argv[2]) if len(sys.argv) > 2 else 100):
            t = T1*(1 + 2*mp.mpf(random.random()))
            n1, bl = blocks(t); tf = float(t/(2*mp.pi))
            for (N, L) in bl:
                x = mp.cbrt(mp.pi/t)*(N+L); W = x**3; Np = N - mp.mpf(1)/2
                mu1 = 2*mp.mpf(N+L)**1.5/L*(1/mp.sqrt(Np) - 1/mp.sqrt(Np+L))
                al = 1/c + 8*mu1/15*mp.sqrt(c + 1/x); be = 8*K/(3*mp.sqrt(c)) + C0/x + E1*(c + 1/x)/(3*W)
                rhs = float((L/x + c)*(al*L + be*x*x))
                p0 = (t/(2*mp.pi))*mp.log(N); p0 = float(p0 - mp.floor(p0))
                n = np.arange(L, dtype=float); d = tf*np.log1p(n/N); ph = p0 + (d - np.floor(d))
                ps = np.cumsum(np.exp(2j*np.pi*ph)); v = float(np.max(np.abs(ps)**2))/rhs; cnt += 1
                if v > mx[0]: mx = (v, (mp.nstr(t, 15), N, L))
        print(f'Lemma C.4 on exact block sums: {cnt} blocks (all blocks of 1..n1 cut, every L\' <= L), max |sum|^2/RHS = {mx[0]:.6f} at (t, N, L) = {mx[1]}')
