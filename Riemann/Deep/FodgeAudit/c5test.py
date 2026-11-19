# FodgeAudit c5test.py: Lemma C.5 as used in Lemmas C.6 / C.8 on exact sums of second / third differences of f(x) = (t/2pi) log(N+x).
#  2nd: S_{m,k} = sum_{n=0}^{L'-1-m-k} e(g~(n)), g~(x) = f(x+m+k) - f(x+k) - f(x+m) + f(x); bound mu'' L (mk/V)^(1/2) + 4.527 (V/mk)^(1/2) + 1.0006 + 2.272 mk/V,
#       V = pi (N+L)^4/(3t), mu'' = (N+L)^2/(N'(N'+L)), N >= 5e4, L <= N/75.
#  3rd: S_{m,k,j} with g^(x) = g~(x+j) - g~(x); bound mu3 L (mkj/V3)^(1/2) + 4.531 (V3/mkj)^(1/2) + 1.001 + 2.31 mkj/V3, V3 = pi (N+L)^5/(12t), N >= 2e4.
import sys, random, math
import mpmath as mp
mp.mp.dps = 45
def run(order, N, L, Lp, sh, lam):
    Nm = mp.mpf(N); Np = Nm - mp.mpf(1)/2; P = 1
    for x in sh: P *= x
    if order == 2:
        V = lam and None
        t = lam*mp.pi*(Nm+L)**4/(3*P); V = mp.pi*(Nm+L)**4/(3*t); mu = (Nm+L)**2/(Np*(Np+L))
        rhs = mu*L*mp.sqrt(P/V) + mp.mpf('4.527')*mp.sqrt(V/P) + mp.mpf('1.0006') + mp.mpf('2.272')*P/V
    else:
        t = lam*mp.pi*(Nm+L)**5/(12*P); V = mp.pi*(Nm+L)**5/(12*t); mu = 2*(Nm+L)**2.5/(3*L)*(Np**-1.5 - (Np+L)**-1.5)
        rhs = mu*L*mp.sqrt(P/V) + mp.mpf('4.531')*mp.sqrt(V/P) + mp.mpf('1.001') + mp.mpf('2.31')*P/V
    c = t/(2*mp.pi); lg = {}
    def f(x):
        if x not in lg: lg[x] = c*mp.log(Nm + x)
        return lg[x]
    def gt(n, m, k): return f(n+m+k) - f(n+k) - f(n+m) + f(n)
    cnt = Lp - sum(sh)
    if cnt <= 0: return None
    if order == 2: ph = [gt(n, sh[0], sh[1]) for n in range(cnt)]
    else: ph = [gt(n+sh[2], sh[0], sh[1]) - gt(n, sh[0], sh[1]) for n in range(cnt)]
    S = abs(mp.fsum(mp.expjpi(2*(p - mp.floor(p))) for p in ph))
    return float(S/rhs), float(S), float(rhs)
if __name__ == '__main__':
    random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 5); ncase = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    for order, N in ((2, 5*10**4), (3, 2*10**4)):
        L = N//75; worst = (0, None); n = 0
        for i in range(ncase):
            Lp = random.choice([L, L, random.randint(3, L)])
            lam = 10**random.uniform(-3.5, math.log10(2.5))
            if order == 2: sh = (random.randint(1, max(1, Lp//3)) if random.random() < .5 else random.randint(1, 5), random.randint(1, 5) if random.random() < .5 else random.randint(1, max(1, Lp//3)))
            else: sh = tuple(random.randint(1, 4) if random.random() < .6 else random.randint(1, max(1, Lp//4)) for _ in range(3))
            r = run(order, N, L, Lp, sh, lam)
            if r is None: continue
            n += 1
            if r[0] > worst[0]: worst = (r[0], (Lp, sh, round(lam, 5), round(r[1], 4), round(r[2], 3)))
        print(f'order {order}, N = {N}, L = {L}: {n} exact sums, max |S|/bound = {worst[0]:.5f} at (L\', shifts, lam, |S|, bound) = {worst[1]}', flush=True)
