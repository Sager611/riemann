# hahn_lemma.py — ceilaudit10 T2b: is there a Polya–Turan–Romik lemma for the continuous Hahn basis g_n (explore63 hahn_d2n.py:
# (2n+3) g_{n+1} = 8x g_n - (2n+1) g_{n-1}, g_0 = 1)?  Test T[z^k] = kappa_k g_k on random real-rooted p (degree 2..9), several kappa.
import random, mpmath as mp
mp.mp.dps = 40
def gcoeffs(nmax):
    G = [[mp.mpf(1)], [mp.mpf(0), mp.mpf(8)/3]]
    for n in range(1, nmax):
        c = [mp.mpf(0)]*(n+2)
        for k, v in enumerate(G[n]): c[k+1] += 8*v
        for k, v in enumerate(G[n-1]): c[k] -= (2*n+1)*v
        G.append([v/(2*n+3) for v in c])
    return G
G = gcoeffs(12)
lead = [G[k][k] for k in range(13)]
kap = {
 'monic': lambda k: 1/lead[k],
 'monic/(3/2)_k': lambda k: 1/(lead[k]*mp.rf(1.5, k)),
 'monic/k!': lambda k: 1/(lead[k]*mp.factorial(k)),
 'monic/(3/2)_k^2': lambda k: 1/(lead[k]*mp.rf(1.5, k)**2),
 'monic*k!/(3/2)_k^2': lambda k: mp.factorial(k)/(lead[k]*mp.rf(1.5, k)**2),
 'MP-like k!/(3/2)_k': lambda k: mp.factorial(k)/mp.rf(1.5, k),
}
random.seed(11)
for name, kf in kap.items():
    for sign in ('any', 'neg'):
        bad = 0; worst = 0; ex = None
        for trial in range(150):
            d = random.randint(2, 9); scale = random.choice([0.3, 1, 3, 10])
            p = [mp.mpf(1)]
            for _ in range(d):
                r = random.uniform(-1, 1)*scale if sign == 'any' else -random.uniform(0, 1)*scale
                q = [mp.mpf(0)]*(len(p)+1)
                for k, v in enumerate(p): q[k+1] += v; q[k] -= r*v
                p = q
            T = [mp.mpf(0)]*(d+1)
            for k, a in enumerate(p):
                w = kf(k)*a
                for j, cc in enumerate(G[k]): T[j] += w*cc
            im = max(abs(mp.im(z)) for z in mp.polyroots(T[::-1], maxsteps=400, extraprec=400))
            if im > 1e-15:
                bad += 1
                if im > worst: worst = im; ex = (d, scale)
        print('%-20s roots %-3s: %3d/150 images with non-real zeros (max|Im| %s, e.g. deg,scale %s)' % (name, sign, bad, mp.nstr(worst, 3), ex), flush=True)
