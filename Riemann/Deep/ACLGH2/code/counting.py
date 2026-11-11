# counting.py: counting function N(kappa) = #{n : kappa_n > kappa} in the computed models, and the independent-rung estimate
# N_rung(kappa) = sum over frequencies n of C_S (S-smooth rationals, weight gamma(n)) of 2n [ |gamma(n)|/sqrt(n) > kappa ]
# (the term 2 gamma(n) cos(2 pi n u w) on [0,1]^2 has about 2n singular values of size |gamma(n)|/sqrt(n): finite cosine transform).
import sys, json, glob, math
KS = [0.3, 0.2, 0.15, 0.1, 0.07, 0.05, 0.03, 0.02, 0.01, 0.003, 0.001]
def terms(primes, nmax):
    t = {1.0: 1.0}
    for p in primes:
        new = {}
        for n, gm in t.items():
            k = -1
            while True:
                m = n*p**k
                if m > nmax: break
                new[m] = new.get(m, 0.0) + gm*((1 - 1/p) if k >= 0 else -1/p)
                k += 1
        t = new
    return t
def nrung(primes, ks):
    nmax = 4.0/ks**2
    return sum(2*n for n, g in terms(primes, nmax).items() if abs(g)/math.sqrt(n) > ks)
if __name__ == '__main__':
    for S in ([2], [2, 3], [2, 3, 5], [2, 3, 5, 7], [2, 3, 5, 7, 11]):
        tag = ''.join(map(str, S))
        print('S = {inf,%s}' % ','.join(map(str, S)))
        print('   %-26s' % 'kappa:', ' '.join('%7g' % k for k in KS))
        files = sorted(glob.glob('../logs/p%s_n*.json' % tag), key=lambda f: int(f.split('_n')[1].split('.')[0])) + sorted(glob.glob('../logs/L_p%s_t*.json' % tag))
        for f in files:
            if 'killed' in f: continue
            d = json.load(open(f)); k = d['kappa']; kmin = k[-1]
            row = [('%7d' % sum(1 for x in k if x > ks)) if ks > kmin else '      -' for ks in KS]
            print('   %-26s' % f.split('/')[-1][:-5], ' '.join(row))
        print('   %-26s' % 'N_rung', ' '.join('%7d' % nrung(S, ks) for ks in KS))
