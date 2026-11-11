# profile.py: kappa_n at sample indices, kappa_n sqrt(n), and partial sums P(n0) = sum_{n<=n0} kappa_n m_n(R) (index from 0).
import sys, json, math
IDX = [0, 2, 5, 10, 20, 50, 100, 200, 300, 400, 600]
for f in sys.argv[1:]:
    d = json.load(open(f)); k = d['kappa']; I = d['info']
    size = ('nmax=%g nodes=%d' % (I['nmax'], I['nodes'])) if 'nmax' in I else ('tmax=%g depth=%g N=%d' % (I['tmax'], I['depth'], I['N']))
    print('%s S=%s %s t=%g nk=%d' % (f.split('/')[-1], I['primes'], size, I.get('t', 1), len(k)))
    print('  kappa:', ' '.join('%d:%.3f' % (i, k[i]) for i in IDX if i < len(k)))
    print('  k*sqrt(n):', ' '.join('%d:%.2f' % (i, k[i]*math.sqrt(i)) for i in IDX[1:] if i < len(k)))
    for R, m in d['mass'].items():
        P = 0; out = []
        for i in range(len(k)):
            P += k[i]*m[i]
            if i in IDX: out.append('%d:%.3f' % (i, P))
        print('  R=%s partial:' % R, ' '.join(out), ' total(nk)=%.3f' % P)
