# windows.py: window sums W(a,b] = sum_{a<n<=b} kappa_n m_n(R) from semi_pair JSON output (index n from 0, as kappa_0 >= kappa_1 >= ...).
# usage: windows.py file.json [file2.json ...]
import sys, json
WIN = [(-1, 4), (4, 10), (10, 20), (20, 50), (50, 100), (100, 200), (200, 400), (400, 800)]
for f in sys.argv[1:]:
    d = json.load(open(f)); k = d['kappa']; info = d['info']
    size = ('nmax=%g' % info['nmax']) if 'nmax' in info else ('tmax=%g depth=%g' % (info['tmax'], info['depth']))
    print(f.split('/')[-1], 'S=', info.get('primes'), size, 'nk=', len(k), 'k_last=%.3g' % k[-1])
    for R, m in d['mass'].items():
        row = []
        for a, b in WIN:
            if b + 1 > len(k): row.append('   -   '); continue
            row.append('%.3g' % sum(k[i]*m[i] for i in range(a + 1, b + 1)))
        print('  R=%-5s' % R, ' '.join('(%d,%d]:%s' % (a, b, r) for (a, b), r in zip(WIN, row)))
