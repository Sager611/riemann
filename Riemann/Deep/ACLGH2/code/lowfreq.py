# lowfreq.py: low-frequency spectral data of the pair. For each file: the modes with kappa*m(R) > 1e-4 (R = 20) or m(20) > 0.2,
# W(kappa*) = sum_{kappa_n <= kappa*} kappa_n m_n(R) at sample thresholds, and kappa_nat(tau) = min{kappa_n : W(kappa_n) > tau}
# (the angle of the last mode whose inclusion keeps the tail above tau; the onset is n* = N(kappa_nat), N = counting function).
import sys, json
KS = [0.3, 0.2, 0.1, 0.05, 0.03, 0.02, 0.01, 0.005, 0.002, 0.001, 1e-4]
for f in sys.argv[1:]:
    d = json.load(open(f)); k = d['kappa']; I = d['info']
    size = ('nmax=%g' % I['nmax']) if 'nmax' in I else ('tmax=%g depth=%g N=%d' % (I['tmax'], I['depth'], I['N']))
    print('%s S=%s %s t=%g' % (f.split('/')[-1], I['primes'], size, I.get('t', 1)))
    m20 = d['mass']['20.0']; m5 = d['mass']['5.0']; xm = d.get('xmean')
    sel = [n for n in range(len(k)) if k[n]*m20[n] > 1e-4 or m20[n] > 0.2]
    print('  modes (n kappa m20 m5 [xmean]):', ' '.join('(%d %.4f %.2f %.2f%s)' % (n, k[n], m20[n], m5[n], (' %.1f' % xm[n]) if xm else '') for n in sel[:40]))
    for R in ('5.0', '20.0'):
        m = d['mass'][R]
        W = lambda ks: sum(k[n]*m[n] for n in range(len(k)) if k[n] <= ks)
        print('  R=%s W(kappa*):' % R, ' '.join('%g:%.4f' % (ks, W(ks)) for ks in KS))
        order = sorted(range(len(k)), key=lambda n: k[n])
        out = []
        for tau in (1e-2, 1e-3):
            acc = 0.0; nat = None
            for n in order:
                acc += k[n]*m[n]
                if acc > tau: nat = n; break
            out.append('tau=%g: kappa_nat=%s (n=%s, m=%.2f)' % (tau, ('%.4f' % k[nat]) if nat is not None else '-', nat, m[nat] if nat is not None else 0))
        print('   ', '; '.join(out))
