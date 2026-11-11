# natsum.py: table of kappa_nat(S, R, tau) (lowfreq.py definition) and W(kappa*) for a list of Laurent-model files.
import sys, json
def nat(k, m, tau):
    acc = 0.0
    for n in sorted(range(len(k)), key=lambda n: k[n]):
        acc += k[n]*m[n]
        if acc > tau: return k[n]
    return None
print('%-28s %3s %6s | %-27s | %-27s | W20(.01) W20(.001) | tot5 tot20 | k0 xm0' % ('file', '|S|', 'pmax', 'R=5: nat(1e-2) nat(1e-3)', 'R=20: nat(1e-2) nat(1e-3)'))
for f in sys.argv[1:]:
    d = json.load(open(f)); k = d['kappa']; I = d['info']; P = I['primes']
    out = []
    for R in ('5.0', '20.0'):
        m = d['mass'][R]
        out.append('%12.5f %12.5f' % (nat(k, m, 1e-2) or -1, nat(k, m, 1e-3) or -1))
    m20 = d['mass']['20.0']; m5 = d['mass']['5.0']
    W = lambda ks, m: sum(k[n]*m[n] for n in range(len(k)) if k[n] <= ks)
    tot = lambda m: sum(k[n]*m[n] for n in range(len(k)))
    print('%-28s %3d %6s | %s | %s | %8.4f %8.4f | %5.2f %5.2f | %.4f %.1f' % (f.split('/')[-1][:-5], len(P) + 1, P[-1] if P else '-', out[0], out[1],
          W(0.01, m20), W(0.001, m20), tot(m5), tot(m20), k[0], d['xmean'][0]))
