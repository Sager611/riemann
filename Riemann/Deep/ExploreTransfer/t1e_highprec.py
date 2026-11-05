# t1e_highprec.py -- explore3: validate the double-precision spectra (t <= 27) and extend to t = 30, 40 with
# multiprecision QR eigenvalues (Arb acb_mat.eig(algorithm='approx')).
import sys, time
import numpy as np
from flint import acb, acb_mat, ctx
sys.path.insert(0, '.')
from t1_flint import build_M
from t1b_branches import build_acb

def top_mp(b, K, prec):
    M = build_acb(b, K, prec)
    ctx.prec = prec
    ev = M.eig(algorithm='approx')
    m = sorted((abs(complex(e)) for e in ev), reverse=True)
    return m

if __name__ == '__main__':
    t0 = time.time()
    print('validation at t=25, K=150: numpy(double) vs Arb QR (300 bits)')
    for s in (0.25, 0.5, 1.0):
        b = complex(s, 25.0)
        mnp = np.sort(abs(np.linalg.eigvals(build_M(b, 150, prec=514))))[::-1][:6]
        mmp = top_mp(b, 150, 300)[:6]
        print('  sigma=%.2f numpy: %s | Arb: %s' % (s, ' '.join('%.4f' % x for x in mnp), ' '.join('%.4f' % x for x in mmp)), flush=True)
    print('extension (Arb QR, 400 bits)')
    for t, Ks in ((30.0, (180, 220)), (40.0, (220, 260))):
        for s in (0.25, 0.5, 1.0):
            row = []
            for K in Ks:
                t1 = time.time()
                m = top_mp(complex(s, t), K, 400)
                n1 = sum(1 for x in m if x > 1); n02 = sum(1 for x in m if x > 0.2)
                row.append('K=%d #>1=%d #>0.2=%d top8=%s (%.0fs)' % (K, n1, n02, ','.join('%.3f' % x for x in m[:8]), time.time() - t1))
            print('  b=%.2f+%.0fi: ' % (s, t) + ' | '.join(row), flush=True)
    print('time %.1fs' % (time.time() - t0))
