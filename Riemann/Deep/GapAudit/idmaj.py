"""gapaudit: (1) exact identity of DERIVATION 3 at random points (mpmath, 50 digits);
(2) sanity check of the tail majorant |Gam_ij| <= 4(i+j+8)^2((i+2r)(j+2r))^(2r) gam^(i+j) on the computed twisted arrays
(certmid's own arrays, r = 2, kappa = 3.324, D = 80): max ratio by anti-diagonal band.  Imports GapBar modules read-only."""
import sys, os, math, random
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'GapBar'))
import mpmath as mp
mp.mp.dps = 50
random.seed(1)
worst = mp.mpf(0)
for _ in range(2000):
    X1, X2, X3, X4 = [mp.mpf(random.random()) * mp.mpf('0.25') for _ in range(4)]
    sg = mp.mpf(random.uniform(-1, 1)); c = mp.mpf(random.uniform(0, 12))
    E = lambda t: mp.e ** (1j * t)
    U = E(c * (X1 - X3)) + E(c * (X2 - X4)) - E(c * (X1 - X4)) * E(c * sg) - E(c * (X2 - X3)) * E(-c * sg)
    lhs = E(c * ((X3 + X4) - (X1 + X2)) / 2) * U
    rhs = 4 * mp.sin(c * (X1 - X2 + sg) / 2) * mp.sin(c * (X3 - X4 + sg) / 2)
    worst = max(worst, abs(lhs - rhs))
print('identity: max |e^{ic(r2-r1)/2} U_L conj U_R - 4 sin sin| over 2000 random points =', mp.nstr(worst, 5))
from flint import arb
import gcert as GC, gtw
kap = GC.dec('3.324'); GC.TH = GC.dec('0.2499'); r = 2
gam = int(float(2 * GC.A_(kap) * arb.pi().mid() + 3)) + 1
D = int(sys.argv[1]) if len(sys.argv) > 1 else 80
GC.set_prec(128 + int(D * math.log2(max(gam / 4, 2))) + 32)
c = GC.A_(kap) * arb.pi()
Gam = gtw.gamma_arrays_tw(c, D, gam, r)
print('gam =', gam, ' 2c+2 =', float((2 * c + 2).mid()), ' D =', D, ' IMCHECK =', gtw.IMCHECK)
for key in 'abc':
    G = Gam[key]; band = {}
    for i in range(D):
        for j in range(D):
            bd = arb(4 * (i + j + 8) ** 2 * ((i + 2 * r) * (j + 2 * r)) ** (2 * r)) * arb(gam) ** (i + j)
            q = float((abs(G[i][j]) / bd).upper())
            k = (i + j) // 20
            band[k] = max(band.get(k, 0.0), q)
    print('array', key, ': max |Gam_ij|/majorant by band i+j in [20k,20k+20):', ' '.join('%d:%.2e' % (k, band[k]) for k in sorted(band)))
