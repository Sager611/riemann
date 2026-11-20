#!/usr/bin/env python3
# alpaut probe P2: case Z (PropZ, AV/HeightAndVanishing/Counts.lean:49-52): V(l, e, trivial) = 0 for l >= 5, e >= 2 even, (l-1) | 2e.
# Classical reduction (PROBES.md P2): V = 0 <=> the omega^(-e)-eigenspace A_{omega^(-e)} of Cl(Q(zeta_l))[l^inf] vanishes;
# omega^(-e) = 1 if (l-1)|e, else omega^((l-1)/2) (quadratic; e even forces l = 1 mod 4; A = Cl(Q(sqrt l))[l^inf]).
# (i) case list for l in {5,7,11,13}, e <= 12, plus h(Q(zeta_l)) by PARI bnfinit (all eigenspaces at once);
# (ii) h(Q(sqrt l)) by PARI quadclassunit for l = 1 mod 4 < LMAX: l does not divide h;
# (iii) independent route: Leopoldt reflection rank A_{omega^{(l-1)/2}} <= rank A_{omega^{(l+1)/2}}; Herbrand-Ribet:
#       A_{omega^{(l+1)/2}} != 0 <=> l | B_{(l-1)/2}.  So l not dividing B_{(l-1)/2} forces the vanishing as well.
import sys, cypari2
pari = cypari2.Pari(); pari.allocatemem(64 * 10**6)
LMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
print('(i) cases l in {5,7,11,13}, 2 <= e <= 12 even, (l-1) | 2e:')
for l in (5, 7, 11, 13):
    hz = pari('bnfinit(polcyclo(%d), 1).no' % l)
    cs = []
    for e in range(2, 13, 2):
        if (2 * e) % (l - 1) == 0:
            cs.append('e=%d:%s' % (e, 'trivial' if e % (l - 1) == 0 else 'quadratic Q(sqrt %d)' % l))
    print('  l=%d h(Q(zeta_l))=%s  %s' % (l, hz, ', '.join(cs)))
rows = bad = badB = 0; maxh = 0
for l in pari('primes([5, %d])' % LMAX):
    l = int(l)
    if l % 4 != 1: continue
    h = int(pari('quadclassunit(%d)[1]' % l))
    B = pari('bernfrac(%d)' % ((l - 1) // 2))
    Bm = int(pari('lift(Mod(numerator(%s), %d) / Mod(denominator(%s), %d))' % (B, l, B, l)))
    rows += 1; maxh = max(maxh, h); bad += (h % l == 0); badB += (Bm == 0)
    if l < 60: print('  l=%d h(Q(sqrt l))=%d B_{(l-1)/2} mod l = %d' % (l, h, Bm))
print('(ii) l = 1 mod 4, 5 <= l < %d: %d primes, max h = %d, primes with l | h: %d' % (LMAX, rows, maxh, bad))
print('(iii) primes with l | B_{(l-1)/2}: %d (Ankeny-Artin-Chowla / Mordell: none known below 2e11)' % badB)
