# independent count: sign changes of the raw Wronskians [u,s](X), [u,c](X) at finite X (no phase unwrapping)
# versus the Pruefer-phase eigenvalues (Theta in pi Z / pi/2 + pi Z).  Usage: python3 signchange.py <eig_even> <eig_odd> <wr_file> <Emax>
import sys
def load_eig(fn, Emax):
    return [(int(l.split()[0]), float(l.split()[1])) for l in open(fn) if not l.startswith('#') and float(l.split()[1]) <= Emax]
Emax = float(sys.argv[4])
e0 = load_eig(sys.argv[1], Emax); e1 = load_eig(sys.argv[2], Emax)
W = [tuple(map(float, l.split())) for l in open(sys.argv[3])]
def signch(col):
    res = []
    for i in range(1, len(W)):
        a, b = W[i-1][col], W[i][col]
        if (a < 0) != (b < 0):
            E0, E1 = W[i-1][0], W[i][0]; res.append(E0 + (E1-E0)*a/(a-b))
    return res
sc = signch(1); cc = signch(2)
print(f'E <= {Emax}: even #eig(Theta) = {len(e0)}, #sign changes [u,s](X) = {len(sc)}; odd #eig(Theta) = {len(e1)}, #sign changes [u,c](X) = {len(cc)}')
for name, ev, sg in (('even', e0, sc), ('odd', e1, cc)):
    m = min(len(ev), len(sg))
    mx = max(abs(sg[i]-ev[i][1]) for i in range(m))
    print(f'{name}: max |signchange - eigenvalue| over the first {m} = {mx:.4f} (finite-X shift ~ (E^2+kappa^2)/(2 omega X) / Theta\'(E))')
allv = sorted([(E, 'e') for _, E in e0] + [(E, 'o') for _, E in e1])
print('even/odd strictly interlace:', all(allv[i][1] != allv[i+1][1] for i in range(len(allv)-1)), ' first labels:', ''.join(x[1] for x in allv[:10]))
