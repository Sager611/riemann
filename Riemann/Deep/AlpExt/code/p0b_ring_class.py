# alpext P0b: ring class characters. For an order O_f of conductor f in F = Q(sqrt-D0), ring class characters psi satisfy
# psi^c = psi^-1, so psi^c = psi <=> psi^2 = 1. Ring class groups Pic(O_f) = primitive forms of discriminant -f^2 D0;
# exponent > 2 <=> a ring class character of order >= 3 exists (field not in any CM field, not a Dirichlet product).
import importlib.util, sys
spec = importlib.util.spec_from_file_location('p0', 'Riemann/Deep/AlpExt/code/p0_cm_obstruction.py')
src = open('Riemann/Deep/AlpExt/code/p0_cm_obstruction.py').read().split('first = None')[0]
ns = {}; exec(src.split('# (3) class groups')[0], ns)  # imports
exec('from math import gcd, isqrt\n' + src.split('# (3) class groups of imaginary quadratic discriminants -D, D <= 400, via reduced forms; exponent > 2 <=> psi with psi^2 != 1')[1], ns)
forms, compose, reduce_form = ns['forms'], ns['compose'], ns['reduce_form']
for D0 in (4, 3, 7, 8, 11):
    row = []
    for f in range(1, 11):
        D = f * f * D0
        Fm = forms(D); I = reduce_form((1, D % 2, (D % 2 + D) // 4), D)
        exp2 = all(reduce_form(compose(g, g, D), D) == I for g in Fm)
        row.append(f'f={f}:h={len(Fm)}{"" if exp2 else "*"}')
    print(f'F = Q(sqrt-{D0}) (disc -{D0}):', ' '.join(row))
print('* = ring class group of exponent > 2: a ring class character with psi^c != psi, outside the CM-host step')
