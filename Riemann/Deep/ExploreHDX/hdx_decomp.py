# hdx_decomp.py (explore64): which part of the Garland remainder R = K(1-m) is negative?  R = R_pole + R_arch + R_atoms
# (hard links, 1-m(z) = min(z/2a',1)).  Also lambda_min of R on the complement of the first k basis modes (high frequency).
import sys, math, numpy as np
sys.path.insert(0, 'Riemann/Deep/ExploreFunctional/code'); sys.path.insert(0, 'Riemann/Deep/ExploreHDX')
from xform import Form
from hdx_garland import onem

def parts_R(F, ap, n, sector, kind='hard'):
    a = F.a; nz = 2*n + 60; nx = 2*n + 30
    Rp = np.zeros((n, n)); Ra = np.zeros((n, n)); Rt = np.zeros((n, n))
    cuts = [0.0] + ([2*ap] if 2*ap < 2*a else []) + [2*a]
    for lo, hi in zip(cuts[:-1], cuts[1:]):
        zs, wz = F.gauss(lo, hi, nz)
        for z, w in zip(zs, wz):
            Cz = F.Cmat(z, n, sector, nx); Cz = (Cz + Cz.T)/2; om = float(onem(z, ap, kind))
            if F.pole: Rp += 2*w*om*2*np.cosh(z/2)*Cz
            Ra -= 2*w*om*(F.zkap(z)/z)*Cz
    for h, wt in F.atoms:
        Ch = F.Cmat(h, n, sector, nx); Rt -= wt*float(onem(h, ap, kind))*(Ch + Ch.T)
    return Rp, Ra, Rt

lm = lambda A: np.linalg.eigvalsh((A + A.T)/2)[0]
if __name__ == '__main__':
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    kind = sys.argv[2] if len(sys.argv) > 2 else 'zeta'
    print(f'# n={n}, {kind}: lmin of R_pole, R_arch, R_atoms, R_pole+R_arch, R; then lmin R on span of modes k>=4 and k>=8 (high freq)')
    for c, frac in [(10, .5), (10, .9), (47, .5), (47, .9)]:
        F = Form(c, kind); ap = math.log(c**frac)/2
        for sec in ['even', 'odd']:
            Rp, Ra, Rt = parts_R(F, ap, n, sec); R = Rp + Ra + Rt
            hf = [lm(R[k:, k:]) for k in (4, 8)]
            print(f'{c:3d} {c**frac:7.3f} {sec:4s} pole {lm(Rp): .3f} arch {lm(Ra): .3f} atoms {lm(Rt): .3f} p+a {lm(Rp+Ra): .3f} R {lm(R): .3f} | hf4 {hf[0]: .3f} hf8 {hf[1]: .3f}', flush=True)
