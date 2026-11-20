# hdx_garland.py (explore64): Garland/IMS localization of Weil's form on I = [-a, a], a = (ln c)/2, into links of half-width a'.
# Identity (exact, translation-invariant links chi_b = chi(. - b), int chi^2 = 1, m = chi*chi~):
#   Q_a(f) = int_R Q(chi_b f) db + R(f),  R(f) = 2 int_0^{2a} (1-m(z)) (2cosh(z/2) - kappa(z)) C_f(z) dz - 2 sum_n w_n (1-m(ln n)) C_f(ln n).
# Hard links chi = 1_[-a',a']/sqrt(2a'): 1 - m(z) = min(z/(2a'), 1).  Garland bound: Q_a >= e1(c') + lambda_min(R).
# float64 exploration (Legendre parity basis, ../ExploreFunctional/code/xform.py); not a certificate.
import sys, math, numpy as np
sys.path.insert(0, 'Riemann/Deep/ExploreFunctional/code')
from xform import Form

def onem(z, ap, kind='hard'):
    z = np.asarray(z, float)
    if kind == 'hard': return np.minimum(z/(2*ap), 1.0)
    # smooth link chi = cos(pi x/(2a'))/sqrt(a') on [-a',a']: m(z) = autocorrelation (closed form), z in [0, 2a']
    t = np.minimum(z, 2*ap); u = math.pi*t/(2*ap)
    m = ((2*ap - t)*np.cos(u) + (2*ap/math.pi)*np.sin(u))/(2*ap)
    return 1 - np.where(z >= 2*ap, 0.0, m)

def Rmat(F, ap, n, sector, kind='hard', nz=None, nx=None):
    a = F.a; nz = nz or 2*n + 60; nx = nx or 2*n + 30
    R = np.zeros((n, n))
    cuts = [0.0] + ([2*ap] if 2*ap < 2*a else []) + [2*a]
    for lo, hi in zip(cuts[:-1], cuts[1:]):
        zs, wz = F.gauss(lo, hi, nz)
        for z, w in zip(zs, wz):
            Cz = F.Cmat(z, n, sector, nx); Cz = (Cz + Cz.T)/2
            kz = onem(z, ap, kind)*((2*np.cosh(z/2) if F.pole else 0.0) - F.zkap(z)/z)
            R += 2*w*kz*Cz
    for h, wt in F.atoms:
        Ch = F.Cmat(h, n, sector, nx)
        R -= wt*float(onem(h, ap, kind))*(Ch + Ch.T)
    return R

def eigs_Q(c, n, sector, kind='zeta'):
    F = Form(c, kind); Q, _ = F.Q(n, sector)
    return np.linalg.eigvalsh((Q + Q.T)/2)

if __name__ == '__main__':
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    kind = sys.argv[2] if len(sys.argv) > 2 else 'zeta'
    print(f'# n={n} per sector, kind={kind}; rows: c, c_link, sector, e1(c) e2(c) [Galerkin], e1(c_link) e2(c_link), lmin R_hard, lmin R_smooth, garland bound (hard)')
    for c in [10, 20, 30, 47]:
        for frac in [0.5, 0.75, 0.9]:
            cl = c**frac; ap = math.log(cl)/2
            F = Form(c, kind)
            for sec in ['even', 'odd']:
                eQ = eigs_Q(c, n, sec, kind); eL = eigs_Q(cl, n, sec, kind)
                Rh = Rmat(F, ap, n, sec, 'hard'); Rs = Rmat(F, ap, n, sec, 'smooth')
                lh = np.linalg.eigvalsh(Rh)[0]; ls = np.linalg.eigvalsh(Rs)[0]
                print(f'{c:3d} {cl:8.3f} {sec:4s}  eQ {eQ[0]: .3e} {eQ[1]: .3e}  eL {eL[0]: .3e} {eL[1]: .3e}  Rh {lh: .4f}  Rs {ls: .4f}  G {eL[0]+lh: .4f}', flush=True)
