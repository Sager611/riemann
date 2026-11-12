# explore44 c6: the small-ball term of the transport form of Littlewood's lemma, A(eps) = E log+(eps/|F|), lift vs [1000,2000], zeta vs DH,
# with the Herbst bound for zeta: A(eps) <= int_{log 1/eps}^inf exp(-r^2/(2 L^2)) dr.
import numpy as np, json, math
D = '/Users/adriansager/playground/riemann/.scratch/explore44/'
L2 = {0.6: 2.77591, 0.7: 1.77818, 0.8: 1.23836}
for sig in [0.6, 0.7, 0.8]:
    lz = np.load(D + 'lift_zeta_%g.npy' % sig).astype(float); ld = np.load(D + 'lift_dh_%g.npy' % sig).astype(float)
    vz, vd = np.load(D + 'vals_T1000_s%g.npy' % sig).astype(float)
    row = []
    for R in [1, 2, 3, 4]:
        A = lambda x: float(np.mean(np.maximum(-R - x, 0)))
        herb = math.sqrt(math.pi*L2[sig]/2)*math.erfc(R/math.sqrt(2*L2[sig]))
        row.append('R=%d: zeta lift %.1e T %.1e Herbst<=%.1e | DH lift %.1e T %.1e' % (R, A(lz), A(vz), herb, A(ld), A(vd)))
    print('sigma', sig); [print('  ' + x) for x in row]
    # empirical W1 (1D) between the laws of log|F| on [T,2T] and on the lift = int |F_T - F_lift|
    for name, a, b in [('zeta', vz, lz), ('DH', vd, ld)]:
        q = np.linspace(0.0005, 0.9995, 2000)
        w1 = float(np.mean(np.abs(np.quantile(a, q) - np.quantile(b, q))))
        print('  W1(log|%s| on [T,2T], lift) ~ %.4f ; mean difference %.4f' % (name, w1, a.mean() - b.mean()))
