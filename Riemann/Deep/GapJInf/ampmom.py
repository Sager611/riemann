"""gapjinf: even moments m_j = <(v-S)^{2j}> kappa^{2j} / <1> of one amplifier from mom.py forms; checks.
usage: ampmom.py forms.npz amp.json [v] -> prints moments, Hankel eigenvalues, explore53's h_J for the form in amp.json"""
import sys, json
from math import comb
from fractions import Fraction as Fr
import numpy as np


def load(npz, ampjson, v=None):
    D = np.load(npz); C = D['C']; R = D['R']; norm = D['norm']; kap = float(D['kap']); v0 = float(D['v0'])
    am = json.load(open(ampjson))
    A = np.array([float(x) for x in am['A']])
    y = R @ (A / norm)
    s = np.array([y @ Ck @ y for Ck in C])
    v = float(am['v']) if v is None else v
    t = -(v - v0) * kap
    K = len(s) - 1
    sv = np.array([sum(comb(k, i) * s[i] * t ** (k - i) for i in range(k + 1)) for k in range(K + 1)])
    return sv / sv[0], am, kap, v


if __name__ == '__main__':
    sv, am, kap, v = load(sys.argv[1], sys.argv[2], float(sys.argv[3]) if len(sys.argv) > 3 else None)
    K = len(sv) - 1
    m = sv[0::2]
    print('kappa=%.4f v=%.6f  odd/even check: s1=%.3e s3=%.3e' % (kap, v, sv[1], sv[3]))
    print('even moments m_j (j<=%d):' % (K // 2), ' '.join('%.10g' % x for x in m))
    J = K // 4
    H0 = np.array([[sv[2 * (i + k)] for k in range(J + 1)] for i in range(J + 1)])
    H1 = np.array([[sv[2 * (i + k + 1)] for k in range(J)] for i in range(J)])
    print('Stieltjes Hankel (x=xi^2) min eig/max eig: H0 %.3e  H1 %.3e' % (np.linalg.eigvalsh(H0)[0] / np.linalg.eigvalsh(H0)[-1], np.linalg.eigvalsh(H1)[0] / np.linalg.eigvalsh(H1)[-1]))
    Hm = np.array([[sv[i + k] for k in range(K // 2 + 1)] for i in range(K // 2 + 1)])
    ev = np.linalg.eigvalsh(Hm)
    print('Hamburger Hankel in xi (order %d): min eig %.3e, max %.3e' % (K // 2, ev[0], ev[-1]))
    a = [float(Fr(x)) for x in am['a']]
    Hs = -sum(a[j] * m[j] for j in range(len(a)))
    print('explore53 form (J=%d): H/c0 = %.8g, h_J = -a0 m0/sum_{j>=1} a_j m_j = %.8f' % (len(a) - 1, Hs, -a[0] / sum(a[j] * m[j] for j in range(1, len(a)))))
