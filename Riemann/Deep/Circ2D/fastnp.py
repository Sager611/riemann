# fastnp.py -- float64/numpy version of the reconstruction (for dense grids and the infinite system)
import numpy as np

A = (0, 1, 3, 4, 7, 9)
b = np.sqrt(3.0) / 2; h = 0.4; B = 4.0 / 3.0
PI = np.pi

def laurent_P():
    rho = np.exp(1j * PI / 6)
    poly = {0: 1 + 0j}
    for a in A:
        f = {0: 2 + 0j, 1: -rho ** (-a), -1: -rho ** a}
        new = {}
        for k1, v1 in poly.items():
            for k2, v2 in f.items():
                new[k1 + k2] = new.get(k1 + k2, 0) + v1 * v2
        poly = new
    return np.array([poly.get(l, 0j) for l in range(-6, 7)])  # index l+6

PL = laurent_P()
TL = np.arange(-6, 7) / 6.0

def nodes(nmax, with0=False):
    ns = [n for n in range(1, int(nmax) + 1) if n % 12 in A]
    return ([0] if with0 else []) + ns

def Pval(s):
    s = np.asarray(s, dtype=float)
    r = np.ones_like(s)
    for a in A:
        r = r * (2 * np.sin(PI * (s - a) / 12)) ** 2
    return r

def node_QD(n):
    e = PL * np.exp(1j * PI * TL * n)
    P2 = np.real(np.sum(e * (1j * PI * TL) ** 2)); P3 = np.real(np.sum(e * (1j * PI * TL) ** 3))
    return P2 / 2, P3 / (3 * P2)

class Quad:
    def __init__(self, gl=48, sub=4):
        x, w = np.polynomial.legendre.leggauss(gl)
        T, W, Jp = [], [], []
        for j in range(1, 7):
            for q in range(sub):
                lo = (j - 1) / 6 + q / (6 * sub); hi = lo + 1 / (6 * sub)
                T.append((lo + hi) / 2 + (hi - lo) / 2 * x); W.append((hi - lo) / 2 * w); Jp.append(np.full(gl, j))
        self.t = np.concatenate(T); self.w = np.concatenate(W); self.j = np.concatenate(Jp)
        self.lam = 1j / (b * (self.t + 1j * h)); self.z = -B / (self.t + 1j * h) - 1j * h
        self.at = np.arange(0, 7) / 6.0
        self.alam = 1j / (b * (self.at + 1j * h)); self.az = -B / (self.at + 1j * h) - 1j * h
        self.amass = np.array([PL[6]] + [2 * PL[6 + j] for j in range(1, 7)])  # e_j P_j for C=1

    def basis(self, n):
        """(c-density, d-density) at quadrature points for node n."""
        al = np.zeros(7, complex); be = np.zeros(7, complex)
        for j in range(1, 7):
            l = np.arange(j, 7)
            e = PL[6 + l] * np.exp(1j * PI * (l / 6.0) * n)
            al[j] = e.sum(); be[j] = (e * l / 6.0).sum()
        ph = np.exp(-1j * PI * self.t * n)
        cd = -2 * PI ** 2 * ph * (be[self.j] - self.t * al[self.j])
        dd = 2j * PI * ph * al[self.j]
        return cd, dd

    def density(self, L):
        f = np.zeros_like(self.t, dtype=complex)
        for n, (c, d) in L.items():
            cd, dd = self.basis(n)
            f += c * cd + d * dd
        return f

    def K(self, dens, C, s, k=0):
        """K^{(k)}(s) for array s."""
        s = np.atleast_1d(np.asarray(s, float))
        out = np.zeros(s.shape)
        for i0 in range(0, len(s), 1000):
            ss = s[i0:i0 + 1000]
            E = np.exp(1j * PI * np.outer(ss, self.z))
            v = (E * (self.w * self.lam * (1j * PI * self.z) ** k * dens)).sum(axis=1)
            if C != 0:
                Ea = np.exp(1j * PI * np.outer(ss, self.az))
                v = v + C * (Ea * (self.alam * (1j * PI * self.az) ** k * self.amass)).sum(axis=1)
            out[i0:i0 + 1000] = v.real
        return out

    def p(self, dens, C, s, k=0):
        s = np.atleast_1d(np.asarray(s, float))
        out = np.zeros(s.shape)
        for i0 in range(0, len(s), 1000):
            ss = s[i0:i0 + 1000]
            E = np.exp(1j * PI * np.outer(ss, self.t))
            v = (E * (self.w * (1j * PI * self.t) ** k * dens)).sum(axis=1)
            if C != 0:
                Ea = np.exp(1j * PI * np.outer(ss, self.at))
                v = v + C * (Ea * ((1j * PI * self.at) ** k * self.amass)).sum(axis=1)
            out[i0:i0 + 1000] = v.real
        return out

def R_val(L, C, s):
    s = np.asarray(s, float)
    r = np.full(s.shape, float(C))
    for n, (c, d) in L.items():
        u = s - n
        r = r + c / u ** 2 + d / u
    return r

def read_table(path):
    import re
    tab = {}
    for line in open(path):
        m = re.findall(r'-?\d+', line)
        if len(m) == 5:
            tab[int(m[0])] = tuple(int(x) * 1e-10 for x in m[1:])
    return tab

FIXED = {1: ((1.0, 0.44), -0.013), 2: ((0.0, -0.368), 0.017)}

