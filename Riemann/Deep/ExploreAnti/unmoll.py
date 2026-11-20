# explore55: unmollified holomorphic moments E zeta'(rho)/L -> 1/2, E zeta'(rho)^2/L^2 -> 1/6 (recipe), E|zeta'|^2/L^3 -> 1/12 (Gonek)
import numpy as np, re, mpmath as mp
from rs import zprime, theta
def stats(g, zp, tag):
    L = np.log(g.max()/(2*np.pi)); D = -1j*np.exp(-1j*thm(g))*zp
    print(tag, 'L=%.3f N=%d  E[z\']/L=%.4f%+.4fi  E[z\'^2]/L^2=%.4f%+.4fi  E|z\'|^2/L^3=%.4f' % (L, len(g), *(lambda a:(a.real,a.imag))(D.mean()/L), *(lambda a:(a.real,a.imag))(np.mean(D**2)/L**2), np.mean(abs(D)**2)/L**3))
thm = theta
g = np.loadtxt('../../../.scratch/explore55/zeros1.txt')
for a, b in [(1000, 2000), (10000, 20000), (50000, 100000)]:
    stats(g[a:b], zprime(g[a:b]), 'zeros1[%d:%d]' % (a, b))
t0 = 267653395647; d = np.array([float(l) for l in open('../../../.scratch/explore55/zeros3.txt') if re.fullmatch(r'\s*[0-9]+\.[0-9]+\s*', l)])
zp, z0 = np.load('../../../.scratch/explore55/zeros3_zp10000.npy'); mp.mp.dps = 40
th = np.array([float(mp.fmod(mp.siegeltheta(mp.mpf(t0)+mp.mpf(str(x))), 2*mp.pi)) for x in d])
thm = lambda gg: th
stats(t0 + d, zp, 'zeros3 (1e12)')
