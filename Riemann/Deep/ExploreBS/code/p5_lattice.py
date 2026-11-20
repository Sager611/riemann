# p5_lattice.py (explore57): Birman-Schwinger rank of the lattice-minorant certificate vs the phase-space (Weyl/Bohr-Jessen) count.
# Record method (Chekmachev, Zenodo 23848597, audited ChekAudit/WeilZeta): positivity of Q on the window of length l = ln c is
# reduced to a dual check lam_max(B^T P^-1 B) <= 1 (his Thm 5.1 = Birman-Schwinger for a finite-rank attractive part) whose
# order is |K-| + bad even nodes (the lattice nodes where his lattice symbol is negative).  Here: (i) Bohr-Jessen phase-space
# count W(c) = (l/2pi)|{t: sigma_c(t) < 0}| ~ 2 l E[exp X_c] (sigma_c = Omega - X_c, X_c = 2 sum_{n<c} Lambda(n) n^-1/2 cos(t ln n)),
# (ii) the number of negative values of sigma_c on the odd (Dirichlet) and even (Neumann) node systems r = pi(2k+1)/l, 2 pi k/l,
# k < kscan; compared with the certificate sizes (Km, order) logged by WeilZeta (val*.log, z*.log) and ChekAudit (c = 47).
import sys, math, numpy as np
from scipy.special import digamma
def atoms(c):
    out = []
    for n in range(2, int(math.ceil(c))):
        m = n; p = next(q for q in range(2, n + 1) if n % q == 0)
        while m % p == 0: m //= p
        if m == 1: out.append((n, math.log(p)/math.sqrt(n)))
    return out
def EexpX(c):
    tot = 1.0; th = np.linspace(0, 2*np.pi, 20001)[:-1]
    ps = sorted({n for n, w in atoms(c) if all(n % q for q in range(2, int(n**0.5) + 1))})
    for p in ps:
        s = np.zeros_like(th); k = 1
        while p**k < c:
            s += 2*math.log(p)*p**(-k/2)*np.cos(k*th); k += 1
        tot *= np.exp(s).mean()
    return tot
def sigma(r, A):
    om = digamma(0.25 + 0.5j*r).real - math.log(math.pi)
    X = np.zeros_like(r)
    for n, w in A: X += 2*w*np.cos(r*math.log(n))
    return om - X
CERT = {5: (8, 9), 10: (21, 33), 16: (105, 133), 20: (214, 311), 25: (423, 558), 29: (None, 1084), 47: (6301, 7170)}
KSCAN = {5: 201, 10: 14797, 16: 313858, 20: 7319040, 25: 29076033}
def main():
  cs = [int(x) for x in sys.argv[1].split(',')]; lat = int(sys.argv[2])
  for c in cs:
      l = math.log(c); A = atoms(c); W = 2*l*EexpX(c)
      line = f'c={c:3d} l={l:.4f} W_BJ={W:10.4g} W_BJ/2={W/2:9.4g}'
      if c in CERT: line += f' | cert Km={CERT[c][0]} order={CERT[c][1]}'
      if c <= lat and c in KSCAN:
          neg_o = 0; neg_e = 0; K = KSCAN[c]; ch = 1_000_000
          for a in range(0, K, ch):
              k = np.arange(a, min(a + ch, K), dtype=float)
              neg_o += int(np.count_nonzero(sigma(np.pi*(2*k + 1)/l, A) < 0))
              neg_e += int(np.count_nonzero(sigma(2*np.pi*k/l, A) < 0))
          line += f' | sigma_c<0 at odd nodes {neg_o}, even nodes {neg_e} (k < {K})'
      print(line, flush=True)

if __name__ == '__main__':
    main()
