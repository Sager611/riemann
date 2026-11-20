# explore62 (b): floor of the zero density over S + M1 only (the fourth-order rows are carried by zeta windows, Thm 3):
# continuum of thin-band Gaussian blocks: band w in [0,pi] carries energy dw/pi, amplitude A(w)=1/s(w)^2, length s^2 dw/pi;
# minimise (1/pi^2) int w s^2 dw  s.t. (1/pi) int s dw = M1*sqrt(pi/2), (1/pi) int s^2 dw <= 1.  KKT: s = a/(w+b).
# Then a discrete J-band LP (b4.py columns without the M4 and R rows) as a check, exact-rational feasibility of its witness.
import numpy as np, sys
from scipy.optimize import linprog, brentq
from fractions import Fraction as Fr
def cont(M1):
    mu = M1*np.sqrt(np.pi/2)
    g = lambda b: np.sqrt(b*(np.pi+b))/np.pi*np.log1p(np.pi/b) - mu
    b = brentq(g, 1e-14, 1e6); a2 = b*(np.pi+b)
    dens = a2/np.pi**2*(np.log1p(np.pi/b) - np.pi/(np.pi+b))
    kurt = ((np.pi+b)**3 - b**3)/(np.pi*a2)   # 3 int p A^2 = (3/pi) int 1/s^2
    return b, dens, kurt
print("M1   beta     floor(S+M1)  background kurtosis  (Gaussian M1=sqrt(2/pi)=0.7979 -> 1/sqrt3=0.5774)")
for M1 in (0.7978, 0.7, 0.6, 0.45, 0.3, 0.2, 0.1, 0.05):
    b, d, k = cont(M1); print(f"{M1:.4f} {b:.5f}  {d:.5f}      {k:9.2f}")
# zeta: M1 = <|Z|>/<Z^2>^(1/2) ~ c L^(-1/4) (Heath-Brown; Heap-Radziwill-Soundararajan), 0.45 at T=3e12 (ceilaudit7, rough), M4 ~ L^2/(2pi^2)
L0 = np.log(3e12/(2*np.pi))
for L in (L0, 1e2, 1e3, 1e4, 1e6):
    M1 = 0.45*(L0/L)**0.25; b, d, k = cont(M1)
    print(f"L={L:.3g}: M1~{M1:.3f} floor~{d:.4f}; background kurtosis {k:.1f} vs zeta M4~L^2/(2pi^2)={L*L/(2*np.pi**2):.1f}")
J = int(sys.argv[1]) if len(sys.argv) > 1 else 32; M1 = 0.45
edges = np.linspace(0, np.pi, J+1); amps = np.exp(np.linspace(np.log(1e-3), np.log(1e3), 121))
meta = []
for j in range(J):
    w1, w2 = edges[j], edges[j+1]; l2 = (w2**3-w1**3)/(3*(w2-w1))
    for A in amps: meta.append((j, A, np.sqrt(l2)/np.pi))
n = len(meta); Aeq = []; beq = []
for j in range(J): Aeq.append([m[1] if m[0] == j else 0.0 for m in meta]); beq.append(1.0/J)
Aeq.append([np.sqrt(2/np.pi)*np.sqrt(m[1]) for m in meta]); beq.append(M1)
res = linprog([m[2] for m in meta], A_ub=[[1.0]*n], b_ub=[1.0], A_eq=np.array(Aeq), b_eq=np.array(beq), bounds=[(0, None)]*n, method='highs')
print(f"discrete LP J={J} M1={M1}: min zero density {res.fun:.5f}, live length {res.x.sum():.4f}; continuum {cont(M1)[1]:.5f}")
