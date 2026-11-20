# explore55: falsification attempt of Conjecture H beyond theta = 1/2 (B - u^2 = -1/12 for every mollifier?) on zeros 50001-100000.
import numpy as np, sys, resource
from moll import feats, bound_complex, bound_real
from rs import zprime
g = np.loadtxt(sys.argv[1])[50000:100000]; zp = zprime(g); T = g[-1]; L = np.log(T/(2*np.pi))
for th in [0.5, 0.75, 1.0]:
    D, Mp, _ = feats(g, zp, T**th, 1, chunk=100); Y = D*Mp[:, 0]; z = Y/(Y.mean()/abs(Y.mean()))
    u = abs(Y.mean())/L; A = np.mean(abs(z)**2)/L**2; B = np.mean(z**2)/L**2; P1 = 1.0
    print(f"theta={th}: y={T**th:.0f} u={u:.4f} (BHB-form {0.5+th/2:.4f}) A={A:.4f} B={B.real:.4f}{B.imag:+.4f}i B-u^2={B.real-u*u:.4f} A-u^2={A-u*u:.4f} CGG={bound_complex(Y[:, None]):.4f} 2D={bound_real(Y[:, None]):.4f} rss={resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20:.0f}MB", flush=True)
