# explore55: same statistics as moll.py at large height t0 + delta (Odlyzko zeros3), phases reduced mod 2pi in mpmath.
import numpy as np, sys, json, mpmath as mp
from moll import mobius, bound_complex, bound_real
mp.mp.dps = 40; TWOPI = 2*np.pi
def prep(t0, nmax):
    ln = np.log(np.arange(1, nmax+1, dtype=float)); r = np.empty(nmax)
    for n in range(1, nmax+1):
        r[n-1] = float(mp.fmod(t0*mp.log(n), 2*mp.pi))
    return ln, r
def main(path, t0, thetas, J, nz=None, chunk=8):
    import re
    d = np.array([float(l) for l in open(path) if re.fullmatch(r'\s*[0-9]+\.[0-9]+\s*', l)])[:nz]; t = t0 + d; T = t.max(); L = np.log(T/TWOPI)
    N = int(np.floor(np.sqrt(t.min()/TWOPI))); assert N == int(np.floor(np.sqrt(T/TWOPI)))
    ymax = int(T**max(thetas)); nmax = max(N, ymax); ln, r = prep(t0, nmax)
    thm = np.array([float(mp.fmod(mp.siegeltheta(mp.mpf(t0)+mp.mpf(str(x))), 2*mp.pi)) for x in d])
    dth = 0.5*np.log(t/TWOPI); nn = np.arange(1, N+1, dtype=float); w = nn**-0.5
    import os; cache = path.replace('.txt', f'_zp{len(t)}.npy')
    zp = np.empty(len(t)); z0 = np.empty(len(t))
    if os.path.exists(cache):
        zp, z0 = np.load(cache)
    else:
      for a in range(0, len(t), chunk):
        ph = thm[a:a+chunk, None] - (r[None, :N] + d[a:a+chunk, None]*ln[None, :N])
        zp[a:a+chunk] = -2*np.sum(w*(dth[a:a+chunk, None] - ln[None, :N])*np.sin(ph), 1)
        z0[a:a+chunk] = 2*np.sum(w*np.cos(ph), 1)
      p = np.sqrt(t/TWOPI) - N; c0 = np.cos(2*np.pi*(p*p-p-1/16))/np.cos(2*np.pi*p)
      z0 += (-1)**(N-1)*(t/TWOPI)**-0.25*c0
      np.save(cache, np.stack([zp, z0]))
    D = -1j*np.exp(-1j*thm)*zp; out = dict(N_rs=N, T=float(T), L=float(L), nzeros=len(t),
        maxZ=float(np.max(abs(z0))), medZp=float(np.median(abs(zp))))
    mu = mobius(ymax); res = []
    for th in thetas:
        y = T**th; n = np.nonzero(mu[1:int(y)+1])[0] + 1; x = np.log(y/n)/np.log(y)
        Wt = np.stack([mu[n]*n**-0.5*x**j for j in range(1, J+1)], 1); Mp = np.empty((len(t), J), complex)
        for a in range(0, len(t), chunk):
            E = np.exp(-1j*(r[None, n-1] + d[a:a+chunk, None]*ln[None, n-1])); Mp[a:a+chunk] = E @ Wt
        for JJ in range(1, J+1):
            YA = D[:, None]*Mp[:, :JJ]
            rr = dict(theta=th, J=JJ, CGG_1D=bound_complex(YA), real_2D=bound_real(YA))
            y1 = YA[:, 0]; z = y1/(y1.mean()/abs(y1.mean()))
            rr['par'] = float(np.var(z.real)/(np.var(z.real)+np.var(z.imag))); rr['mu/L'] = float(abs(y1.mean())/L)
            rr['A/L2'] = float(np.mean(abs(z)**2)/L**2); rr['ReB/L2'] = float(np.mean(z**2).real/L**2); res.append(rr)
    out['res'] = res; print(json.dumps(out), flush=True)
if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]), [float(x) for x in sys.argv[3].split(',')], int(sys.argv[4]), int(sys.argv[5]) if len(sys.argv) > 5 else None)
