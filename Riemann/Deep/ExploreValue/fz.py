# F and F' for F = zeta or Davenport-Heilbronn (DH), vectorised numpy; chi'/chi; Hardy-type Z.
# Dirichlet series to M terms + Euler-Maclaurin tail (K Bernoulli terms), as in ExploreSymm/code/levin.py.
import numpy as np, math
from scipy.special import bernoulli, digamma
K = 12; B2 = bernoulli(2*K)[2::2]; FAC = [math.factorial(2*k) for k in range(1, K+1)]
KAP = (math.sqrt(10-2*math.sqrt(5)) - 2)/(math.sqrt(5) - 1)
class Fn:
    def __init__(self, kind, tmax):
        self.kind = kind; self.N = int(tmax/math.pi) + 60
        if kind == 'zeta':
            self.q = 1; self.classes = [(1.0, 0.0)]; M = self.N - 1; coef = np.ones(M)
        else:
            self.q = 5; per = [1.0, KAP, -KAP, -1.0, 0.0]
            self.classes = [(per[c-1], c/5) for c in range(1, 5)]; M = 5*self.N - 1
            coef = np.array([per[(m-1) % 5] for m in range(1, M+1)])
        self.coef = coef; self.logm = np.log(np.arange(1, M+1, dtype=float))
    def _htail(self, s, X):
        lX = math.log(X); Xs = np.exp(-s*lX)
        h = X*Xs/(s-1) + Xs/2; hs = -lX*X*Xs/(s-1) - X*Xs/(s-1)**2 - lX*Xs/2
        poch = s.copy(); dpoch = np.ones_like(s)
        for k in range(1, K+1):
            if k > 1:
                for i in (2*k-3, 2*k-2):
                    dpoch = dpoch*(s+i) + poch; poch = poch*(s+i)
            term = B2[k-1]/FAC[k-1]*X**(-(2*k-1))*Xs
            h = h + term*poch; hs = hs + term*(dpoch - lX*poch)
        return h, hs
    def __call__(self, s, blk=200):
        s = np.atleast_1d(np.asarray(s, complex)); F = np.empty_like(s); Fp = np.empty_like(s)
        for i in range(0, len(s), blk):
            sb = s[i:i+blk]; E = np.exp(-np.outer(sb, self.logm))
            F[i:i+blk] = E @ self.coef; Fp[i:i+blk] = E @ (-self.coef*self.logm)
            for w, al in self.classes:
                h, hs = self._htail(sb, self.N + al)
                if self.q > 1:
                    qs = np.exp(-sb*math.log(self.q)); h, hs = qs*h, qs*(hs - math.log(self.q)*h)
                F[i:i+blk] += w*h; Fp[i:i+blk] += w*hs
        return F, Fp
    def chilog(self, s):   # chi'/chi, F(s) = chi(s) F(1-s)
        s = np.asarray(s, complex)
        if self.kind == 'zeta':
            return math.log(math.pi) - 0.5*digamma((1-s)/2) - 0.5*digamma(s/2)
        return -math.log(5/math.pi) - 0.5*digamma(1 - s/2) - 0.5*digamma((1+s)/2)
    def theta(self, t):    # chi(1/2+it) = exp(-2 i theta(t)); Z(t) = exp(i theta) F(1/2+it) real
        from scipy.special import loggamma
        t = np.asarray(t, float)
        if self.kind == 'zeta':
            return np.imag(loggamma(0.25 + 0.5j*t)) - 0.5*t*math.log(math.pi)
        return np.imag(loggamma(0.75 + 0.5j*t)) + 0.5*t*math.log(5/math.pi)
def winding(f, pts):
    """pts: closed polygon samples (complex array, counterclockwise); f values there -> (winding, max phase step)."""
    ph = np.unwrap(np.angle(f)); st = np.max(np.abs(np.diff(ph)))
    return (ph[-1] - ph[0])/(2*math.pi), st
def rect(sa, sb, t1, t2, dt_line, dt_right=0.02, ds=0.002):
    """counterclockwise boundary of [sa,sb]x[t1,t2]: bottom, right, top, left."""
    b = np.arange(sa, sb, ds) + 1j*t1
    r = sb + 1j*np.arange(t1, t2, dt_right)
    tp = np.arange(sb, sa, -ds) + 1j*t2
    lf = sa + 1j*np.arange(t2, t1, -dt_line)
    return np.concatenate([b, r, tp, lf, [sa + 1j*t1]])
def line_eval(f, sigma, t0, n, dt, Bk=100):
    """F, F' at sigma + i(t0 + dt*k), k < n, via block phase rotation (explore45's scheme)."""
    lm = f.logm; W = f.coef*np.exp(-sigma*lm); Wp = -f.coef*lm*np.exp(-sigma*lm)
    Om = np.exp(-1j*dt*np.outer(lm, np.arange(Bk))); F = np.empty(n, complex); Fp = np.empty(n, complex)
    for j0 in range(0, n, Bk):
        nb = min(Bk, n - j0); tb = t0 + dt*(j0 + np.arange(nb))
        ph = np.exp(-1j*tb[0]*lm)[:, None]*Om[:, :nb]
        F[j0:j0+nb] = W @ ph; Fp[j0:j0+nb] = Wp @ ph
        s = sigma + 1j*tb
        for w, al in f.classes:
            h, hs = f._htail(s, f.N + al)
            if f.q > 1:
                qs = np.exp(-s*math.log(f.q)); h, hs = qs*h, qs*(hs - math.log(f.q)*h)
            F[j0:j0+nb] += w*h; Fp[j0:j0+nb] += w*hs
    return F, Fp
