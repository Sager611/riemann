# Feasibility prototype (float, NOT rigorous): point estimates of the per-block bound at one t, combined vs split square root.
# uv run --no-project --with mpmath --with numpy python3 feasibility_proto.py 540 700 875   (h = 1.01, 1.018, 1.03, and split)
# float prototype (log-space) of per-block bounds; NOT rigorous
import numpy as np, mpmath as mp, sys
from sanity_paper import AkBk
pi = np.pi
b1 = float((1/mp.pi)*(53632*mp.sqrt(2)*mp.pi**4/729)**(mp.mpf(1)/30)*(76545*mp.sqrt(2)/107264)**(mp.mpf(1)/8))
b2 = float(mp.pi**(-mp.mpf(7)/8)*(729/(53632*mp.sqrt(2)*mp.pi**4))**(mp.mpf(1)/30))
ETAS = np.linspace(0.8, 2.6, 19)
_cache = {}
def A5B5(eta, h):
    key = (round(eta, 6), round(h, 9))
    if key not in _cache:
        h5 = 76545*np.sqrt(2)/107264*h**9
        A5, B5 = AkBk(5, mp.mpf(eta), mp.mpf(h5))
        _cache[key] = (float(A5), float(B5))
    return _cache[key]
_c4 = {}
def A4B4(eta, H):
    key = (round(eta, 6), round(H, 9))
    if key not in _c4:
        A4, B4 = AkBk(4, mp.mpf(eta), mp.mpf(H))
        _c4[key] = (float(A4), float(B4))
    return _c4[key]

def blk_aba3b(la, lt, h, etas=ETAS, combined=True, qopt=True, th1=None):
    """log of bound on |sum_{a<n<=b} n^{-1/2-it}| for block a<n<=ha. returns min over eta, q."""
    best = np.inf
    lN = la + np.log(h - 1)
    E6 = 6 + 5/pi*np.log(2)
    E4 = 4.686*np.sqrt(pi)*h**1.5
    E5 = 0.893*np.exp(-la/3)*5/pi*((h-1)/pi)**(1/3) + 6/pi*h**3*(h-1)
    for eta in etas:
        A5, B5 = A5B5(eta, h)
        c3 = b1*np.sqrt(pi)*A5*h**(21/8)*(h-1)
        c4 = b2*np.sqrt(pi)*B5*h**1.5*(h-1)**(7/8)
        # q grid around th1 a^{36/41} t^{-11/41}
        lq0 = 36/41*la - 11/41*lt
        if th1 is not None:
            lqs = np.array([lq0 + np.log(th1)])
        else:
            lqs = lq0 + np.log(np.linspace(0.2, 3.0, 57))
        lqs = np.maximum(lqs, 0.0)
        # second factor terms (logs)
        T = np.stack([lN - lqs,
                      np.log(1800/2911*c3) - la/5 + 11/30*(lqs + lt),
                      np.log(28800/54481*c4) - 11/20*la + 61/120*(lqs + lt),
                      np.log(8/3*E4) + 1.5*la - 0.5*(lt + lqs),
                      np.log(9/14*E5) + (lt + lqs)/3 - la/3,
                      np.full_like(lqs, np.log(E6))])
        lfirst = np.logaddexp(lN, lqs)  # N - 1 + q <= N + q
        if combined:
            ltot = lfirst + np.logaddexp.reduce(T, axis=0)
            lb = 0.5*ltot
        else:
            # paper style: sqrt split into (T0+T1), T2, T3, T4, T5
            lb = np.logaddexp.reduce(np.stack([0.5*(lfirst + np.logaddexp(T[0], T[1]))] + [0.5*(lfirst + T[i]) for i in range(2, 6)]), axis=0)
        m = lb.min()
        if m < best: best = m
    return best - la/2

def blk_4th(la, lt, h, etas=np.linspace(0.5, 2.0, 16)):
    best = np.inf
    lN = la + np.log(h - 1)
    llam = np.log(3/pi) + lt - 4*(la + np.log(h))
    for eta in etas:
        A4, B4 = A4B4(eta, h**4)
        v = np.logaddexp(np.log(A4*h) + lN + llam/14, np.log(B4) + 0.75*lN - llam/14)
        best = min(best, v)
    return best - la/2

def total(lt, h=1.018, lamin_frac=0.30, combined=True, qopt=True, use4=True, verbose=False):
    la_top = 0.5*(lt - np.log(2*pi))
    la_min = lamin_frac*lt
    nb = int(np.ceil((la_top - la_min)/np.log(h)))
    s = 0.0; parts = {'aba3b': 0.0, '4th': 0.0, 'triv': 0.0}
    for k in range(1, nb + 1):
        la = la_top - k*np.log(h)
        v1 = blk_aba3b(la, lt, h, combined=combined) if la > 0.36*lt else np.inf
        v2 = blk_4th(la, lt, h) if use4 else np.inf
        vt = np.log(np.sqrt(h) - 1) + la/2 + np.log(2)  # 2(sqrt(b)-sqrt(a))
        v = min(v1, v2, vt)
        key = 'aba3b' if v == v1 else ('4th' if v == v2 else 'triv')
        e = np.exp(v - 27/164*lt)
        parts[key] += e; s += e
    s1 = 2*np.exp(la_min/2 - 27/164*lt)
    s += s1; parts['S1'] = s1
    A = 2*s
    if verbose: print({k: round(2*x, 4) for k, x in parts.items()})
    return A


import numpy as np, mpmath as mp, sys
pi = np.pi
def blocks_aba3b(la, lt, h, etas, combined=True, qgrid=np.linspace(0.15, 3.0, 58)):
    la = la[:, None]
    lN = la + np.log(h - 1)
    E6 = 6 + 5/pi*np.log(2); E4 = 4.686*np.sqrt(pi)*h**1.5
    E5 = 0.893*np.exp(-la/3)*5/pi*((h-1)/pi)**(1/3) + 6/pi*h**3*(h-1)
    lq = np.maximum(36/41*la - 11/41*lt + np.log(qgrid)[None, :], 0.0)
    lfirst = np.logaddexp(lN, lq)
    best = np.full(la.shape[0], np.inf)
    for eta in etas:
        A5, B5 = A5B5(eta, h)
        c3 = b1*np.sqrt(pi)*A5*h**(21/8)*(h-1); c4 = b2*np.sqrt(pi)*B5*h**1.5*(h-1)**(7/8)
        T = [lN - lq, np.log(1800/2911*c3) - la/5 + 11/30*(lq + lt), np.log(28800/54481*c4) - 11/20*la + 61/120*(lq + lt),
             np.log(8/3*E4) + 1.5*la - 0.5*(lt + lq), np.log(9/14*E5) + (lt + lq)/3 - la/3, np.full_like(lq, np.log(E6))]
        if combined:
            lb = 0.5*(lfirst + np.logaddexp.reduce(np.stack(T), axis=0))
        else:
            lb = np.logaddexp.reduce(np.stack([0.5*(lfirst + np.logaddexp(T[0], T[1]))] + [0.5*(lfirst + T[i]) for i in range(2, 6)]), axis=0)
        best = np.minimum(best, lb.min(axis=1))
    return best - la[:, 0]/2

def blocks_4th(la, lt, h, etas=np.linspace(0.5, 2.0, 16)):
    lN = la + np.log(h - 1); llam = np.log(3/pi) + lt - 4*(la + np.log(h))
    best = np.full(la.shape, np.inf)
    for eta in etas:
        A4, B4 = A4B4(eta, h**4)
        best = np.minimum(best, np.logaddexp(np.log(A4*h) + lN + llam/14, np.log(B4) + 0.75*lN - llam/14))
    return best - la/2

def total(lt, h=1.018, combined=True, etas=np.linspace(0.8, 2.6, 19), verbose=False, lamin_frac=0.30):
    la_top = 0.5*(lt - np.log(2*pi)); la_min = lamin_frac*lt
    nb = int(np.ceil((la_top - la_min)/np.log(h)))
    la = la_top - np.arange(1, nb + 1)*np.log(h)
    v1 = np.where(la > 0.36*lt, blocks_aba3b(la, lt, h, etas, combined), np.inf)
    v2 = blocks_4th(la, lt, h)
    vt = np.log(np.sqrt(h) - 1) + la/2 + np.log(2)
    v = np.minimum(np.minimum(v1, v2), vt)
    e = np.exp(v - 27/164*lt)
    s1 = 2*np.exp(la_min/2 - 27/164*lt)
    A = 2*(e.sum() + s1)
    if verbose:
        p = {'aba3b': 2*e[(v == v1)].sum(), '4th': 2*e[(v == v2) & (v != v1)].sum(), 'triv': 2*e[(v == vt) & (v != v1) & (v != v2)].sum(), 'S1': 2*s1}
        print(lt, {k: round(float(x), 3) for k, x in p.items()})
    return A

if __name__ == '__main__':
    hs = [1.018]
    for lt in [float(x) for x in sys.argv[1:]]:
        out = []
        for h in [1.01, 1.018, 1.03]:
            out.append((h, round(float(total(lt, h=h)), 3)))
        print(lt, 'combined', out, 'split(h=1.018)', round(float(total(lt, combined=False)), 3)); sys.stdout.flush()
