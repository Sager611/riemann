#!/usr/bin/env python3
# dbnrecaudit falsify_dcoef.py: imported input I4 (U.tex Lemma src:coeff l.1266-1272, l.3311-3320): d_t(n) >= 0 and nondecreasing in t,
# 2 <= n <= 4096, on a dense t grid of [0, 0.11] (float64, all n at once); and the derivative identity
# d/dt d_t(n) = (log^2 n/4) d_t(n) + (log n/4) sum_{ab=n, a,b>=2} d_t(a) d_t(b) against central differences. Falsification only.
import numpy as np, json
K = 4096; ts = np.linspace(0.0, 0.11, 2201)
ln = np.log(np.arange(K + 1, dtype=float)); ln[0] = ln[1] = 0.0
def dco(t):
    b = np.exp(np.outer(ln**2/4, t)); b[0] = 0; b[1] = 1
    d = ln[:, None]*b
    for r in range(2, K//2 + 1):
        m = np.arange(2, K//r + 1); d[r*m] -= d[r][None, :]*b[m]
    return d
neg, negd = np.inf, np.inf
for c in range(0, len(ts) - 1, 200):   # chunks of 201 t values (overlapping by one), < 30 MB per array
    tc = ts[c:c + 201]; d = dco(tc)
    scale = np.log(np.arange(2, K + 1))[:, None]*np.exp(np.outer(ln[2:]**2/4, tc))
    neg = min(neg, float((d[2:]/scale).min()))
    dd = np.diff(d[2:], axis=1)/np.diff(tc)[None, :]
    negd = min(negd, float((dd/scale[:, 1:]).min()))
    d0col = d[:, 0] if c == 0 else d0col
# derivative identity at a few t
h = 1e-6; res = []
for t0 in (0.0025, 0.03, 0.1):
    dp, dmn, d0 = dco(np.array([t0 + h]))[:, 0], dco(np.array([t0 - h]))[:, 0], dco(np.array([t0]))[:, 0]
    fd = (dp - dmn)/(2*h)
    rhs = ln**2/4*d0
    conv = np.zeros(K + 1)
    for a in range(2, K//2 + 1):
        bb = np.arange(2, K//a + 1); conv[a*bb] += d0[a]*d0[bb]
    rhs = rhs + ln/4*conv
    rel = np.abs(fd[2:] - rhs[2:])/(np.abs(rhs[2:]) + 1e-12)
    res.append(dict(t=t0, max_rel_err=float(rel.max())))
nonpp = [n for n in range(2, K + 1) if len({p for p in range(2, n + 1) if n % p == 0 and all(p % q for q in range(2, int(p**0.5) + 1))}) > 1]
print(json.dumps(dict(n_max=K, t_grid=[0.0, 0.11, len(ts)], min_d_over_scale=float(neg), min_dtd_over_scale=float(negd),
                      derivative_identity=res, d_at_0_nonprimepower_maxabs=float(np.abs(d0col[nonpp]).max()))))

