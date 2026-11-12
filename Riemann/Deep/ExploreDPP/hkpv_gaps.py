
# HKPV/DPP side: gap laws of the sine process and the a-discrete sine processes vs real zeta and DH zeros.
import numpy as np, json, sys, resource
from numpy.polynomial.legendre import leggauss
LAMS = [0.5, 7/11, 0.76, 0.7851, 1.0]
def sine_E0(s, m=40):
    x, w = leggauss(m); x = 0.5*s*(x+1); w = 0.5*s*w
    d = x[:,None]-x[None,:]
    K = np.sinc(d)                     # np.sinc(t) = sin(pi t)/(pi t)
    A = np.sqrt(w)[:,None]*K*np.sqrt(w)[None,:]
    lam = np.linalg.eigvalsh(A)
    return np.prod(1-lam), np.sort(lam)[::-1]
def sine_G(l, e=1e-4):
    # spacing CDF G(l) = 1 + dE(0;l)/dl
    return 1 + (sine_E0(l+e)[0]-sine_E0(l-e)[0])/(2*e)
def dsine_gap(a, mmax=12):
    # a-discrete sine on aZ, K_ij = a*sinc(a(i-j)); Palm law of the gap: P(gap = m a)
    out = {}
    for m in range(1, mmax+1):
        idx = np.arange(m+1); K = a*np.sinc(a*(idx[:,None]-idx[None,:]))
        Bm = np.zeros(m+1); Bm[1:m] = 1.0
        out[m] = abs(np.linalg.det(K - np.diag(Bm)))/a
    return out
def D_from_gaps(g, l): return float(np.mean(g <= l))
def parse_gw(fn):
    gaps = []; nz = 0; tcs = []
    lines = open(fn).read().split('\n')
    tc = None
    for L in lines:
        if L.startswith('ZG'):
            tc = float(L.split()[2].split('=')[1]); tcs.append(tc)
        elif L.startswith('ZZ'):
            z = np.array([float(v) for v in L.split()[2:]])
            nz += len(z)
            dens = np.log(tc/(2*np.pi))/(2*np.pi)
            gaps.append(np.diff(np.sort(z))*dens)
    return np.concatenate(gaps), nz, np.mean(tcs)
rows = []
# sine process
Gs = {l: sine_G(l) for l in LAMS}
rows.append(('sine process (GUE limit), HKPV/Fredholm', {l: Gs[l] for l in LAMS}, None))
for a in (0.5, 7/11):
    pg = dsine_gap(a); tot = sum(pg.values())
    rows.append((f'{a:.4f}-discrete sine (DPP on aZ), Palm gap law', {l: sum(p for m,p in pg.items() if m*a <= l+1e-12) for l in LAMS}, f'sum P(gap)={tot:.6f}; P(gap=a)={pg[1]:.6f}, P(gap=2a)={pg[2]:.6f}'))
# zeta low height
z = np.load(sys.argv[1]); z = z[z > 100]
g = np.diff(z)*np.log(z[:-1]/(2*np.pi))/(2*np.pi)
rows.append((f'zeta zeros 100<t<{z[-1]:.0f} ({len(g)} gaps)', {l: D_from_gaps(g, l) for l in LAMS}, f'mean gap {g.mean():.4f}'))
for fn in sys.argv[2:4]:
    gg, nz, tcm = parse_gw(fn)
    rows.append((f'zeta zeros t~{tcm:.3g} (gram_win -zeros, {len(gg)} gaps)', {l: D_from_gaps(gg, l) for l in LAMS}, f'mean gap {gg.mean():.4f}'))
# DH zeros on the line
d = json.load(open(sys.argv[4]))
on = np.sort(np.array(d['online'], dtype=float)); on = on[on > 20]
gd = np.diff(on)*np.log(on[:-1]/(2*np.pi))/(2*np.pi)   # same smooth density as zeta (DH: conductor 5, use log(5t/2pi)?)
gd5 = np.diff(on)*np.log(5*on[:-1]/(2*np.pi))/(2*np.pi)
rows.append((f'DH zeros on the line 20<t<{on[-1]:.0f} ({len(gd)} gaps), density log(5t/2pi)/2pi', {l: D_from_gaps(gd5, l) for l in LAMS}, f'mean gap {gd5.mean():.4f}; offline zeros in file: {len(d["offline"])}'))
print('| configuration | ' + ' | '.join(f'D({l:.4g})' for l in LAMS) + ' | note |')
print('|---|' + '---|'*len(LAMS) + '---|')
for name, vals, note in rows:
    print(f'| {name} | ' + ' | '.join(f'{vals[l]:.4f}' for l in LAMS) + f' | {note or ""} |')
E0, lamk = sine_E0(0.76)
print('HKPV: sine kernel on [0,0.76]: Bernoulli parameters', np.round(lamk[:4], 6), 'E(0;0.76)=', round(E0, 6), 'P(N>=2)=', round(1 - np.prod(1-lamk) - sum(lamk[i]*np.prod(np.delete(1-lamk, i)) for i in range(len(lamk))), 8))
print('maxrss MB', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6)
