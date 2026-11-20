# ceilaudit8 B2: Laguerre orders vs off-line pairs (explore29 barrier). |f(x+iy)|^2 = sum_n L_n(f)(x) y^{2n}/(2n)!,
# L_n/f(x)^2 = (2n)! [Y^n] prod_real(1+Y/(x-g)^2) * prod_pairs (1 + 2Y(w^2-d^2)/A^2 + Y^2/A^2), Y=y^2, w=x-g0, A=w^2+d^2.
# Also L_1(f^(m)) (Turan-type: Fleischer applied to |f^(m)|) via log-derivative sums and Bell polynomials.
import numpy as np, json, sys, math
Z = np.load('Riemann/Deep/ExploreHelio/results/zeta_zeros.npy'); T1 = Z[-1]
NMAX = 40
def coeffs(x, real, pairs, tail_p):  # x: grid; real: real zeros (incl. negatives); pairs: list (g0,d); returns c[n, ix]
    near = []; far = []
    c = np.zeros((NMAX+1, x.size)); c[0] = 1.0
    xm = x.mean(); sel = np.abs(real - xm) < 400; rn = real[sel]; rf = real[~sel]
    for g in rn:
        u = 1.0/(x-g)**2
        c[1:] = c[1:] + u*c[:-1]
    # far zeros: log prod = sum_m (-1)^{m+1} p_m Y^m/m, m<=3
    p = [None]+[np.full(x.size, tail_p[m-1]) for m in (1,2,3)]
    for j in range(0, rf.size, 1000):
        u = 1.0/(x[:,None]-rf[None,j:j+1000])**2
        for m in (1,2,3): p[m] += (u**m).sum(1)
        del u
    lg = np.zeros((NMAX+1, x.size)); lg[1]=p[1]; lg[2]=-p[2]/2; lg[3]=p[3]/3
    ex = np.zeros_like(c); ex[0]=1.0  # exp of series via e' = l' e
    for n in range(1, NMAX+1):
        ex[n] = sum(k*lg[k]*ex[n-k] for k in range(1, min(n,3)+1))/n
    c = conv(c, ex)
    for (g0, d) in pairs:
        w = x-g0; A = w*w+d*d; f = np.zeros((3, x.size)); f[0]=1; f[1]=2*(w*w-d*d)/A**2; f[2]=1/A**2
        c = conv(c, f)
    return c
def conv(a, b):
    out = np.zeros((NMAX+1, a.shape[1]))
    for i in range(b.shape[0]):
        out[i:] += b[i]*a[:NMAX+1-i]
    return out
def tailsum(x0, dens, Tlo):  # int_Tlo^inf dens(g) [ (g-x0)^-2m + (g+x0)^-2m ] dg, m=1,2,3
    g = np.geomspace(Tlo, Tlo*1e4, 20001); w = np.gradient(g)
    return [float((dens(g)*((g-x0)**(-2*m)+(g+x0)**(-2*m))*w).sum()) for m in (1,2,3)]
def first_fail(c, x, gbar):
    out = []
    for n in range(1, NMAX+1):
        mn = c[n].min()
        if mn < 0: out.append(n)
    return out
def zeta_model(ia, ib, dg_list, tag):
    ga, gb = Z[ia], Z[ib]; g0 = (ga+gb)/2; gbar = 2*math.pi/math.log(g0/(2*math.pi))
    keep = np.ones(Z.size, bool); keep[[ia, ib]] = False; zz = Z[keep]; real = np.concatenate([zz, -zz])
    tp = tailsum(g0, lambda g: np.log(g/(2*math.pi))/(2*math.pi), T1)
    print(f"== {tag}: pair {ga:.4f},{gb:.4f} g0={g0:.4f} gbar={gbar:.4f}")
    for dg in dg_list:
        d = dg*gbar; L = 2*d+3*gbar; x = np.linspace(g0-L, g0+L, 3001)+1e-7
        c = coeffs(x, real, [(g0, d)], tp)
        # check: sum_n c_n(g0) d^{2n} ~ 0 (f(g0+id)=0)
        i0 = np.argmin(np.abs(x-g0)); terms = np.array([c[n,i0]*d**(2*n) for n in range(NMAX+1)])
        ff = first_fail(c, x, gbar)
        # L1 normalised like explore29: max (log|f|)'' gbar^2 = -min c1 gbar^2 on merged interval
        lo = Z[ia-1]; hi = Z[ib+1]; mi = (x>lo)&(x<hi)
        l1 = -(c[1][mi]).min()*gbar**2
        # where does the first failure occur? |x-g0|/d
        nstar = ff[0] if ff else None
        loc = (abs(x[np.argmin(c[nstar])]-g0)/d) if ff else None
        print(f"delta/gbar={dg:5.2f} L1max(log)''gbar2={l1:+9.3f} failing orders={ff[:12]} n*={nstar} pred pi*d/(2gbar)={math.pi*dg/2:.2f} argmin|x-g0|/d={loc if loc is None else round(loc,3)} check sum c_n d^2n / sum|.| = {terms.sum()/np.abs(terms).sum():.2e}")
        sys.stdout.flush()
def dh():
    d = json.load(open('Riemann/Deep/ExploreHelio/results/dh_zeros.json')); on = np.array(d['online']); off = d['offline']; TM = d['TMAX']
    real = np.concatenate([on, -on]); pairs = [(g, b-0.5) for (b, g) in off] + [(-g, b-0.5) for (b, g) in off]
    print("== DH: 9 off-line pairs below 420")
    for (b, g) in off:
        gbar = 2*math.pi/math.log(5*g/(2*math.pi)); dd = b-0.5
        tp = tailsum(g, lambda t: np.log(5*t/(2*math.pi))/(2*math.pi), TM)
        L = 2*dd+3*gbar; x = np.linspace(g-L, g+L, 3001)+1e-7
        c = coeffs(x, real, pairs, tp); ff = first_fail(c, x, gbar)
        print(f"DH pair g={g:.3f} delta={dd:.4f} delta/gbar={dd/gbar:.3f} failing orders={ff[:8]} min c1*gbar^2={c[1].min()*gbar**2:+.3f}")
if __name__ == '__main__':
    zeta_model(6708, 6709, [0.2, 0.5, 0.75, 1.0, 1.1, 1.25, 1.5, 2.0, 3.0, 4.0, 6.0], 'Lehmer #6709/#6710')
    zeta_model(6606, 6607, [0.5, 1.0, 1.25, 1.5, 2.0, 3.0, 4.0, 6.0], 'generic #6607/#6608')
    dh()
