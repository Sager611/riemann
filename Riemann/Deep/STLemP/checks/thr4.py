
import numpy as np, sys, signal
signal.alarm(175)

NU = np.arange(0.0, 8.0, 0.0025)

def E_sieves(eta, m, v, opts):
    nu = NU
    vals = [1 + np.maximum(nu, 0.5) + nu*(2*eta-1)]          # per-D (LS1)
    vals.append(np.maximum(nu, 3.5) + nu*(2*eta-1))          # HT (LS2)
    if opts.get('PT', True):
        vals.append(np.maximum(nu+0.75, 2.5) + nu*(2*eta-1)) # PT (LS3)
    if opts.get('lift', False):
        best = np.full_like(nu, np.inf)
        for k in range(0, 41):
            x = np.maximum(1.0, nu) + k*0.05
            E0 = np.maximum(nu - x/2, 0)
            if m <= 1:
                main = 0.5*x + v
            else:
                crude = v + (m-1)/2*E0
                refined = v*2/(m+1) + nu*(m-1)/(m+1)
                if opts.get('crude', False): refined = 1e9
                main = 0.5*x + np.minimum(crude, refined)
            fam = np.maximum(x, main)
            val = 0.5*x + fam + nu*(2*eta-1)
            best = np.minimum(best, val)
        vals.append(best)
    E = np.min(np.vstack(vals), axis=0)
    E = np.where(nu >= v - 1e-12, E, -np.inf)   # class I lengths start at V
    return E

def F_H(eta, m, i2, opts):
    best = 1e9
    vs = [0.5] if not opts.get('varyV', False) else list(np.arange(0.05, 0.5001, 0.025))
    for v in vs:
        E = E_sieves(eta, m, v, opts)
        runmax = np.maximum.accumulate(E)
        c2 = (i2 + 1.5)/2 - NU*(0.5 - eta)
        val = np.where(NU >= v, np.maximum(runmax, c2), np.inf)
        best = min(best, float(np.min(val)))
    return best

def F_F(eta, i2):
    return (7 + 2*i2)*eta/(1+2*eta)

def ok(A, m, i2, opts, useC):
    etamax = 1/(2*A)
    for k in range(1, 201):
        eta = etamax*k/200
        bar = 1 + A*eta
        if useC and A > (m+2)/(1+2*eta): continue
        if F_F(eta, i2) < bar: continue
        if F_H(eta, m, i2, opts) < bar: continue
        return False
    return True

def thr(m, i2, opts, useC):
    lo, hi = 1.0, 6.0
    for _ in range(22):
        mid = (lo+hi)/2
        if ok(mid, m, i2, opts, useC): hi = mid
        else: lo = mid
    return hi

cfg = sys.argv[1]; m = int(sys.argv[2])
i2 = min((m+3)/2, 3*(m+1)/4)
O = {'lift':True, 'crude':('crude' in cfg), 'varyV':('varyV' in cfg), 'PT':True}
print(cfg, m, round(thr(m, i2, O, m==1), 4), flush=True)
