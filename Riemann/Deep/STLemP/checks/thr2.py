
import numpy as np, sys, signal
signal.alarm(170)

NU = np.arange(0.5, 12.0, 0.0025)

def E_sieves(eta, m, opts):
    nu = NU
    vals = [1 + np.maximum(nu, 0.5) + nu*(2*eta-1)]
    if opts.get('HT', True):
        vals.append(np.maximum(nu, 3.5) + nu*(2*eta-1))
    if opts.get('PT', True):
        vals.append(np.maximum(nu+0.75, 2.5) + nu*(2*eta-1))
    if opts.get('lemP', False):
        # Lemma P family sieve for nu <= 1 (m=1 detector: main fine); for m>=2 crude main loss
        if m <= 1 or opts.get('mainfree', False):
            v = 1.5 + nu*(2*eta-1)
        else:
            E0 = np.maximum(nu - 0.5, 0)
            main = 0.5 + 0.5 + np.minimum((m-1)/2*E0, -0.5 + 0.5*2/(m+1) + nu*(m-1)/(m+1) )
            v = 0.5 + np.maximum(1.0, main) + nu*(2*eta-1)
        v = np.where(nu <= 1.0, v, np.inf)
        vals.append(v)
    if opts.get('lift', False):
        best = np.full_like(nu, np.inf)
        for k in range(0, 61):
            x = np.maximum(1.0, nu) + k*0.05
            E0 = np.maximum(nu - x/2, 0)
            if m <= 1 or opts.get('mainfree', False):
                main = 0.5*x + 0.5
            else:
                crude = 0.5 + (m-1)/2*E0
                refined = 0.5*2/(m+1) + nu*(m-1)/(m+1)
                main = 0.5*x + np.minimum(crude, refined)
            fam = np.maximum(x, main)
            val = 0.5*x + fam + nu*(2*eta-1)
            best = np.minimum(best, val)
        vals.append(best)
    return np.min(np.vstack(vals), axis=0)

def F_H(eta, m, i2, opts):
    E = E_sieves(eta, m, opts)
    runmax = np.maximum.accumulate(E)
    c2 = (i2 + 1.5)/2 - NU*(0.5 - eta)
    return np.min(np.maximum(runmax, c2))

def F_F(eta, i2):
    return (7 + 2*i2)*eta/(1+2*eta)

def ok(A, m, i2, opts, useC):
    etamax = 1/(2*A)
    for k in range(1, 301):
        eta = etamax*k/300
        bar = 1 + A*eta
        if useC and A > (m+2)/(1+2*eta): continue
        if F_F(eta, i2) < bar: continue
        if F_H(eta, m, i2, opts) < bar: continue
        return False
    return True

def thr(m, i2, opts, useC):
    lo, hi = 1.0, 20.0
    for _ in range(30):
        mid = (lo+hi)/2
        if ok(mid, m, i2, opts, useC): hi = mid
        else: lo = mid
    return hi

cfg = sys.argv[1]
ms = [int(a) for a in sys.argv[2].split(',')]
for m in ms:
    i2_perD = (m+3)/2
    i2_lift = min(i2_perD, 3*(m+1)/4)
    if m == 1 and 'lemP' in cfg:
        i2_lift = 1.5
    O = {'HT':True,'PT':('noPT' not in cfg),'lift':('lift' in cfg and 'liftI2only' not in cfg),'lemP':('lemP' in cfg),'mainfree':('mainfree' in cfg)}
    i2 = i2_lift if ('lift' in cfg or 'lemP' in cfg) else i2_perD
    t = thr(m, i2, O, useC=(m==1))
    print(cfg, m, round(t,4), flush=True)
