import os
BASE = os.path.abspath(os.environ.get('RAMAUDIT', '.scratch/ramaudit')) + '/'
import sys, os, time, json, math, pickle, glob
R=BASE + 'repo/ramanujans-inequality-HEAD'
sys.path.insert(0, R+'/certify/arb'); sys.path.insert(0, R+'/certify')
os.chdir(R+'/certify/arb')
import cells_arb as CA
import scan_arb as SA
from flint import arb
PI = arb.pi()
OUT = BASE + 'cellmod_%s.json' % sys.argv[1]
LKS = [float(x) if '.' in x else int(x) for x in sys.argv[2].split(',')]
# rigorous crude bound for sup|G(1/2+it)| <= 16 * PROD * W / log^2 W, PROD = prod_p (1 + (sqrt p - 1)/p^2)
def prodG():
    import math
    N = 10**6
    sieve = bytearray([1])*(N+1); sieve[0]=sieve[1]=0
    for i in range(2,int(N**0.5)+1):
        if sieve[i]: sieve[i*i::i] = bytearray(len(sieve[i*i::i]))
    lp = 0.0
    for p in range(2,N+1):
        if sieve[p]: lp += math.log1p((math.sqrt(p)-1)/p/p)
    # tail: sum_{p>N} (sqrt p)/p^2 <= sum_{n>N} n^{-1.5} <= 2/sqrt(N)
    return math.exp(lp + 2/math.sqrt(N))
PG = prodG()
print('PROD', PG, flush=True)
tab = {}
for f in sorted(glob.glob(os.path.join(CA.CODE, 'lemmaD_tableGH0.01_*.pkl'))): tab.update(pickle.load(open(f, 'rb')))
try: out = json.load(open(OUT))
except Exception: out = {}
t0 = time.time(); n = 0
for key in sorted(k for k in tab if tab[k] is not None and k[1] in LKS):
    k, lk = key[0], key[1]; nm = '%d,%s' % (k, lk)
    if nm in out: continue
    fv, p = tab[key]
    de = arb(k)/10**4; L = arb(lk)
    dd, den, d4 = CA.consts_arb(p, de, L)
    u, t, w, x = [arb(float(v)) for v in p]; al = 1 - de
    d5 = arb('1.0061')*(1 + (arb('1.333') + 4*(-w*L/2).exp())/(w*L))
    J = arb('4.72')*arb('0.618')*(L/6).exp()*L
    c5 = d4*d5/(2*PI)*((w*L).exp()/(w*L))*((L*x*(arb('1.5') - 2*al)).exp() + ((L*u - 2*L.log())*(arb('1.5') - 2*al)).exp())*J
    c1sq = den + c5
    rG = arb(16*PG)/(d5*(w*L))          # crude/Bellotti ratio for the G(1/2+it) factor
    rG = rG if rG > 1 else arb(1)
    den_rep = c1sq - c5*rG
    pair, nt = CA.pair_arb(k, float(lk))
    far = 2*(2*PI).sqrt()*arb('2.4')*CA.Iprime_arb(k, lk)
    nu = SA.lam_at(float(lk))/L
    YX = ((u - x)*L - 2*L.log())          # log(Y/X), Y = U/L^2
    kap2 = (1 + (YX*2*nu).exp())/2
    base = 2*arb(float(p[3]))*de*L - (de*L).log()
    v_old = (dd*(pair + far)*nt/den).log() + base
    v_far2 = (dd*(pair + 2*far)*nt/den).log() + base
    v_rep = (dd*kap2*(pair + 2*far)*nt/den_rep).log() + base if den_rep > arb('0.2') else None
    out[nm] = {'old': repr(SA.fup(v_old)), 'far2': repr(SA.fup(v_far2)), 'rep': (repr(SA.fup(v_rep)) if v_rep is not None else None),
               'pair': float(pair.mid()), 'far': float(far.mid()), 'kap2': float(kap2.mid()), 'rG': float(rG.mid()), 'c5': float(c5.mid()), 'den': float(den.mid())}
    n += 1
    if n % 20 == 0:
        json.dump(out, open(OUT + '.tmp', 'w')); os.replace(OUT + '.tmp', OUT)
        print(n, nm, out[nm], '%.0fs' % (time.time() - t0), flush=True)
json.dump(out, open(OUT + '.tmp', 'w')); os.replace(OUT + '.tmp', OUT)
print('done', n, '%.0fs' % (time.time() - t0), flush=True)

