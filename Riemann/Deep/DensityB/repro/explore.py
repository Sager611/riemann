
# densityb explore.py: numeric exploration of zero-density exponents with/without the half-plane inputs.
# Own code (not ANTEDB code); reads best_beta.json produced by best_beta.py (ANTEDB beta bounds).
import json, math, sys, os
from fractions import Fraction as F
B = 0.87338196661320752
beta = [(float(F(a)),float(F(b)),float(F(c)),float(F(m))) for a,b,c,m in json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'best_beta.json')))]
def beta_low(al):  # al in (0,1/2]
    best = math.inf
    for x0,x1,c,m in beta:
        if x0-1e-12 <= al <= x1+1e-12: best = min(best, c+m*al)
    return best
def g_old(tau):  # exponent g with sum_{n~N} n^{-it} << N^{g}, |t| = N^tau  (tau>=1)
    al = 1.0/tau
    if al <= 0.5: return tau*beta_low(al)
    return 1 - tau/2 + tau*beta_low(1-al)
def g_new(tau):
    v = g_old(tau)
    if 1 <= tau <= 2: v = min(v, 1-B+(B-0.5)*tau)
    if tau >= 2: v = min(v, B)
    return v
def mc_range_hm(s, g):  # sup tau s.t. sup_{1<=t<=tau} g(t) < 2s-1
    lvl = 2*s-1; t = 1.0; step = 0.0005
    if g(1.0) >= lvl: return 1.0
    while t < 6 and g(t+step) < lvl: t += step
    return t
def tauB(s): return (2*s-2+B)/(B-0.5)
def rows_max(rows, s, t): return max(c+a*s+b*t for c,a,b in rows)
LIT = []
LIT.append(('L2', lambda s,t: max(2-2*s, 1+t-2*s)))
LIT.append(('Huxley', lambda s,t: max(2-2*s, 4+t-6*s)))
LIT.append(('HB', lambda s,t: max(2-2*s, 10+t-13*s)))
for k in range(1,20):
    LIT.append((f'Jutila{k}', (lambda k: lambda s,t: max(2-2*s, 4-2/k-(6-2/k)*s+t, 6*k-8*k*s+t))(k)))
LIT.append(('GM', lambda s,t: max(2-2*s, 3.6-4*s, 2.4-4*s+t)))
for k in range(1,20):
    rows = [(2,-2,0),(5,-6,0),(4*k/(k+1),-6*k/(k+1),k/(k+1)),(20*k/(4*k+3),-24*k/(4*k+3),2/(4*k+3)),(2,-4,4/3),(3,-4,0.5),(4.5,-7,1),(72/19,-112/19,18/19)]
    LIT.append((f'GM2_{k}', (lambda rows: lambda s,t: rows_max(rows,s,t))(rows)))
for k in range(1,11):
    LIT.append((f'HB2_{k}', (lambda k: lambda s,t: max(2-2*s, k*t+k*(2-4*s), 2*t/3+k*(12-16*s)/3) if s>0.75 else math.inf)(k)))
LIT.append(('Ivic2', lambda s,t: max(2-2*s, t+9-12*s, t-(84*s-65)/6)))
# Bourgain optimized LV (ANTEDB literature.py polytopes): rows [b,a_s,a_t,a_r] meaning b+a_s s+a_t t+a_r r>=0
BP = [
 [[16/3,-20/3,1/3,-1],[0,0,0,1],[10,-14,1,0],[-1,0,1,0],[4,4,-5,0],[-11,16,-1,0]],
 [[5,-7,3/4,-1],[0,0,0,1],[8,-8,-1,0],[-16,20,1/3,0],[-6,10,-7/6,0],[-4,-4,5,0]],
 [[3,-5,1,-1],[0,0,0,1],[-8,8,1,0],[2,-6,2,0],[-10,14,-2/3,0],[6,-2,-2,0]],
 [[0,-4,2,-1],[0,0,0,1],[-6,2,2,0],[-12,12,1,0],[1e6,0,-1,0],[1,-1,0,0]],
 [[8,-12,4/3,0],[0,0,0,1],[-0.5,1,0,0],[15,-21,1,0],[12,-12,-1,0],[-1.5,0,1,0],[6,-10,7/6,0]],
 [[2,-2,0,-1],[0,0,0,1],[1,-1,0,0],[-0.5,1,0,0],[-10,14,-1,0],[-1,0,1,0],[1e6,0,-1,0],[-2,6,-2,0]],
 [[9,-12,2/3,-1],[0,0,0,1],[1.5,0,-1,0],[-0.5,1,0,0],[-1,0,1,0],[11,-16,1,0],[16,-20,-1/3,0]],
]
def bourgain(s,t):
    best = -math.inf; covered=False
    for P in BP:
        lo, hi, ok = 0.0, math.inf, True
        for b,a1,a2,a3 in P:
            v = b+a1*s+a2*t
            if a3 == 0:
                if v < -1e-12: ok=False; break
            elif a3 < 0: hi = min(hi, v/(-a3))
            else: lo = max(lo, -v/a3)
        if ok and hi >= lo-1e-12:
            covered=True; best=max(best,hi)
    return best if covered else math.inf
LIT.append(('Bourgain', bourgain))
def L_base(s,t,new):
    v = min(f(s,t) for _,f in LIT)
    v = min(v, max(2-2*s, 2-2*s+t-mc_hm_cache[(s,new)]))
    return v
mc_hm_cache = {}
def Lstar_factory(s,new):
    mc_hm_cache[(s,new)] = mc_range_hm(s, g_new if new else g_old)
    # subdivision closure
    lvl = 2-2*s; t=0.0; step=0.0005
    while t < 6 and L_base(s,t+step,new) <= lvl+1e-12: t+=step
    tM = t
    def L1(t): return min(L_base(s,t,new), lvl+max(0.0,t-tM))
    def Ls(t): return min(k*L1(t/k) for k in range(1,9))
    return Ls, tM
def Z(s,t,Ls,new):  # zeta LV bound, t>=2
    if s > t*beta_low(1/t): return -math.inf
    if new and s > B: return -math.inf
    return min(Ls(t), 2*t-12*s+6, min(t-6*(s-0.5), 5*t-32*(s-0.5)))
def A_bound(s,new):
    Ls, tM = Lstar_factory(s,new)
    ts = [1.0+0.004*i for i in range(0,1500)]
    Lr = {t: Ls(t)/t for t in ts}
    Zr = {t: (Z(s,t,Ls,new)/t) for t in ts if t>=2}
    best = math.inf; bt=None
    for t0 in [1.0+0.004*i for i in range(0,700)]:
        sup1 = max(Lr[t] for t in ts if t0-1e-9 <= t <= 2*t0+1e-9)
        sup2 = max([Zr[t] for t in Zr if t < t0-1e-9] or [-math.inf])
        v = max(sup1,sup2)
        if v < best: best, bt = v, t0
    return best/(1-s), bt, tM, mc_hm_cache[(s,new)]
if __name__ == '__main__':
    sig = [float(x) for x in sys.argv[1:]] or [0.80,0.82,0.84,0.85,0.855,0.86,0.865,0.87,0.873]
    for s in sig:
        a0 = A_bound(s,False); a1 = A_bound(s,True)
        print(f"s={s:.4f} rec3/(2s)={1.5/s:.5f} | old: A<={a0[0]:.5f} t0={a0[1]:.3f} tMC={a0[2]:.4f} tHM={a0[3]:.4f} | new: A<={a1[0]:.5f} t0={a1[1]:.3f} tMC={a1[2]:.4f} tHM={a1[3]:.4f} tauB={tauB(s):.4f} | cor3 new={3/(max(tauB(s),3*s-1)+1-s):.5f}", flush=True)
