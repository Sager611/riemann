# b3nu task (1): cancellation in the dual total of Cor B' (B3Cert FULLPROOF l. 56-63), zero shifts, Levinson mollifier P(x) = x.
# Dual total  D = sum_{g,p,q} b_{gq} b_{gp} (4pi/(gpq)) S(p,q),  S(p,q) = sum_n d(n) e(s pbar n/q) T(n/(pq)),  b_h = mu(h) P(log(y/h)/log y).
# s = +1 is the oscillatory (J/Y) branch twist of E1Audit/voronoi_check.py (e(-Abar n/q), A = -p); s = -1 reported as a control.
# Ratios: Rn = sum|w||S| / sum|w| sum_n d|T| (n-sum cancellation), Rpq = |D| / sum|w||S| ((g,p,q)-sum cancellation), Rtot = |D| / sum|w| sum_n d|T|.
import sys, math, time, numpy as np
from b3nu_T import Tinterp
def mobius(N):
    mu = np.ones(N + 1, int); mu[0] = 0; pr = np.ones(N + 1, bool); pr[:2] = False
    for p in range(2, N + 1):
        if pr[p]:
            pr[2 * p::p] = False; mu[p::p] *= -1; mu[p * p::p * p] = 0
    return mu
def divisors(N):
    d = np.zeros(N + 1, np.int64)
    for k in range(1, N + 1): d[k::k] += 1
    return d
T, th = float(sys.argv[1]), float(sys.argv[2]); ys = [int(v) for v in sys.argv[3].split(',')]; rmax = 4.0
# optional argv[4] = params_exp010.json: frozen SL/EXP-010 mollifier polynomial P(x) = S(r_P x)/S(r_P), S = sum_{i<7} x^{2i+1}/(2i+1)! (as b3c_arb.py)
Pm = lambda x: x
if len(sys.argv) > 4:
    import json
    from fractions import Fraction
    ps = [q for q in json.load(open(sys.argv[4]))['parameter_sets'] if q['nu'] == '999/10000'][0]; rP = float(Fraction(ps['r_P']))
    S7 = lambda x: sum(x ** (2 * i + 1) / math.factorial(2 * i + 1) for i in range(7)); Pm = lambda x: S7(rP * x) / S7(rP)
    print(f"mollifier: frozen EXP-010 P, r_P = {rP:.6f}, P(1) = {Pm(1.0):.6f}, P(0.5) = {Pm(0.5):.6f}", flush=True)
H = T ** th; tau = T / H ** 2; Y = max(ys); t0 = time.time()
mu = mobius(Y); dn = divisors(int(rmax * tau * Y * Y) + 2)
TI = Tinterp(T, H, 0.9 / (Y * Y), rmax * tau * 1.02)
pairs = [(p, q) for p in range(1, Y + 1) for q in range(1, Y + 1) if math.gcd(p, q) == 1]
S = {}; A = {}; Rand = {}
for (p, q) in pairs:
    nm = int(rmax * tau * p * q)
    if nm < 1: S[p, q] = {1: 0j, -1: 0j}; A[p, q] = 0.0; Rand[p, q] = 0.0; continue
    n = np.arange(1, nm + 1); Tv = TI(n / (p * q)); dd = dn[1:nm + 1]
    pb = pow(p, -1, q) if q > 1 else 0
    S[p, q] = {s: np.sum(dd * np.exp(2j * np.pi * s * ((pb * n) % q) / q) * Tv) for s in (1, -1)}
    A[p, q] = float(np.sum(dd * np.abs(Tv))); Rand[p, q] = float(np.sqrt(np.sum(dd ** 2 * np.abs(Tv) ** 2)))
print(f"T={T:.0e} theta={th} H={H:.1f} tau={tau:.3f} nu_crit={(3*th-1)/4:.4f} y_crit={T**((3*th-1)/4):.2f}  pairs={len(pairs)}  setup+S {time.time()-t0:.0f}s", flush=True)
for y in ys:
    L = math.log(y) if y > 1 else 1.0
    b = lambda h: mu[h] * Pm(math.log(y / h) / L if y > 1 else 1.0) if h <= y else 0.0
    tot = {1: 0j, -1: 0j}; absall = 0.0; abspq = {1: 0.0, -1: 0.0}; rnd = 0.0; npair = 0; perpair = []
    for g in range(1, y + 1):
        for (p, q) in pairs:
            if g * p > y or g * q > y: continue
            wt = b(g * q) * b(g * p) * 4 * math.pi / (g * p * q)
            if wt == 0 or A[p, q] == 0: continue
            npair += 1; absall += abs(wt) * A[p, q]; rnd += abs(wt) * Rand[p, q]
            for s in (1, -1):
                tot[s] += wt * S[p, q][s]; abspq[s] += abs(wt) * abs(S[p, q][s])
            perpair.append(abs(S[p, q][1]) / A[p, q])
    main = math.sqrt(math.pi) * H
    pp = np.array(perpair) if perpair else np.array([0.0])
    print(f" y={y:3d} nu={math.log(y)/math.log(T):.4f} terms(g,p,q)={npair:4d}: sum|.| / (sqrt(pi)H) = {absall/main:.3e};  "
          f"s=+1: Rn={abspq[1]/absall:.3f} Rpq={abs(tot[1])/abspq[1]:.3f} Rtot={abs(tot[1])/absall:.4f} |D|/(sqrt(pi)H)={abs(tot[1])/main:.3e};  "
          f"s=-1: Rtot={abs(tot[-1])/absall:.4f};  random-n model Rn={rnd/absall:.3f};  per-pair |S|/A median {np.median(pp):.3f} max {pp.max():.3f}", flush=True)
print(f"done {time.time()-t0:.0f}s", flush=True)
