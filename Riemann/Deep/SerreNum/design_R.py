#!/usr/bin/env python3
"""Design (floating point, NOT part of the proof): chains of intervals of rho = log R (R = rad(2 N_E N_E'))
with kernels, for the Chebotarev step with excluded primes (DERIVATION section 3).  Output JSON for cert.py.

For an interval [ra, rb]:  x = (a*ra + b)^(1/delta),  Lambda <= Lb := (nmax-1)*rb + nmax*log(nmax),  log m <= rb,
LHS(n) = P1 - x^-delta M (Lb/2 + cB - n g) - x^-1 [(n-1)P0 + n P0 eps + 2 Lb Rf + n pw + n kn rb] > margin*P1
for n in {2, nmax}.  Reuses ../ChebExp/design.py (kernel constants, Nelder-Mead).
  usage: python3 design_R.py NAME DELTA NMAX A B RHO0 h,ss,s0 [margin]
         python3 design_R.py --f DELTA NMAX rho1,rho2,...   (profile: least x^delta per rho)
"""
import sys, os, math, json
sys.dont_write_bytecode = True     # leave no __pycache__ in ../ChebExp
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'ChebExp'))
import design as D   # noqa: E402  (ChebExp floating-point design helpers)

LOG2 = math.log(2.0)

def _triv_ratio_log(K, lx, h, ss, J1=30):
    """float version of tau_x computed in log space (ChebExp design.triv_ratio overflows for wide kernels)."""
    tot = 0.0
    for j in range(1, J1 + 1):
        z = j + ss                       # Phi(-j) = ((e^{h z} - 1)/z)^3
        lphi = 3 * (h * z + math.log1p(-math.exp(-h * z)) - math.log(z))
        tot += math.exp(min(700.0, -j * lx + lphi))
    lq = 3 * h - lx
    tail = math.exp(min(700.0, (J1 + 1) * lq)) / (1 - math.exp(lq)) if lq < 0 else float('inf')
    return tot / K['P0'] + tail

D.triv_ratio = _triv_ratio_log

def lam_max(nmax, rho):
    return (nmax - 1) * rho + nmax * math.log(nmax)

def rho_b(delta, nmax, lx, h, ss, s0, margin):
    """largest rb such that LHS >= margin*P1 at x = e^lx, for n in {2, nmax} (affine in rb)."""
    if h <= 0 or s0 <= 1.0001 or s0 > 6: return -1e9
    K = D.consts(delta, h, ss, s0)
    U = K['U']
    if lx <= U + 0.01: return -1e9
    xd = math.exp(-delta * lx); M = max(K['A1'], math.exp((2 * delta - 1) * lx) * K['A2'])
    ex = math.exp(-lx); eps = D.triv_ratio(K, lx, h, ss)
    Rf = min(K['kn'] * lx, K['kn'] * (U + LOG2), K['P0'] + K['kn'] * LOG2) / LOG2
    pw = D.POW(lx, h, ss, K['kn'])
    nl = nmax * math.log(nmax)
    res = []
    for n in (2.0, float(nmax)):
        a0 = K['P1'] * (1 - margin) - xd * M * (nl / 2 + K['cB'] - n * K['g']) \
            - ex * ((n - 1) * K['P0'] + n * K['P0'] * eps + 2 * nl * Rf + n * pw)
        a1 = -xd * M * (nmax - 1) / 2 - ex * (2 * (nmax - 1) * Rf + n * K['kn'])
        res.append(-a0 / a1)
    return min(res)

def best_kernel(delta, nmax, lx, p, margin):
    f = lambda q: -rho_b(delta, nmax, lx, q[0], q[1], q[2], margin)
    starts = [np.array(p, float)]
    if 3 * p[0] > 0.9 * lx:                      # previous kernel too wide for this x: add scaled starts
        starts += [np.array([0.22 * lx, 0.85, 1.06]), np.array([0.25 * lx, 0.8, 1.1])]
    best = None
    for s in starts:
        q, v = D.nm(f, s, [0.2 * s[0], 0.1, 0.2 * (s[2] - 1)], 120)
        if best is None or v < best[1]: best = (q, v)
    return best[0], -best[1]

def tail_val(delta, nmax, a, b, rN, h, ss, s0, margin):
    """relative margin of the tail check (valid for all rho >= rN), mirrors cert.py lhs_tail_R."""
    if h <= 0 or s0 <= 1.0001 or s0 > 6: return -1e9
    K = D.consts(delta, h, ss, s0)
    XN = a * rN + b; lx = math.log(XN) / delta; U = K['U']
    if lx <= U + 0.5 or math.exp(lx) <= D.XEXACT: return -1e9
    M = max(K['A1'], math.exp((2 * delta - 1) * lx) * K['A2'])
    eps = D.triv_ratio(K, lx, h, ss)
    nl = nmax * math.log(nmax)
    al = (nmax - 1) / 2.0; be = nl / 2 + K['cB'] - 2 * K['g']
    zr = max(al / a, (al * rN + be) / XN)
    Rfc = min(K['kn'] * (U + LOG2), K['P0'] + K['kn'] * LOG2) / LOG2
    ex = math.exp(-lx)
    u = np.linspace(0, 3 * h, 6001)
    S2 = float(np.max(D.F(u, h, ss) * np.exp(-u / 2))); S3 = float(np.max(D.F(u, h, ss) * np.exp(-u / 3)))
    P2 = D.Phi(np.array([0.5 + 0j]), h, ss)[0].real; P3 = D.Phi(np.array([1 / 3 + 0j]), h, ss)[0].real
    powx = 1.01624 * (math.exp(-lx / 2) * (S2 + P2 / 2) + math.exp(-2 * lx / 3) * (S3 + P3 / 3)
                      + 2 * lx * math.exp(-3 * lx / 4) * (K['kn'] + K['P0'] / 4))
    LbN = lam_max(nmax, rN)
    val = K['P1'] * (1 - margin) - M * zr - ex * (nmax * K['P0'] * (1 + eps) + 2 * LbN * Rfc + nmax * K['kn'] * rN) - nmax * powx
    return val / K['P1']

def chain(delta, nmax, a, b, rho0, p0, margin=3e-3, maxsteps=400):
    out = []; ra = rho0; p = np.array(p0, float)
    for step in range(maxsteps):
        ft = lambda q: -tail_val(delta, nmax, a, b, ra, q[0], q[1], q[2], margin)
        qt, vt = D.nm(ft, p, [0.2 * p[0], 0.1, 0.2 * (p[2] - 1)], 100)
        if -vt > 0:
            out.append(dict(kind='tail', rho_a=ra, h=float(qt[0]), ss=float(qt[1]), s0=float(qt[2])))
            print(f'  tail from rho={ra:.5f} h={qt[0]:.4f} ss={qt[1]:.4f} s0={qt[2]:.4f} slack={-vt:.4f}', flush=True)
            return out
        lx = math.log(a * ra + b) / delta
        q, rb = best_kernel(delta, nmax, lx, p, margin)
        if rb <= ra + 1e-6:
            print(f'  STUCK at rho={ra:.5f} best rb={rb}', flush=True); return None
        rb = min(rb, 1e6)
        out.append(dict(kind='int', rho_a=ra, rho_b=rb * (1 - 1e-9), h=float(q[0]), ss=float(q[1]), s0=float(q[2])))
        print(f'  [{ra:.5f}, {rb:.5f}] h={q[0]:.4f} ss={q[1]:.4f} s0={q[2]:.4f} logx={lx:.3f}', flush=True)
        ra = rb * (1 - 1e-9); p = q
    return None

def profile(delta, nmax, rhos, p0=(30.0, 0.85, 1.06)):
    """for each rho: least X = x^delta (bisection) such that some kernel gives rho_b >= rho."""
    p = np.array(p0, float)
    for rho in rhos:
        lo_, hi_ = 1.0, 1e5
        # find feasible hi
        for it in range(40):
            mid = math.sqrt(lo_ * hi_)
            q, rb = best_kernel(delta, nmax, math.log(mid) / delta, p, 0.0)
            if rb >= rho: hi_ = mid; p = q
            else: lo_ = mid
            if hi_ / lo_ < 1.002: break
        L = lam_max(nmax, rho)
        print(f'rho={rho:8.3f} Lambda_max={L:10.2f} X=x^delta={hi_:10.3f} X/(Lmax+log3)={hi_/(L+math.log(3)):.5f} X/rho={hi_/rho:.3f} kernel={p}', flush=True)

if __name__ == '__main__':
    if sys.argv[1] == '--f':
        delta = eval(sys.argv[2]); nmax = int(sys.argv[3]); rhos = [float(v) for v in sys.argv[4].split(',')]
        p0 = tuple(float(v) for v in sys.argv[5].split(',')) if len(sys.argv) > 5 else (30.0, 0.85, 1.06)
        profile(delta, nmax, rhos, p0)
        sys.exit(0)
    name = sys.argv[1]; dexpr = sys.argv[2]; delta = eval(dexpr); nmax = int(sys.argv[3])
    a = float(sys.argv[4]); b = float(sys.argv[5]); rho0 = float(sys.argv[6])
    p0 = [float(v) for v in sys.argv[7].split(',')]
    margin = float(sys.argv[8]) if len(sys.argv) > 8 else 3e-3
    res = chain(delta, nmax, a, b, rho0, p0, margin)
    json.dump(dict(name=name, delta_expr=dexpr, nmax=nmax, a=a, b=b, rho0=rho0, intervals=res),
              open(f'chain_{name}.json', 'w'), indent=1)
    print('done', name, 'ok' if res else 'FAILED', len(res) if res else 0)
