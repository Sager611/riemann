#!/usr/bin/env python3
"""Independent branch-and-bound for the smaller strata (atom coordinates, no KKT):
prove sup |a_n| < 2/e over a box of x = (lam_1..lam_N, theta_2..theta_N), theta_1 = 0, inside
D_N = {lam_1 >= lam_j, theta ordered, theta_2+theta_N <= 2pi, T_lo <= t <= T_hi}.
Own code (python-flint Arb): tests S, T, V (centred mean-value bound of a_n = e^{-t} u_n with
u = exp(sum h_k z^k), h_k = -2 sum lam_j e^{ik theta_j}), P (Cauchy/Poisson upper Riemann sum).
usage: kz_atoms.py root n N [limit]      (whole stratum)
       kz_atoms.py jobs n N logfile K maxnodes seed   (K random completed job boxes from a job log)"""
import sys, math, json, random, time
from flint import arb, acb, ctx
import kz_check as KC

ctx.prec = 96
PI = arb.pi()
TWO_PI = 2 * PI
TWO_E_LO = (2 / arb(1).exp()).lower()
TLO = KC.TLO
THI = KC.THI

def unit_ranges(lo, hi, n):
    out = []
    for k in range(1, n + 1):
        c, s = KC.trig_range(lo, hi, k)
        out.append(acb(c, s))
    return out

def u_coeffs(lams, E, n):
    """lams: list of arb; E[j][k-1] = e^{ik theta_j} (acb). returns u_0..u_n (acb)."""
    N = len(lams)
    h = [acb(0)] + [sum((-2 * lams[j] * E[j][k - 1] for j in range(N)), acb(0)) for k in range(1, n + 1)]
    u = [acb(1)]
    for k in range(1, n + 1):
        u.append(sum((m * h[m] * u[k - m] for m in range(1, k + 1)), acb(0)) / k)
    return u

def an_bound(B, n, N):
    lam = [KC.iv(B[j][0], B[j][1]) for j in range(N)]
    E = [[acb(1)] * n] + [unit_ranges(B[N + j - 1][0], B[N + j - 1][1], n) for j in range(1, N)]
    t = sum(lam, arb(0)); et = (-t).exp()
    u = u_coeffs(lam, E, n)
    nat = (et * u[n]).abs_upper()
    # centre
    c = [(B[i][0] + B[i][1]) / 2 for i in range(len(B))]
    r = [(B[i][1] - B[i][0]) / 2 for i in range(len(B))]
    lamc = c[:N]
    Ec = [[acb(1)] * n] + [[acb((k * c[N + j - 1]).cos(), (k * c[N + j - 1]).sin()) for k in range(1, n + 1)] for j in range(1, N)]
    ac = ((-sum(lamc, arb(0))).exp() * u_coeffs(lamc, Ec, n)[n]).abs_upper()
    contrib = []
    for j in range(N):        # d/d lam_j
        d = sum((-2 * E[j][k - 1] * u[n - k] for k in range(1, n + 1)), acb(0)) - u[n]
        contrib.append((et * d).abs_upper() * r[j])
    for j in range(1, N):     # d/d theta_j
        d = sum((-2 * acb(0, k) * lam[j] * E[j][k - 1] * u[n - k] for k in range(1, n + 1)), acb(0))
        contrib.append((et * d).abs_upper() * r[N + j - 1])
    mv = ac + sum(contrib, arb(0))
    return nat.min(mv), contrib

def circ_maxdist(a, b, c, d):
    """max over phi in [a,b], th in [c,d] of circular distance (as arb upper bound, <= pi)."""
    L, U = a - d, b - c
    k0 = int(math.floor(float(((L - PI) / TWO_PI).lower()))) - 1
    k1 = int(math.ceil(float(((U - PI) / TWO_PI).upper()))) + 1
    for k in range(k0, k1 + 1):
        x = PI + 2 * k * PI
        if not (x.upper() < L.lower() or x.lower() > U.upper()):
            return PI
    def dist(v):
        m = v - TWO_PI * round(float(v / TWO_PI))
        return abs(m)
    return dist(L).max(dist(U)).min(PI)

def poisson_bound(B, n, N, rs=(0.80, 0.88, 0.93, 0.96), M=96):
    lamlo = [B[j][0].max(arb(0)) for j in range(N)]
    ths = [(arb(0), arb(0))] + [(B[N + j - 1][0], B[N + j - 1][1]) for j in range(1, N)]
    best = None
    for rr in rs:
        r = arb(rr); s = arb(0)
        for i in range(M):
            a = TWO_PI * i / M; b = TWO_PI * (i + 1) / M
            e = arb(0)
            for j in range(N):
                dm = circ_maxdist(a, b, ths[j][0], ths[j][1])
                P = (1 - r * r) / (1 - 2 * r * dm.cos() + r * r)
                e += lamlo[j] * P.lower()
            s += (-e).exp()
        v = (s / M / r ** n).upper()
        best = v if best is None else best.min(v)
    return best

def test(B, n, N):
    T = arb(THI[n])
    for j in range(1, N):
        if B[0][1] < B[j][0]: return 'S'
    for j in range(1, N - 1):
        if B[N + j - 1][0] > B[N + j][1]: return 'S'
    if N >= 2 and (B[N][0] + B[2 * N - 2][0]).lower() > TWO_PI.upper(): return 'S'
    tlo = sum((B[j][0] for j in range(N)), arb(0)); thi = sum((B[j][1] for j in range(N)), arb(0))
    if thi < TLO or tlo > T: return 'T'
    bnd, contrib = an_bound(B, n, N)
    if bnd < TWO_E_LO: return 'V'
    if tlo >= 1.5:
        if poisson_bound(B, n, N) < TWO_E_LO: return 'P'
    return contrib

def split(B, contrib, n, N):
    # coordinate of largest mean-value contribution; fall back to widest scaled
    scores = [float(c) for c in contrib]
    i = max(range(len(B)), key=lambda k: scores[k])
    if scores[i] <= 0:
        i = max(range(len(B)), key=lambda k: float(B[k][1] - B[k][0]))
    mid = (B[i][0] + B[i][1]) / 2
    a = [list(x) for x in B]; b = [list(x) for x in B]
    a[i][1] = mid; b[i][0] = mid
    return a, b

def prove(B, n, N, limit):
    stack = [B]; nodes = 0; cnt = {}
    while stack:
        X = stack.pop(); nodes += 1
        res = test(X, n, N)
        if isinstance(res, str):
            cnt[res] = cnt.get(res, 0) + 1; continue
        if nodes >= limit: return False, nodes, cnt
        stack.extend(split(X, res, n, N))
    return True, nodes, cnt

def root(n, N):
    T = THI[n]
    B = [[arb(0), arb(T)]] + [[arb(0), arb(T / 2)] for _ in range(N - 1)]
    if N >= 2:
        B += [[arb(0), PI.upper()]] + [[arb(0), TWO_PI.upper()] for _ in range(N - 2)]
    return B

def jobs(logfile, N):
    J = {}; D = []
    for line in open(logfile):
        f = line.split()
        if f[0] == 'J':
            vals = [float.fromhex(x) for x in f[4:]]; d = len(vals) // 2
            J[int(f[1])] = (vals[:d], vals[d:])
        elif f[0] == 'D':
            D.append((int(f[1]), int(f[2])))
    return J, D

if __name__ == '__main__':
    mode = sys.argv[1]; n = int(sys.argv[2]); N = int(sys.argv[3])
    if mode == 'root':
        limit = int(sys.argv[4]) if len(sys.argv) > 4 else 200000
        t0 = time.time(); ok, nodes, cnt = prove(root(n, N), n, N, limit)
        print(json.dumps(dict(mode='root', n=n, N=N, proved=ok, nodes=nodes, closed=cnt, sec=round(time.time() - t0, 1))))
    else:
        logfile, K, maxnodes, seed = sys.argv[4], int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7])
        J, D = jobs(logfile, N); rng = random.Random(seed)
        pool = [x for x in D if x[1] <= maxnodes]
        pick = rng.sample(pool, min(K, len(pool)))
        res = []
        for jid, pn in pick:
            lo, hi = J[jid]; B = [[arb(lo[i]), arb(hi[i])] for i in range(len(lo))]
            t0 = time.time(); ok, nodes, cnt = prove(B, n, N, 40 * pn + 2000)
            res.append(dict(job=jid, prod_nodes=pn, proved=ok, my_nodes=nodes, closed=cnt, sec=round(time.time() - t0, 1)))
            print(json.dumps(res[-1]), flush=True)
        print(json.dumps(dict(mode='jobs', n=n, N=N, jobs_in_log=len(D), eligible=len(pool), checked=len(res), proved=sum(r['proved'] for r in res))))
