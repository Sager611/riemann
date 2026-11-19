#!/usr/bin/env python3
"""chebaudit: independent 256-bit recomputation (mpmath, own code; cert.py not imported) of LHS(x; Lb, n) of
ChebExp DERIVATION (5.1)/(5.3) for extreme and random kernels of chains.json. Sups A1, A2 by dense float grid +
256-bit golden-section refinement (heuristic, so values should sit just below cert's rigorous upper bounds).
Also: direct evaluation of the tail inequality at many Lambda >= Lambda_N (no monotonicity argument used)."""
import json, math, random, re, cmath, sys, time
import mpmath as mp
mp.mp.prec = 256
D = '/Users/adriansager/playground/riemann/Riemann/Deep/ChebExp/'
chains = json.load(open(D + 'chains.json'))
cert = {}; cur = None; idx = 0
for line in open(D + 'logs/cert_output.txt'):
    m = re.match(r'=== chain (\S+):', line)
    if m: cur = m.group(1); idx = 0; continue
    if cur and (line.startswith('  [') or line.startswith('  tail [')):
        g = lambda pat: (lambda mm: float(mm.group(1)) if mm else None)(re.search(pat, line))
        cert[(cur, idx)] = dict(margin=g(r'LHS/P1 >= ([-0-9.e+]+)'), A1=g(r'A1=([-0-9.e+]+)'), A1lo=g(r'\(lo ([-0-9.e+]+)'),
                                M=g(r' M=([-0-9.e+]+)'), P1=g(r'P1=([-0-9.e+]+)'), lx=g(r'log x>?=([-0-9.e+]+)'))
        idx += 1
LOG3, LOG2, LP2 = mp.log(3), mp.log(2), mp.log(mp.pi / 2)
DEL = {'1/32': mp.mpf(1) / 32, '1-(777-1081**0.5)/852': 1 - (777 - mp.sqrt(1081)) / 852}

def M3(v):
    if v <= 0 or v >= 3: return mp.mpf(0)
    if v <= 1: return v * v / 2
    if v <= 2: return (-2 * v * v + 6 * v - 3) / 2
    return (3 - v) ** 2 / 2
class Ker:
    def __init__(s, h, ss, s0):
        s.h, s.ss, s.s0 = mp.mpf(h), mp.mpf(ss), mp.mpf(s0); s.U = 3 * s.h
        s.hf, s.ssf, s.s0f = h, ss, s0
    def F(s, u): return mp.exp(s.ss * u) * s.h ** 2 * M3(u / s.h)
    def Phi(s, z0):
        z = z0 - s.ss
        if abs(z) < mp.mpf(10) ** -60: return s.h ** 3
        return ((1 - mp.exp(-s.h * z)) / z) ** 3
    def maxF(s, shift): return gmax(lambda u: mp.exp(-shift * u) * s.F(u), mp.mpf(0), s.U)
def gmax(f, a, b, it=170):
    gr = (mp.sqrt(5) - 1) / 2; c = b - gr * (b - a); d = a + gr * (b - a); fc, fd = f(c), f(d)
    for _ in range(it):
        if fc > fd: b, d, fd = d, c, fc; c = b - gr * (b - a); fc = f(c)
        else: a, c, fc = c, d, fd; d = a + gr * (b - a); fd = f(d)
    return max(fc, fd, f(a), f(b))
def supline(K, beta):
    bf = float(beta); hf, sf, s0f = K.hf, K.ssf, K.s0f; d = s0f - bf; e = math.exp(hf * (sf - bf))
    def Rfl(g):
        z = complex(bf - sf, g)
        ph = hf ** 3 if abs(z) < 1e-12 else abs(((1 - cmath.exp(-hf * z)) / z) ** 3)
        return ph * (d * d + g * g) / d
    R0 = Rfl(0.0); G = 1.0
    while (1 + e) ** 3 * (d * d + G * G) / (d * G ** 3) > R0 / 200: G *= 1.25
    step = min(0.002, (2 * math.pi / hf) / 400); n = int(G / step) + 2
    if n > 600000: step = G / 600000; n = 600002
    vals = [Rfl(i * step) for i in range(n)]
    cand = [i for i in range(n) if (i == 0 or vals[i] >= vals[i - 1]) and (i == n - 1 or vals[i] >= vals[i + 1])]
    cand.sort(key=lambda i: -vals[i])
    dm = mp.mpf(K.s0) - beta
    Rmp = lambda g: abs(K.Phi(mp.mpc(beta, g))) * (dm * dm + g * g) / dm
    best = max(gmax(Rmp, mp.mpf(max(0, (i - 1) * step)), mp.mpf((i + 1) * step), 110) for i in cand[:6])
    return best, n
PR = None
def primes_upto(N):
    b = bytearray([1]) * (N + 1); b[0] = b[1] = 0
    for i in range(2, int(N ** 0.5) + 1):
        if b[i]: b[i * i::i] = bytearray(len(range(i * i, N + 1, i)))
    return [i for i in range(N + 1) if b[i]]
def Pk_exact(K, lx):
    global PR
    xf = float(mp.exp(lx)); tot = mp.mpf(0); j = 2
    if PR is None: PR = primes_upto(1415000)
    while 2.0 ** j <= xf * 1.0001:
        lo_p = math.exp((float(lx) - float(K.U)) / j) * 0.999 - 2; hi_p = math.exp(float(lx) / j) * 1.001 + 2
        assert hi_p < PR[-1] or j > 2
        for p in PR:
            if p > hi_p: break
            if p < lo_p: continue
            u = lx - j * mp.log(p)
            if 0 < u < K.U: tot += mp.log(p) * K.F(u)
        j += 1
    return tot
def consts(K):
    C = dict(P1=mp.re(K.Phi(mp.mpf(1))), P0=mp.re(K.Phi(mp.mpf(0))), Ph2=mp.re(K.Phi(mp.mpf(1) / 2)), Ph3=mp.re(K.Phi(mp.mpf(1) / 3)),
             kn=K.maxF(0), S2=K.maxF(mp.mpf(1) / 2), S3=K.maxF(mp.mpf(1) / 3))
    C['g'] = min((mp.log(mp.pi) - mp.digamma(K.s0 / 2)) / 2, (mp.log(2 * mp.pi) - mp.digamma(K.s0)) / 2)
    C['cB'] = 1 / K.s0 + 1 / (K.s0 - 1)
    return C
def tau(K, C, lx):
    t = mp.mpf(0); j = 1
    while True:
        term = mp.exp(-j * lx) * mp.re(K.Phi(-j)) / C['P0']; t += term; j += 1
        if term < mp.mpf(10) ** -75 or j > 4000: return t
def Pk_smooth(K, C, lx):
    J = int(mp.floor(lx / LOG2))
    return mp.mpf('1.01624') * (mp.exp(lx / 2) * (C['S2'] + C['Ph2'] / 2) + mp.exp(lx / 3) * (C['S3'] + C['Ph3'] / 3)
                                 + max(J - 3, 0) * mp.exp(lx / 4) * (C['kn'] + C['P0'] / 4))
def LHS(delta, C, K, A1, A2, lx, Lam, n, Pk):
    x = mp.exp(lx); Mx = max(A1, mp.exp((2 * delta - 1) * lx) * A2)
    B = Lam / 2 + C['cB'] - n * C['g']
    Rf = min(C['kn'] * lx, C['kn'] * (K.U + LOG2), C['P0'] + C['kn'] * LOG2) / LOG2
    return C['P1'] - mp.exp(-delta * lx) * Mx * B - ((n - 1) * C['P0'] + n * C['P0'] * tau(K, C, lx) + 2 * Lam * Rf + n * Pk) / x

def eval_int(ch, k):
    iv = ch['intervals'][k]; delta = DEL[ch['delta_expr']]; c = mp.mpf(ch['c'])
    K = Ker(iv['h'], iv['ss'], iv['s0']); C = consts(K)
    A1, n1 = supline(K, 1 - delta); A2, n2 = supline(K, delta)
    La, Lb = mp.mpf(iv['Lam_a']), mp.mpf(iv['Lam_b'])
    lx = mp.log(c * (La + LOG3)) / delta; x = mp.exp(lx)
    assert lx > K.U
    sm = Pk_smooth(K, C, lx)
    if x <= mp.mpf('2e12'): Pk = Pk_exact(K, lx); mode = 'exact'; ratio = sm / Pk
    else: Pk = sm; mode = 'smooth'; ratio = None
    v = min(LHS(delta, C, K, A1, A2, lx, Lb, mp.mpf(2), Pk), LHS(delta, C, K, A1, A2, lx, Lb, Lb / LP2, Pk))
    Mx = max(A1, mp.exp((2 * delta - 1) * lx) * A2)
    return dict(margin=v / C['P1'], A1=A1, M=Mx, P1=C['P1'], lx=lx, mode=mode, sm_over_exact=ratio, g=C['g'])
def eval_tail(ch, k):
    iv = ch['intervals'][k]; delta = DEL[ch['delta_expr']]; c = mp.mpf(ch['c'])
    K = Ker(iv['h'], iv['ss'], iv['s0']); C = consts(K)
    A1, n1 = supline(K, 1 - delta); A2, n2 = supline(K, delta)
    LN = mp.mpf(iv['Lam_a']); ellN = LN + LOG3; lx = mp.log(c * ellN) / delta; x = mp.exp(lx)
    Mx = max(A1, mp.exp((2 * delta - 1) * lx) * A2)
    add = max(mp.mpf(0), C['cB'] - LOG3 / 2 - 2 * C['g']); nmax = LN / LP2
    powx = mp.mpf('1.01624') * (mp.exp(-lx / 2) * (C['S2'] + C['Ph2'] / 2) + mp.exp(-2 * lx / 3) * (C['S3'] + C['Ph3'] / 3)
                                + 2 * lx * mp.exp(-3 * lx / 4) * (C['kn'] + C['P0'] / 4))
    val = C['P1'] - Mx * (1 / (2 * c) + add / (c * ellN)) - (nmax * C['P0'] * (1 + tau(K, C, lx)) + 2 * C['kn'] * LN * lx / LOG2) / x - nmax * powx
    conds = ((1 - delta) * ellN > LOG3, 2 * delta * ellN / (ellN - LOG3) < 1, delta * (ellN / (ellN - LOG3) + 1 / mp.log(c * ellN)) < mp.mpf(3) / 4)
    direct = []
    for f in [1, 1.01, 1.1, 1.3, 1.6, 2, 3, 5, 10, 30, 100, 1e3, 1e4, 1e6, 1e9, 1e15]:
        Lam = LN * f; lxx = mp.log(c * (Lam + LOG3)) / delta
        Pk = Pk_smooth(K, C, lxx)
        w = min(LHS(delta, C, K, A1, A2, lxx, Lam, mp.mpf(2), Pk), LHS(delta, C, K, A1, A2, lxx, Lam, Lam / LP2, Pk)) / C['P1']
        direct.append((f, float(w)))
    return dict(margin=val / C['P1'], A1=A1, M=Mx, P1=C['P1'], lx=lx, conds=conds, direct=direct, g=C['g'])
def main():
    t0 = time.time(); random.seed(20261119)
    allidx = [(ci, k) for ci, ch in enumerate(chains) for k in range(len(ch['intervals']))]
    sel = set()
    for ci, ch in enumerate(chains):
        ivs = ch['intervals']; ints = [k for k, iv in enumerate(ivs) if iv['kind'] == 'int']
        sel |= {(ci, 0), (ci, len(ivs) - 1), (ci, ints[-1]), (ci, min(ints, key=lambda k: cert[(ch['name'], k)]['margin']))}
    sel |= set(random.sample(allidx, int(sys.argv[1]) if len(sys.argv) > 1 else 15))
    if len(sys.argv) > 2: sel = set(list(sorted(sel))[:int(sys.argv[2])])
    worst = 1e9; bad = 0
    print('chain idx kind  log_x  cert_margin  my_margin  my-cert  myA1/certA1  certA1lo<=myA1  myM/certM  pow  smooth/exact')
    for ci, k in sorted(sel):
        ch = chains[ci]; iv = ch['intervals'][k]; ce = cert[(ch['name'], k)]
        r = eval_int(ch, k) if iv['kind'] == 'int' else eval_tail(ch, k)
        mg = float(r['margin']); d = mg - ce['margin']
        flag = 'OK' if (mg > 0 and d > -6e-6) else 'CHECK'
        if flag != 'OK': bad += 1
        worst = min(worst, mg)
        print('%-8s %3d %-4s %8.3f  %.5f  %.6f  %+.2e  %.6f  %s  %.6f  %s  %s  %s' % (ch['name'], k, iv['kind'], float(r['lx']), ce['margin'], mg, d,
              float(r['A1']) / ce['A1'], (ce['A1lo'] is None) or (ce['A1lo'] <= float(r['A1']) * (1 + 1e-4)), float(r['M']) / ce['M'],
              r.get('mode', 'tail'), ('%.3f' % float(r['sm_over_exact'])) if r.get('sm_over_exact') else '-', flag), flush=True)
        if iv['kind'] == 'tail':
            print('    tail conds', r['conds'], ' direct min margin over Lambda = LN*f:', min(w for f, w in r['direct']), 'at f =', min(r['direct'], key=lambda t: t[1])[0], flush=True)
    print('entries', len(sel), 'flagged', bad, 'min my margin %.5f' % worst, '(%.0f s)' % (time.time() - t0))
main()
