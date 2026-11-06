#!/usr/bin/env python3
"""Independent floating-point cross-checks (not part of the proof).
  uv run --no-project --with numpy --with mpmath python3 xcheck.py > logs/xcheck_output.txt
1. maximum principle (DERIVATION 3.4): the max of R over the interior of the strip does not exceed the max on the
   two boundary lines (dense grid), for the first kernel of each uniform chain;
2. an independent re-evaluation (mpmath quadrature for Phi, dense grids for the sups, a direct prime-power sum)
   of LHS(x; Lambda_b, 2) / Phi(1) for those intervals, to compare with logs/cert_output.txt."""
import json, math, numpy as np, mpmath as mp
mp.mp.dps = 30
ch = {c['name']: c for c in json.load(open('chains.json'))}
B2 = (777 - math.sqrt(1081)) / 852

def M3(v):
    v = np.asarray(v, float)
    return np.where((v >= 0) & (v < 1), v * v / 2, np.where((v >= 1) & (v < 2), (-2 * v * v + 6 * v - 3) / 2,
                    np.where((v >= 2) & (v <= 3), (3 - v) ** 2 / 2, 0.0)))

for name, delta in [('CTu', 1 / 32), ('Lu', 1 - B2)]:
    iv = ch[name]['intervals'][0]; c = ch[name]['c']
    h, ss, s0 = iv['h'], iv['ss'], iv['s0']; U = 3 * h
    Phi = lambda s: ((1 - np.exp(-h * (s - ss))) / (s - ss)) ** 3
    lx = math.log(c * (iv['Lam_a'] + math.log(3))) / delta; x = math.exp(lx)
    # 1. maximum principle
    bet = np.linspace(delta, 1 - delta, 301); gam = np.linspace(1e-6, 30, 6001)
    Bm, Gm = np.meshgrid(bet, gam)
    S = Bm + 1j * Gm
    R = np.exp((Bm - 1 + delta) * lx) * np.abs(Phi(S)) * ((s0 - Bm) ** 2 + Gm ** 2) / (s0 - Bm)
    print(f'{name}: max R interior grid = {R[:, 1:-1].max():.6g}; max on beta=1-delta: {R[:, -1].max():.6g}; on beta=delta: {R[:, 0].max():.6g}'
          f'  -> interior <= boundary: {R[:, 1:-1].max() <= max(R[:, -1].max(), R[:, 0].max()) * (1 + 1e-9)}')
    # 2. independent LHS
    F = lambda u: mp.e ** (ss * u) * h * h * mp.mpf(float(M3(u / h)))
    P1 = mp.quad(lambda u: F(u) * mp.e ** (-u), [0, h, 2 * h, 3 * h])
    P0 = mp.quad(F, [0, h, 2 * h, 3 * h])
    print(f'   Phi(1): quadrature {mp.nstr(P1, 12)} vs closed form {Phi(1.0).real:.12g};  Phi(0): {mp.nstr(P0, 12)} vs {Phi(1e-300 + 0j).real:.12g}')
    uu = np.linspace(0, U, 300001); Fv = np.exp(ss * uu) * h * h * M3(uu / h); kn = Fv.max()
    gg = np.linspace(0, 60, 1200001)
    A1 = (np.abs(Phi((1 - delta) + 1j * gg)) * ((s0 - 1 + delta) ** 2 + gg ** 2) / (s0 - 1 + delta)).max()
    A2 = (np.abs(Phi(delta + 1j * gg)) * ((s0 - delta) ** 2 + gg ** 2) / (s0 - delta)).max()
    M = max(A1, math.exp((2 * delta - 1) * lx) * A2)
    g = min((math.log(math.pi) - float(mp.digamma(s0 / 2))) / 2, (math.log(2 * math.pi) - float(mp.digamma(s0))) / 2)
    Lb = iv['Lam_b']; n = 2
    B = Lb / 2 + 1 / s0 + 1 / (s0 - 1) - n * g
    # prime powers directly
    N = int(math.isqrt(int(x))) + 2
    sv = np.ones(N + 1, bool); sv[:2] = False
    for i in range(2, int(N ** 0.5) + 1):
        if sv[i]: sv[i * i::i] = False
    ps = np.nonzero(sv)[0].astype(float)
    Pk = 0.0; j = 2
    while 2.0 ** j <= x:
        u = lx - j * np.log(ps); ok = (u >= 0) & (u <= U)
        Pk += float(np.sum(np.log(ps[ok]) * np.exp(ss * u[ok]) * h * h * M3(u[ok] / h))); j += 1
    tau = sum(math.exp(-j * lx + 3 * (h * (j + ss) + math.log(-math.expm1(-h * (j + ss))) - math.log(j + ss))) for j in range(1, 80)) / float(P0)
    Rf = min(kn * lx, kn * (U + math.log(2)), float(P0) + kn * math.log(2)) / math.log(2)
    low = ((n - 1) * float(P0) + n * float(P0) * tau + 2 * Lb * Rf + n * Pk) / x
    lhs = float(P1) - math.exp(-delta * lx) * M * B - low
    print(f'   independent: A1={A1:.6g} kn={kn:.6g} P_k(x)={Pk:.6g} at log x={lx:.4f}: LHS(x; Lb, 2)/Phi(1) = {lhs / float(P1):.5f} (cert: see logs/cert_output.txt)')
