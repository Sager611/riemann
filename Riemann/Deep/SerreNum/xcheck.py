#!/usr/bin/env python3
"""Independent floating-point re-evaluation (mpmath, 40 digits; no code shared with cert.py or ChebExp) of the
interval inequality f(n) of DERIVATION §3.3 at the first interval of every chain and the tail inequality of §3.4.
Not rigorous (sups on a dense grid); its job is to catch implementation slips in cert.py.
  uv run --no-project --with mpmath python3 xcheck.py > logs/xcheck_output.txt
"""
import json, sys
sys.dont_write_bytecode = True
import mpmath as mp
mp.mp.dps = 40
B2 = (777 - mp.sqrt(1081)) / 852
DEL = {'1/32': mp.mpf(1) / 32, '1-(777-1081**0.5)/852': 1 - B2}

def Phi(s, h, ss):                      # ((1 - e^{-h(s-ss)})/(s-ss))^3
    z = s - ss
    return ((1 - mp.e ** (-h * z)) / z) ** 3 if z != 0 else mp.mpf(h) ** 3

def F(u, h, ss):                        # e^{ss u} B3(u), B3 = h^2 M3(u/h)
    v = u / h
    if v < 0 or v > 3: return mp.mpf(0)
    m = v * v / 2 if v <= 1 else ((-2 * v * v + 6 * v - 3) / 2 if v <= 2 else (3 - v) ** 2 / 2)
    return mp.e ** (ss * u) * h * h * m

def supF(h, ss, shift, N=3000):
    return max(F(3 * h * i / N, h, ss) * mp.e ** (-shift * 3 * h * i / N) for i in range(N + 1))

def supline(beta, h, ss, s0):           # sup_{gamma >= 0} |Phi(beta+i gamma)| ((s0-beta)^2+gamma^2)/(s0-beta)
    d = s0 - beta; best = mp.mpf(0)
    for k in range(4000):
        g = mp.mpf(k) / 100 if k < 2000 else mp.mpf(20) * mp.mpf(1.003) ** (k - 2000)
        v = abs(Phi(mp.mpc(beta, g), h, ss)) * (d * d + g * g) / d
        if v > best: best = v
    return best

def gfun(s0):
    a = (mp.log(mp.pi) - mp.digamma(s0 / 2)) / 2
    b = mp.log(2 * mp.pi) - mp.digamma(s0)
    return min(a, b / 2)

def pieces(h, ss, s0, delta, lx):
    U = 3 * h; x = mp.e ** lx
    P1 = Phi(mp.mpf(1), h, ss); P0 = Phi(mp.mpf(0), h, ss)
    A1 = supline(1 - delta, h, ss, s0); A2 = supline(delta, h, ss, s0)
    M = max(A1, x ** (2 * delta - 1) * A2)
    kn = supF(h, ss, 0); S2 = supF(h, ss, mp.mpf(1) / 2); S3 = supF(h, ss, mp.mpf(1) / 3)
    tau = sum(x ** (-j) * Phi(mp.mpf(-j), h, ss) for j in range(1, 60)) / P0
    Rf = min(kn * lx, kn * (U + mp.log(2)), P0 + kn * mp.log(2)) / mp.log(2)
    J = int(mp.floor(lx / mp.log(2)))
    Pk = mp.mpf('1.01624') * (x ** 0.5 * (S2 + Phi(mp.mpf(1) / 2, h, ss) / 2) + x ** (mp.mpf(1) / 3) * (S3 + Phi(mp.mpf(1) / 3, h, ss) / 3)
                             + max(J - 3, 0) * x ** 0.25 * (kn + P0 / 4))
    cB = 1 / s0 + 1 / (s0 - 1)
    return dict(U=U, x=x, P1=P1, P0=P0, M=M, kn=kn, S2=S2, S3=S3, tau=tau, Rf=Rf, Pk=Pk, cB=cB, g=gfun(s0))

def main():
    chains = json.load(open('chains.json'))
    worst = None
    for ch in chains:
        delta = DEL[ch['delta_expr']]; a = mp.mpf(ch['a_str']); b = mp.mpf(ch['b_str']); n_max = ch['nmax']
        nl = n_max * mp.log(n_max)
        iv = ch['intervals'][0]
        h, ss, s0 = mp.mpf(iv['h']), mp.mpf(iv['ss']), mp.mpf(iv['s0'])
        lx = mp.log(a * iv['rho_a'] + b) / delta
        P = pieces(h, ss, s0, delta, lx)
        Lb = (n_max - 1) * mp.mpf(iv['rho_b']) + nl
        vals = []
        for n in (2, n_max):
            B = Lb / 2 + P['cB'] - n * P['g']
            v = P['P1'] - P['x'] ** (-delta) * P['M'] * B - ((n - 1) * P['P0'] + n * P['P0'] * P['tau'] + 2 * Lb * P['Rf']
                                                             + n * P['Pk'] + n * P['kn'] * iv['rho_b']) / P['x']
            vals.append(v / P['P1'])
        r1 = min(vals)
        tv = ch['intervals'][-1]
        h, ss, s0 = mp.mpf(tv['h']), mp.mpf(tv['ss']), mp.mpf(tv['s0'])
        rN = mp.mpf(tv['rho_a']); XN = a * rN + b; lx = mp.log(XN) / delta
        P = pieces(h, ss, s0, delta, lx)
        al = mp.mpf(n_max - 1) / 2; be = nl / 2 + P['cB'] - 2 * P['g']
        zr = max(al / a, (al * rN + be) / XN)
        Rfc = min(P['kn'] * (P['U'] + mp.log(2)), P['P0'] + P['kn'] * mp.log(2)) / mp.log(2)
        x = P['x']
        powx = mp.mpf('1.01624') * (x ** -0.5 * (P['S2'] + Phi(mp.mpf(1) / 2, h, ss) / 2) + x ** (-mp.mpf(2) / 3) * (P['S3'] + Phi(mp.mpf(1) / 3, h, ss) / 3)
                                    + 2 * lx * x ** (-mp.mpf(3) / 4) * (P['kn'] + P['P0'] / 4))
        t = P['P1'] - P['M'] * zr - (n_max * P['P0'] * (1 + P['tau']) + 2 * ((n_max - 1) * rN + nl) * Rfc + n_max * P['kn'] * rN) / x - n_max * powx
        r2 = t / P['P1']
        print(f"{ch['name']:5s} first interval [{iv['rho_a']:.4f}, {iv['rho_b']:.4f}]: LHS/P1 = {mp.nstr(r1, 6)};  tail from {tv['rho_a']:.3f}: LHS/P1 = {mp.nstr(r2, 6)}", flush=True)
        for r in (r1, r2):
            worst = r if worst is None or r < worst else worst
    print('all positive:', worst > 0, ' smallest:', mp.nstr(worst, 6))

if __name__ == '__main__':
    main()
