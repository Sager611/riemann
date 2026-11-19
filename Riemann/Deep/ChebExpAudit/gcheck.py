#!/usr/bin/env python3
"""chebaudit: for all 405 kernels: g(s0) = min(a, b/2) (mpmath digamma) range, and Phi(1) by plain closed form AND by
mpmath quadrature of int_0^U F(u) e^{-u} du, compared with cert's printed P1."""
import json, re
import mpmath as mp
mp.mp.dps = 30
D = '/Users/adriansager/playground/riemann/Riemann/Deep/ChebExp/'
chains = json.load(open(D + 'chains.json'))
P1c = {}; cur = None; idx = 0
for line in open(D + 'logs/cert_output.txt'):
    m = re.match(r'=== chain (\S+):', line)
    if m: cur = m.group(1); idx = 0; continue
    if cur and (line.startswith('  [') or line.startswith('  tail [')):
        P1c[(cur, idx)] = float(re.search(r'P1=([-0-9.e+]+)', line).group(1)); idx += 1
def M3(v):
    if v <= 0 or v >= 3: return mp.mpf(0)
    if v <= 1: return v * v / 2
    if v <= 2: return (-2 * v * v + 6 * v - 3) / 2
    return (3 - v) ** 2 / 2
gs = []; worst = 0; worstq = 0; n = 0
for ch in chains:
    for k, iv in enumerate(ch['intervals']):
        h, ss, s0 = mp.mpf(iv['h']), mp.mpf(iv['ss']), mp.mpf(iv['s0'])
        g = min((mp.log(mp.pi) - mp.digamma(s0 / 2)) / 2, (mp.log(2 * mp.pi) - mp.digamma(s0)) / 2); gs.append((float(g), float(s0)))
        cf = ((1 - mp.exp(-h * (1 - ss))) / (1 - ss)) ** 3
        q = mp.quad(lambda u: mp.exp((ss - 1) * u) * h * h * M3(u / h), [0, h, 2 * h, 3 * h])
        worst = max(worst, abs(float(cf) / P1c[(ch['name'], k)] - 1)); worstq = max(worstq, abs(float(q / cf) - 1)); n += 1
print('kernels', n, '; g(s0) range [%.4f, %.4f] (at s0 = %.4f and %.4f)' % (min(gs)[0], max(gs)[0], min(gs)[1], max(gs)[1]))
print('max |closed-form Phi(1) / cert printed P1 - 1| = %.2e (cert prints 5 digits)' % worst)
print('max |quadrature Phi(1) / closed form - 1| = %.2e' % worstq)
