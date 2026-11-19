# compare.py (bothachain): every chain step against Botha's lower bounds mu and the float scan at his windows (../BothaAudit/logs).
# mu_odd = Uo/(Uo/mu_odd), mu_l2 = U2/(U2/mu_l2) from points.jsonl (our rigorous Galerkin upper bounds Uo, U2, point.py) and compare2.out's ratios (4 digits).
# Checks (a violation refutes a certificate): for every window c_B <= c_{i+1}: eps_i <= lambda(c_{i+1}) <= lambda(c_B) <= Uo(c_B), U2(c_B);
# for every window c_B >= c_i: mu(c_B) <= lambda_1^odd(c_B) and lambda_1^even(c_B) <= U(c_i) (no relation forced; printed: mu(c_B)/U(c_i) > 1 is Botha's
# own step condition with our U, i.e. his record and ours agree on the ordering). Float scan: Uo/U and U2/U at the nearest window (no crossing: >> 1).
import json, math, os
from fractions import Fraction as F
H = os.path.dirname(os.path.abspath(__file__)); A = os.path.join(H, '..', 'BothaAudit', 'logs')
def lg(s): s = str(s).strip('[').split(' ')[0]; m, _, e = s.partition('e'); return math.log10(float(m)) + int(e)
P = {d['c']: d for d in map(json.loads, open(os.path.join(A, 'points.jsonl')))}
W = []
for l in open(os.path.join(A, 'compare2.out')):
    x = [t.strip() for t in l.split('|')]
    if len(x) > 7 and x[1] in P:
        d = P[x[1]]; W.append(dict(c=F(x[1]), cs=x[1], L=x[0], Uo=lg(d['Uo_up']), U2=lg(d['U2_up']), U=lg(d['U_up']),
                                   mo=lg(d['Uo_up']) - math.log10(float(x[6])), m2=lg(d['U2_up']) - math.log10(float(x[7]))))
R = [d for d in map(json.loads, open(os.path.join(H, 'logs', 'chain.jsonl'))) if d.get('ok')]
print('| step c_i → c_{i+1} | eps_i (our lower bound at c_{i+1}) | U(c_i) | Botha window c_B ≤ c_{i+1} (nearest) | Uo(c_B)/eps, U2(c_B)/eps (≥ 1 required) | next window c_B′ ≥ c_i | mu_odd(c_B′), mu_l2(c_B′) | min mu(c_B′)/U(c_i) | scan Uo/U, U2/U at c_B |')
print('|---|---|---|---|---|---|---|---|---|')
bad = 0
for i, d in enumerate(R):
    a, b = F(d['cprev']), F(d['c']); e = lg(d['eps']); u = lg(d['Uprev'])
    lo = [w for w in W if w['c'] <= b]; hi = [w for w in W if w['c'] >= a]
    for w in lo:   # every window below c_{i+1}: eps <= lambda(c_B) <= Uo(c_B), U2(c_B)
        if not (w['Uo'] >= e and w['U2'] >= e): bad += 1; print('VIOLATION', d['c'], w['cs'])
    s1 = s2 = s3 = '-'
    if lo:
        w = lo[-1]; s1 = '%s (L = %s)' % (w['cs'], w['L']); s2 = '%.3g, %.3g' % (10**(w['Uo'] - e), 10**(w['U2'] - e)); s3 = '%.3g, %.3g' % (10**(w['Uo'] - w['U']), 10**(w['U2'] - w['U']))
    t1 = t2 = t3 = '-'
    if hi:
        w = hi[0]; t1 = w['cs']; t2 = '%.4ge%d, %.4ge%d' % (10**(w['mo'] % 1), math.floor(w['mo']), 10**(w['m2'] % 1), math.floor(w['m2'])); t3 = '%.3g' % 10**(min(w['mo'], w['m2']) - u)
    print('| %s → %s | %s | %s | %s | %s | %s | %s | %s | %s |' % (d['cprev'], d['c'], d['eps'], d['Uprev'], s1, s2, t1, t2, t3, s3))
print()
print('violations: %d (windows checked: all c_B <= c_{i+1} for each of %d steps)' % (bad, len(R)))
