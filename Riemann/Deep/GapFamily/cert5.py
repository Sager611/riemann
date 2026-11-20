"""explore53: Arb certificate of per-gap validity (Lemma 3, exact5.py) for a J = 5 form; writes logs/cert5_<tag>.json.
Form = SDP form (logs/scanJ5_r2_d64_k3.355.json) rationalised (1e-12) and repaired: a0 += D0 (monotone: raises every
min_b), phi12 += E12, phi32 += E32, phi14 -= E32 (telescoping coefficients: no global effect).
Checks (128 bits): (1) p > 0 on [1, inf); (2) C = [[a3+phi32, phi34/2],[phi34/2, a5]] > 0 (exact rationals);
(3) leading coefficients of N11, N22, DET > 0 and X* = largest real root of N11, N22, DET (and 1): R(x) > 0 for x > X*,
so every gap g < pi/sqrt(X*) is valid termwise (S_odd >= 0, S_even + C/2 >= C/2 > 0);
(4) g in [pi/sqrt(X*), pi]: adaptive Arb cells, k <= K with K > sqrt(X*) (dropped terms are PSD), 2x2 PSD tests.
usage: cert5.py tag D0 E12 E32 [K]"""
import sys, json
from fractions import Fraction as Fr
from flint import arb, fmpq, fmpq_poly, ctx
sys.path.insert(0, '.')
from exact5 import numer5
from exact6 import numer6
import os
JJ = int(os.environ.get('JJ', '5'))
ctx.prec = 128
tag = sys.argv[1]; D0 = Fr(sys.argv[2]); E12 = Fr(sys.argv[3]); E32 = Fr(sys.argv[4])
K = int(sys.argv[5]) if len(sys.argv) > 5 else 61
import os; l = json.load(open(os.environ.get('SRC', 'logs/scanJ5_r2_d64_k3.355.json')))[0]
a = [Fr(x).limit_denominator(10 ** 12) for x in l['a%d' % JJ]]
ph = {k.strip('()').replace(', ', ''): Fr(v).limit_denominator(10 ** 12) for k, v in l['phi%d' % JJ].items()}
a[0] += D0; ph['12'] += E12; ph['32'] += E32; ph['14'] -= E32
p, N11, N12, N22, DET, C = (numer6 if JJ == 6 else numer5)(a, ph)
q = lambda x: fmpq(x.numerator, x.denominator)
P = fmpq_poly([q(x) for x in p])
ok = True
p1 = sum(p)
rp = P.complex_roots()
bad = [r for r, m in rp if r.imag.contains(0) and not (r.real < 1)]
print('(1) p(1) = %.6g, real roots >= 1: %d' % (float(p1), len(bad)))
ok &= p1 > 0 and not bad
detC = C[0][0] * C[1][1] - C[0][1] ** 2
print('(2) C11 = %.4g, det C = %.4g' % (float(C[0][0]), float(detC)))
ok &= C[0][0] > 0 and detC > 0
Xs = [arb(1)]
for name, Pn in (('N11', N11), ('N22', N22), ('DET', DET)):
    lead = Pn[-1]
    rr = [r for r, m in fmpq_poly([q(x) for x in Pn]).complex_roots() if r.imag.contains(0)]
    Xs += [r.real for r in rr]
    print('(3) %s: degree %d, leading %.4g, real roots %s' % (name, len(Pn) - 1, float(lead), [r.real.str(8) for r in rr]))
    ok &= lead > 0
Xstar = max(x.upper() for x in Xs)
pi = arb.pi()
gstar = pi / arb(Xstar).sqrt()
K = max(K, int(float(arb(Xstar).sqrt().upper())) + 2)
aa = lambda P: [arb(q(x)) for x in P]
pA, A11, A12, A22 = aa(p), aa(N11), aa(N12), aa(N22)
Cm = [[arb(q(C[0][0])), arb(q(C[0][1]))], [arb(q(C[1][0])), arb(q(C[1][1]))]]
def hv(P, x):
    out = arb(0)
    for c in reversed(P):
        out = out * x + c
    return out
def sums(G):
    So = [arb(0)] * 3; Se = [arb(0)] * 3
    for k in range(1, K + 1):
        x = (k * pi / G) ** 2
        d = 4 * hv(pA, x)
        t = (hv(A11, x) / d, hv(A12, x) / d, hv(A22, x) / d)
        if k % 2:
            So = [So[i] + t[i] for i in range(3)]
        else:
            Se = [Se[i] + t[i] for i in range(3)]
    Se = [Se[0] + Cm[0][0] / 2, Se[1] + Cm[0][1] / 2, Se[2] + Cm[1][1] / 2]
    return So, Se
def psd(S):
    return S[0] > 0 and S[2] > 0 and S[0] * S[2] - S[1] ** 2 > 0
lo = gstar.lower(); hi = pi.upper()
stack = [(lo + (hi - lo) * i / 100, lo + (hi - lo) * (i + 1) / 100, 0) for i in range(100)]
used = 0; worst = [None, None]
while stack and ok:
    g0, g1, dep = stack.pop()
    G = arb((g0 + g1) / 2, (g1 - g0) / 2)
    So, Se = sums(G)
    if psd(So) and psd(Se):
        used += 1
        for j, S in enumerate((So, Se)):
            m = (S[0] * S[2] - S[1] ** 2).lower()
            if worst[j] is None or m < worst[j][0]:
                worst[j] = (m, float(G.mid()))
        continue
    if dep >= 28:
        ok = False; print('  FAIL at g =', G.str(10), So, Se); break
    m = (g0 + g1) / 2
    stack.append((g0, m, dep + 1)); stack.append((m, g1, dep + 1))
print('(4) X* = %.6f, g* = %.6f, K = %d, cells %d; min det lower bounds: S_odd %s, S_even+C/2 %s' % (float(Xstar), float(gstar.mid()), K, used, worst[0], worst[1]))
print('CERTIFIED per-gap validity (J=%d):' % JJ, bool(ok))
json.dump({'a': [str(x) for x in a], 'phi': {k: str(v) for k, v in ph.items()}, 'D0': str(D0), 'E12': str(E12), 'E32': str(E32),
           'Xstar': float(Xstar), 'cells': used, 'K': K, 'ok': bool(ok)}, open('logs/cert5_%s.json' % tag, 'w'), indent=1)
