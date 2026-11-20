"""C4 adversarial/random configurations (model and units of c4_heal.py). Bad pair at G0, offset BETA = 1.6 Delta.
D3: S(T)-worst gap, zeros beyond doubled (multiplicity 2 at spacing 2 Delta, same count);
D5: same gap, bad pair off-centre (3 Delta from the gap middle);
R1..R8 (seed 65): Wigner-surmise spacings (mean Delta), random gap half-width in [0, 7.42] Delta around the pair,
each neighbour off-line with probability 1/2 at a uniform offset in [0, BETA]; 20 zeros per side."""
from c4_heal import *
from fractions import Fraction as Fr
import math, time, random
M = 30; H = Fr(742, 100)
th = 2*math.atan(float(BETA/G0)); dA = 1/math.sin(th)**2
def run(name, rh, ph):
    t0 = time.time(); g = build(rh, ph); dm = dmax(g, int(40*dA))
    print('%-62s zeros %3d  d_max/d_iso = %.3f  (%.0f s)' % (name, len(rh)+2*len(ph), dm/dA, time.time()-t0), flush=True)
    return dm/dA
res = []
dbl = [G0 + s*(H + 2*m)*DEL for s in (-1, 1) for m in range(M//2)]
res.append(run('D3 worst gap, zeros beyond doubled', dbl + dbl, [G0]))
res.append(run('D5 worst gap, pair 3 Delta off-centre', [G0 + Fr(3)*DEL + s*(H + m)*DEL for s in (-1, 1) for m in range(M)], [G0]))
rng = random.Random(65)
def wigner():  # GUE-like spacing, mean 1: density (32/pi^2) s^2 exp(-4 s^2/pi), sampled by rejection
    while True:
        s = rng.uniform(0, 4); 
        if rng.uniform(0, 1.1) < (32/math.pi**2)*s*s*math.exp(-4*s*s/math.pi): return s
for k in range(8):
    h = Fr(rng.randint(0, 742), 100)
    rh, ph = [], [G0]
    for sgn in (-1, 1):
        pos = h
        for m in range(20):
            pos += Fr(int(1000*wigner()), 1000) if m else 0
            z = G0 + sgn*(pos + Fr(1, 2))*DEL
            if rng.random() < 0.5: rh.append(z)
            else:
                b = BETA*Fr(rng.randint(0, 100), 100)
                ph.append((z, b))
    # build with per-pair offsets
    P = [Fr(1)]
    for gz in rh: P = poly_mul(P, [Fr(1), 1/(gz*gz)])
    for item in ph:
        gz, b = (item, BETA) if not isinstance(item, tuple) else item
        re = gz*gz - b*b; im = 2*gz*b; ab2 = re*re + im*im
        P = poly_mul(P, [Fr(1), 2*re/ab2, 1/ab2])
    g = [P[j]*math.factorial(j) for j in range(len(P))]
    t0 = time.time(); dm = dmax(g, int(40*dA))
    print('R%d gap half-width %.2f Delta, %d real + %d off-line pairs            d_max/d_iso = %.3f  (%.0f s)' % (k+1, float(h), len(rh), len(ph)-1, dm/dA, time.time()-t0), flush=True)
    res.append(dm/dA)
print('min ratio over these configurations: %.3f' % min(res))
