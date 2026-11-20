"""C4 extra configurations (same model and units as c4_heal.py): multiplicity, one-sided gap, spread deficit, small beta/Delta (DH-like)."""
from c4_heal import *
from fractions import Fraction as Fr
import math, time
M = 40; H = Fr(742, 100)
th = 2*math.atan(float(BETA/G0)); dA = 1/math.sin(th)**2
cfg = {
  'F double pair isolated': ([], [G0, G0], BETA),
  'G double pair in S(T)-worst gap (4 zeros: half-width 8.42), real beyond': ([G0 + s*(H + 1 + m)*DEL for s in (-1, 1) for m in range(M)], [G0, G0], BETA),
  'H one-sided gap 14.84 Delta above, real below at spacing Delta': ([G0 - m*DEL for m in range(1, M+1)] + [G0 + (2*H + m)*DEL for m in range(M)], [G0], BETA),
  'I deficit spread: real zeros at spacing 2 Delta, |h| <= 25.7 Delta': ([G0 + 2*m*DEL for m in range(-13, 14) if m != 0] + [G0 + s*(Fr(257, 10) + m)*DEL for s in (-1, 1) for m in range(20)], [G0], BETA),
}
for name, (rh, ph, b) in cfg.items():
    t0 = time.time(); g = build(rh, ph, b)
    dm = dmax(g, int(40*dA))
    print('%-72s zeros %3d  d_max/d_iso = %.3f  (%.0f s)' % (name, len(rh)+2*len(ph), dm/dA, time.time()-t0), flush=True)
b = Fr(21, 100)*DEL; th2 = 2*math.atan(float(b/G0)); dA2 = 1/math.sin(th2)**2
for name, rh in (('J beta/Delta = 0.21, isolated', []), ('K beta/Delta = 0.21, real neighbours both sides', [G0 + m*DEL for m in range(-M, M+1) if m != 0])):
    t0 = time.time(); g = build(rh, [G0], b); dm = dmax(g, int(40*dA2))
    print('%-72s zeros %3d  d_max/d_iso = %.3f  (%.0f s)' % (name, len(rh)+2, dm/dA2, time.time()-t0), flush=True)
