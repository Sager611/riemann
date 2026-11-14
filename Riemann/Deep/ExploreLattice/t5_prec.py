# t5_prec.py (explore48): precision regression test after the radius checks (y = 40 at t = 176 must give |beta| tiny).
import sys, math, time
sys.path.insert(0, '.')
from epstein_cs import F_line, Estar, Kscaled
for (t, u) in ((176.0, 2 * math.pi * 40), (130.0, 2 * math.pi * 30), (12.6, 2 * math.pi * 10)):
    t0 = time.time(); print('t=%g u=%.2f  e^{pi t/2}K=%.3e  [%.3fs]' % (t, u, Kscaled(t, u), time.time() - t0))
for t in (170.0, 176.0, 200.0, 240.0, 250.0, 255.0):
    Fv, c, b = F_line(t, 0.0, 40.0, return_parts=True)[:3]
    print('y=40 t=%g beta=%.3e' % (t, b))
print(Estar(complex(0.6, 176.0), 0.0, 40.0))
