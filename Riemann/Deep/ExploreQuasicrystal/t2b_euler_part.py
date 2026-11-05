# t2b_euler_part.py -- zeros of the DH "Euler part" Q_X(s) = c1 P_X(s,chi) + c2 P_X(s,chibar) with Re s > 1/2
# (for zeta the Euler part P_X = exp(...) never vanishes).  Compare with half the off-line zeros of the DH model in t2_out.txt.
import numpy as np, math, sys
sys.argv = ['x']
src = open('t2_models.py').read().split("t0 = time.time()")[0]
g = {}; exec(compile(src, 't2_models.py', 'exec'), g)
def wind_Q(pp, s1, s2, T1, T2, dt=0.002, ds=0.0005):
    bot = s1 + np.arange(0, s2 - s1, ds) + 1j * T1; rig = s2 + 1j * (T1 + np.arange(0, T2 - T1, dt))
    top = s2 - np.arange(0, s2 - s1, ds) + 1j * T2; lef = s1 + 1j * (T2 - np.arange(0, T2 - T1, dt))
    z = np.concatenate([bot, rig, top, lef, [bot[0]]])
    Dc = g['dirpoly'](z, pp, 'c'); Dcb = np.conj(g['dirpoly'](np.conj(z), pp, 'c'))
    m = np.maximum(Dc.real, Dcb.real); v = g['C1'] * np.exp(Dc - m) + g['C2'] * np.exp(Dcb - m)
    ph = np.unwrap(np.angle(v)); return (ph[-1] - ph[0]) / (2 * math.pi), np.max(np.abs(np.diff(ph)))
print("T       X   zeros of Q_X in 0.503<=sigma<=0.98   zeros in 0.503<=sigma<=0.70   [max steps]")
for T in [1000.0, 10000.0, 50000.0, 70000.0]:
    for X in [5, 10, 30]:
        pp = g['prime_powers'](X)
        a, ja = wind_Q(pp, 0.503, 0.98, T, T + 200.0); b, jb = wind_Q(pp, 0.503, 0.70, T, T + 200.0)
        print("%7g %3d   %8.3f   %8.3f   [%.2f %.2f]" % (T, X, a, b, ja, jb)); sys.stdout.flush()
