# W(g') on R_X (Littlewood): 2pi W(g') = I(log|g'|) - J(log|g'|) (edges O(1)); W(g') = W(g'=0) - 2 W(g=inf). Same data as c3.py.
import numpy as np, math, sys
from fz import Fn, line_eval
kind = sys.argv[1]; T1 = float(sys.argv[2]); T2 = float(sys.argv[3]); dt = float(sys.argv[4])
f = Fn(kind, T2 + 20); L = math.log((T1+T2)/2/(2*math.pi)); H = T2 - T1
def gdata(sig, step):
    n = int(round(H/step)) + 1; t = T1 + step*np.arange(n); F, Fp = line_eval(f, sig, T1, n, step, Bk=50)
    g = (Fp/F)/f.chilog(sig + 1j*t); gp = -1j*np.gradient(g, step); return t, g, gp
t0, g0, gp0 = gdata(0.5, dt); I0 = np.trapezoid(np.log(np.abs(gp0)), t0)
for X in (0.5, 1, 2, 4):
    tX, gX, gpX = gdata(0.5 + X/L, 0.01)
    print(f"{kind} X={X}: W(g') = {(I0 - np.trapezoid(np.log(np.abs(gpX)), tX))/(2*math.pi):.3f}", flush=True)
