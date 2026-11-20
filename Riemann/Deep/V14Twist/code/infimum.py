# Infimum of the scalar system with the improved count window (line or general gain), b re-optimized.
import numpy as np, sys
from system import *
DOM = ((11 / 50, 13 / 20), (6 / 25, 0.5))
def Bstar(b, Rfun, nd=1201, nx=61):
    d = np.linspace(*DOM[0], nd); x = np.linspace(*DOM[1], nx)
    dd, xx = np.meshgrid(d, x)
    R = Rfun(dd, xx)
    Bv = Broot(b, dd, xx, R)
    i = np.unravel_index(np.argmax(Bv), Bv.shape)
    return float(Bv[i]), float(dd[i]), float(xx[i])
def golden(f, a, b, tol=1e-12):
    g = (5 ** 0.5 - 1) / 2
    c, d_ = b - g * (b - a), a + g * (b - a)
    fc, fd = f(c), f(d_)
    while b - a > tol:
        if fc < fd: b, d_, fd = d_, c, fc; c = b - g * (b - a); fc = f(c)
        else: a, c, fc = c, d_, fd; d_ = a + g * (b - a); fd = f(d_)
    return (a + b) / 2
def Rnat(dd, xx): return Rd_native(dd, xx)
def mk_line(g0, g1, win=WIN):
    def Rf(dd, xx):
        v, c = R_line(dd, xx, g0, g1)
        inwin = (dd >= win[0]) & (dd <= win[1]) & (xx >= win[2]) & (xx <= win[3])
        return np.where(inwin, np.minimum(v, Rd_native(dd, xx)), Rd_native(dd, xx))
    return Rf
def report(name, Rf, blo=0.05, bhi=0.16):
    bopt = golden(lambda b: Bstar(b, Rf)[0], blo, bhi, 1e-9)
    B, d, x = Bstar(bopt, Rf, nd=4001, nx=121)
    l = 11 / 3 - 4 * B; h = (1 + bopt + 3 * l) / 2
    print("%-34s B*=%.12f  b=%.6f h=%.6f  contact delta=%.5f x=%.4f" % (name, B, bopt, h, d, x)); sys.stdout.flush()
    return B, bopt, h, d, x
if __name__ == '__main__':
    Bn = report("native (PR103 system)", Rnat)
    Bl = report("PR137 line, PR137 window", mk_line(float(G0), float(G1)))
    Bw = report("PR137 line, window = all DOM", mk_line(float(G0), float(G1), (0, 1, 0, 1)))
    # PR137's own parameter point: theta=2183367/2500000, h=81447/100000
    th = 2183367 / 2500000; l = 11 / 3 - 4 * th; bh = 2 * 0.81447 - 1 - 3 * l
    print("PR137 point: b(h=0.81447)=%.7f  B*(b)=%.12f (theta=%.7f)" % (bh, Bstar(bh, mk_line(float(G0), float(G1)), 4001, 121)[0], th))
