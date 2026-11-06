# t1_lattice.py (explore21): grid-code optimality of 2D lattices vs the zeros of their Epstein zeta functions.
# A) table: packing density, theta energy theta(1), Epstein Z(2) at covolume 1, class number, Euler product or not.
# B) D = -15 non-principal form Q2 = 2x^2+xy+2y^2: Z_Q2 = zeta*L(chi_-15) - L(chi_-3)*L(chi_5) (checked at s = 3 by a direct lattice sum),
#    zeros: on-line sign changes vs the argument-principle total; off-line zeros located by minima of |Z_Q2| and Newton.
# Usage: t1_lattice.py A | B T
import sys, math, time
import numpy as np
from flint import acb, arb, dirichlet_char, ctx
import mpmath as mp

def classno(D):
    # reduced primitive forms (a,b,c), b^2-4ac = D < 0
    h = 0; forms = []
    a = 1
    while 3*a*a <= -D:
        for b in range(-a+1, a+1):
            if (b*b - D) % (4*a) == 0:
                c = (b*b - D)//(4*a)
                if c < a: continue
                if b < 0 and (a == c): continue
                if math.gcd(math.gcd(a, abs(b)), c) != 1: continue
                forms.append((a, b, c)); h += 1
        a += 1
    return h, forms

def lattice_vals(a, b, c, R=60):
    m, n = np.meshgrid(np.arange(-R, R+1), np.arange(-R, R+1))
    q = (a*m*m + b*m*n + c*n*n).astype(float).ravel()
    return q[q > 0]

def table():
    forms = [((1,1,1),'hexagonal'), ((1,0,1),'square'), ((1,1,2),''), ((1,0,2),''), ((1,1,3),''), ((1,0,3),'D=-12'),
             ((1,1,4),'principal'), ((2,1,2),'non-principal'), ((1,0,5),'principal'), ((2,2,3),'non-principal'),
             ((1,0,6),'principal'), ((2,0,3),'non-principal'), ((1,1,6),'principal'), ((2,1,3),'non-principal'),
             ((1,1,5),''), ((1,1,41),'D=-163')]
    print('form       D    h  Euler  pack.dens  theta(1)    Z(2)   [covolume 1; hexagonal optimal: max density, min theta, min Z]')
    for (a,b,c), note in forms:
        D = b*b - 4*a*c; h, _ = classno(D)
        cov = math.sqrt(-D)/2
        R = 60 if -D < 30 else 200
        q = lattice_vals(a, b, c, R)/cov
        dens = math.pi*q.min()/4
        th = 1 + np.exp(-math.pi*q).sum()
        # Z(2): direct sum inside the largest disc contained in the box + tail pi/r^2 (covolume 1, s = 2)
        mb = np.arange(-R, R+1)
        edge = np.concatenate([a*R*R + b*R*mb + c*mb*mb, a*mb*mb + b*mb*R + c*R*R]).astype(float)/cov
        rc2 = edge.min()
        Z2 = (q[q < rc2]**-2.0).sum() + math.pi/rc2
        print('%-10s %4d %3d  %-5s  %.4f    %.6f  %.5f  %s' % ((a,b,c), D, h, 'yes' if h == 1 else 'no', dens, th, Z2, note))

ctx.prec = 64
chi15 = dirichlet_char(15, 14); chi3 = dirichlet_char(3, 2); chi5 = dirichlet_char(5, 4)
assert [int(round(float(chi15(n).real))) for n in range(1, 16)] == [1,1,0,1,0,0,-1,1,0,0,-1,0,-1,-1,0]
assert [int(round(float(chi3(n).real))) for n in range(1, 4)] == [1,-1,0]
assert [int(round(float(chi5(n).real))) for n in range(1, 6)] == [1,-1,-1,1,0]

def ZQ2(s):
    s = acb(s)
    return s.zeta()*s.dirichlet_l(chi15) - s.dirichlet_l(chi3)*s.dirichlet_l(chi5)

def ZQ1(s):
    s = acb(s)
    return s.zeta()*s.dirichlet_l(chi15) + s.dirichlet_l(chi3)*s.dirichlet_l(chi5)

def Zhex(s):
    s = acb(s)
    return 6*s.zeta()*s.dirichlet_l(chi3)

KIND = {'Q2': (ZQ2, 15, 2), 'Q1': (ZQ1, 15, 1), 'hex': (Zhex, 3, 1)}
CUR = ['Q2']
def Lam(s, Z=None):
    Zf, D, n0 = KIND[CUR[0]]
    s = acb(s)
    return (acb(arb(D).sqrt()/(2*arb.pi()))**s)*s.gamma()*Zf(s)

def tocomplex(z): return complex(float(z.real.mid()), float(z.imag.mid()))

def argchange(pts, Z=ZQ2, maxjump=0.8):
    # robust: recursive bisection on each segment
    def val(z): return tocomplex(Lam(z, Z))
    def seg(a, fa, b, fb, depth=0):
        d = np.angle(fb/fa)
        if abs(d) <= maxjump or depth > 30: return d
        m = (a + b)/2; fm = val(m)
        return seg(a, fa, m, fm, depth+1) + seg(m, fm, b, fb, depth+1)
    tot = 0.0; fa = val(pts[0])
    for a, b in zip(pts[:-1], pts[1:]):
        fb = val(b); tot += seg(a, fa, b, fb); fa = fb
    return tot

def hardy_real(t, Z=ZQ2):
    s = acb(0.5, t)
    v = Lam(s, Z)
    return float(v.real.mid()), float(v.imag.mid())

def partB(T, kind='Q2'):
    CUR[0] = kind; Zf, Dk, n0 = KIND[kind]
    t0 = time.time()
    # check decomposition at s = 3 by a direct sum
    q = lattice_vals(2, 1, 2, 300)
    mb = np.arange(-300, 301)
    edge = np.concatenate([2*300*300 + 300*mb + 2*mb*mb, 2*mb*mb + mb*300 + 2*300*300]).astype(float); rc = edge.min()
    direct = (q[q < rc]**-3.0).sum() + (2*math.pi/math.sqrt(15))/(2*rc**2)  # tail: (2 pi/sqrt|D|) int_rc^inf u^-3 du (lattice-point density 2/sqrt|D|)
    if kind == 'Q2': print('check s=3: direct %.12f  L-formula %.12f' % (direct, float(ZQ2(3).real.mid())), flush=True)
    ta = 0.5
    # total zero count of Lam in [-2,3] x [ta, T] (both sides of the line)
    nb = np.linspace(0.5, 3, 60); bottom = [complex(x, ta) for x in nb]
    right = [complex(3, y) for y in np.arange(ta, T + 1e-9, 0.25)]
    if right[-1].imag < T: right.append(complex(3, T))
    top = [complex(x, T) for x in np.linspace(3, 0.5, 120)]
    D = argchange(bottom) + argchange(right) + argchange(top)
    Ntot = D/math.pi
    # on-line zeros: sign changes of the real Lam(1/2+it)
    ts = np.arange(ta, T, 0.01)
    vals = np.array([hardy_real(t)[0] for t in ts])
    imax = max(abs(hardy_real(t)[1]) for t in ts[::500])
    sc = np.where(np.sign(vals[:-1]) != np.sign(vals[1:]))[0]
    print('%s T=%g  total zeros (arg principle) %.3f   on-line sign changes %d   max|Im Lam| on line %.2e   smooth main term (T/pi)log(sqrt|D| T/(2 pi e n0)) = %.2f   %.0fs'
          % (kind, T, Ntot, len(sc), imax, T/math.pi*math.log(math.sqrt(Dk)*T/(2*n0*math.pi*math.e)), time.time()-t0), flush=True)
    # off-line candidates: local minima of |Z| at sigma = 0.55..1.0
    cands = []
    for sg in (0.56, 0.62, 0.7, 0.8, 0.9, 1.0):
        tt = np.arange(ta, T, 0.02)
        av = np.array([abs(tocomplex(Zf(acb(sg, t)))) for t in tt])
        mins = [i for i in range(1, len(av)-1) if av[i] < av[i-1] and av[i] < av[i+1] and av[i] < 0.15]
        for i in mins: cands.append((sg, tt[i], av[i]))
    print('candidates', len(cands), flush=True)
    mp.mp.dps = 20
    def f(s):
        z = Zf(acb(str(mp.re(s)), str(mp.im(s))))
        return mp.mpc(float(z.real.mid()), float(z.imag.mid()))
    roots = []
    for sg, t, a in cands:
        try:
            r = mp.findroot(f, mp.mpc(sg, t), tol=1e-12, maxsteps=50)
        except Exception:
            continue
        if mp.re(r) > 0.505 and abs(f(r)) < 1e-8 and not any(abs(r - x) < 1e-6 for x in roots):
            roots.append(r)
    roots.sort(key=lambda z: float(mp.im(z)))
    for r in roots: print('off-line zero  sigma = %.10f  t = %.10f' % (float(mp.re(r)), float(mp.im(r))), flush=True)
    print('done %.0fs' % (time.time()-t0))

if __name__ == '__main__':
    if sys.argv[1] == 'A': table()
    else: partB(float(sys.argv[2]), sys.argv[3] if len(sys.argv) > 3 else 'Q2')
