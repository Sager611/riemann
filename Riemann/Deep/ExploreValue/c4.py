# C4: S1 for complex values a and for a = 0 (zeros of F', Speiser / Levinson-Montgomery): count zeros of F' - a (chi'/chi) F
# (= a-points of g = (F'/F)/(chi'/chi)) in [1/2,3] x [T1,T2]. Prediction: N+ (Re a > 1/2), N+ + N0bar (Re a < 1/2, a != 0),
# N+ + N0bar - (log 2/2pi) H (a = 0: far side F'/F ~ -c 2^{-s} winds).
import numpy as np, math, sys
from fz import Fn, line_eval
kind = sys.argv[1]; T1 = float(sys.argv[2]); T2 = float(sys.argv[3]); dt = float(sys.argv[4])
f = Fn(kind, T2 + 20); n = int(round((T2 - T1)/dt)) + 1; t = T1 + dt*np.arange(n)
F, Fp = line_eval(f, 0.5, T1, n, dt, Bk=50); ch = f.chilog(0.5 + 1j*t)
Z = np.real(np.exp(1j*f.theta(t))*F); n0 = int(np.sum(np.sign(Z[1:]) != np.sign(Z[:-1])))
sg = np.arange(0.5, 3.0, 0.0005)
Fb, Fpb = f(sg + 1j*T1, blk=40); Ft, Fpt = f(sg + 1j*T2, blk=40); chb = f.chilog(sg + 1j*T1); cht = f.chilog(sg + 1j*T2)
nr = int(round((T2 - T1)/0.02)) + 1; tr = T1 + 0.02*np.arange(nr); Fr, Fpr = line_eval(f, 3.0, T1, nr, 0.02); chr_ = f.chilog(3.0 + 1j*tr)
Ntot = (f.theta(np.array([T2]))[0] - f.theta(np.array([T1]))[0])/math.pi + (np.unwrap(np.angle(Ft[::-1]))[-1] - np.unwrap(np.angle(Fb[::-1]))[-1])/math.pi
npl = (Ntot - n0)/2
def wind(G_b, G_r, G_t, G_l):
    seq = np.concatenate([G_b, G_r, G_t[::-1], G_l[::-1], G_b[:1]]); ph = np.unwrap(np.angle(seq))
    return (ph[-1] - ph[0])/(2*math.pi), np.max(np.abs(np.diff(ph)))
print(f"{kind} [{T1},{T2}]: N0bar={n0} Ntot={Ntot:.2f} N+={npl:.2f}; (log2/2pi)H = {math.log(2)/(2*math.pi)*(T2-T1):.2f}")
for a in [0.5+0.3j, 0.5-0.7j, 0.8+0.1j, 2+2j, 0.45+0j, 0.2+0.5j, -1-1j, 0.1-0.3j, 0j]:
    V = lambda F_, Fp_, ch_: Fp_ - a*ch_*F_
    w, st = wind(V(Fb, Fpb, chb), V(Fr, Fpr, chr_), V(Ft, Fpt, cht), V(F, Fp, ch))
    pred = npl if a.real > 0.5 else (npl + n0 if a != 0 else npl + n0 - math.log(2)/(2*math.pi)*(T2-T1))
    print(f"  a={a}: #a-points right of 1/2 = {w:.2f} (max phase step {st:.2f}); predicted {pred:.2f}")
