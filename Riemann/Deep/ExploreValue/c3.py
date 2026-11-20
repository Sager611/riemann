# C3: Nevanlinna's second main theorem in Littlewood form on the half-strip [1/2, 1/2+X/L] x [T1,T2] for
# g = (F'/F)/(chi'/chi) (Re g = 1/2 on the line; zeros of V_c = (1/2+c)-points of g).  For q values a_j and scale lam:
#  2pi[sum_j W(g=a_j) - W(g'=0) - (q-2) W(g=inf)] >= RHS := sum_j I+(g-a_j) - I+(lam g') - I(log+ sum_j |lam g'/(g-a_j)|)
#     - (q-1) H log+(2/delta) - sum_j J(g-a_j) + J(lam g')      (I = line integral, J = far-side integral; edges O(1))
# W = Littlewood weight sum of min(beta-1/2, X/L).  Bad values Re a > 1/2: W(g=a) counts off-line zeros only.
import numpy as np, math, sys, itertools
from fz import Fn, line_eval
kind = sys.argv[1]; T1 = float(sys.argv[2]); T2 = float(sys.argv[3]); dt = float(sys.argv[4])
f = Fn(kind, T2 + 20); L = math.log((T1+T2)/2/(2*math.pi)); H = T2 - T1
def gdata(sig, step):
    n = int(round(H/step)) + 1; t = T1 + step*np.arange(n); F, Fp = line_eval(f, sig, T1, n, step, Bk=50)
    g = (Fp/F)/f.chilog(sig + 1j*t); gp = -1j*np.gradient(g, step); return t, g, gp
t0, g0, gp0 = gdata(0.5, dt)
far = {X: gdata(0.5 + X/L, 0.01) for X in (0.5, 1, 2, 4)}
I = lambda y, t: np.trapezoid(y, t)
lp = lambda x: np.maximum(np.log(np.abs(x)), 0)
bad = [0.5 + c for c in (0.1, 0.25, 0.5, 1, 2)] + [0.5 + 0.5*np.exp(1j*p) for p in (-1.2, -0.6, 0.6, 1.2)]
good = [0.5 - c for c in (0.1, 0.25, 0.4)] + [0.5 - 0.3*np.exp(1j*p) for p in (-1.0, 1.0)] + [-1.0, 0.0]
best = {}
for label, vals in (('bad', bad), ('good', good), ('mixed', bad[:3] + good[:3])):
    for q in (2, 3, 4, 6):
        for A in itertools.islice(itertools.combinations(vals, q), 60):
            A = np.array(A, complex); dl = min(abs(x - y) for x, y in itertools.combinations(A, 2))
            for lam in (0.25/L, 1/L, 4/L):
                lin = sum(I(lp(g0 - a), t0) for a in A) - I(lp(lam*gp0), t0) - I(lp(sum(np.abs(lam*gp0/(g0 - a)) for a in A)), t0) - (q-1)*H*max(0, math.log(2/dl))
                for X, (tX, gX, gpX) in far.items():
                    rhs = lin - sum(I(np.log(np.abs(gX - a)), tX) for a in A) + I(np.log(np.abs(lam*gpX)), tX)
                    k = (label, q)
                    if k not in best or rhs > best[k][0]: best[k] = (rhs, X, lam*L, np.round(A, 3))
# actual Littlewood weights from line integrals (identity), for comparison
for label, q in sorted(best):
    rhs, X, lamL, A = best[(label, q)]; tX, gX, gpX = far[X]
    Wa = [ (I(np.log(np.abs(g0 - a)), t0) - I(np.log(np.abs(gX - a)), tX))/(2*math.pi) for a in A ]   # = W(g=a) - W(g=inf)
    print(f"{kind} {label:5s} q={q}: best SMT lower bound RHS/2pi = {rhs/(2*math.pi):9.2f} (X={X}, lam*L={lamL}); actual sum_j [W(a_j)-W(inf)] = {sum(Wa):9.3f}; a={A}")
print(f"{kind}: H={H} L={L:.4f} N~{H*L/(2*math.pi):.0f}; slack scale H/2pi = {H/(2*math.pi):.1f}")
