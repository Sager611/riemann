# ceilaudit9 E30: is explore30's threshold (adversary admissible only for a >= 1.7, c = e^{2a} >= 30 at T = 3e12; a >= 1.925 / c >= 47 at every
# height) a property of its placement family ("full gaps") or of the test-function class? Periodic block of zeros, period P in t; W/(D 2pi||f||^2) per period.
# (1) single harmonic (pi/a < P <= 2pi/a): |g(y)| <= g(0)/2 for |y| >= a (Cauchy-Schwarz on disjoint supports), so two bumps at +-a' are optimal:
#     W_norm = 1 + (1/n) sum_rho cosh(2a' d_rho) cos(phi_rho), n = D P; best placement = LP over the phase distribution under |S| <= Smax, |S_1| <= S1max.
# (2) multi-harmonic f = sum_k c_k bump(x - k s), k = -K..K, s = a'/K, period 2pi/s: alternate min-eigenvector in c and LP in the placement.
import sys, math, numpy as np
from scipy.optimize import linprog
def lp(u, n, D, P, Smax, S1max, J, anchor=False):
    # variables mu_0..mu_{J-1} >= 0 (atoms at t_j = j*P/J), s0, s1 free. maximise sum mu_j u_j.
    dt = P/J; tj = dt*np.arange(J); nv = J + 2
    A = []; bnd = []
    cum = np.tril(np.ones((J, J)))            # S_j^+ = sum_{i<=j} mu_i - D t_j - s0
    for j in range(J):
        row = np.zeros(nv); row[:J] = cum[j]; row[J] = -1.0
        A.append(row.copy()); bnd.append(Smax + D*tj[j])          # S_j^+ <= Smax
        A.append(-row); bnd.append(Smax - D*tj[j])                # -S_j^+ <= Smax
        row2 = row.copy(); row2[j] -= 1.0                         # S_j^- = S_j^+ - mu_j
        A.append(row2.copy()); bnd.append(Smax + D*tj[j]); A.append(-row2); bnd.append(Smax - D*tj[j])
    # integral of S over [t_j, t_{j+1}) = S_j^+ dt - D dt^2/2 ; S1(t_j) = sum_{i<j} (...) - s1
    Icoef = np.zeros((J+1, nv)); Iconst = np.zeros(J+1)
    for j in range(J):
        Icoef[j+1] = Icoef[j]; Icoef[j+1, :J] += cum[j]*dt; Icoef[j+1, J] += -dt
        Iconst[j+1] = Iconst[j] - D*tj[j]*dt - D*dt*dt/2
    for j in range(J+1):
        row = Icoef[j].copy(); row[J+1] = -1.0
        A.append(row.copy()); bnd.append(S1max - Iconst[j]); A.append(-row); bnd.append(S1max + Iconst[j])
    Aeq = [np.concatenate([np.ones(J), [0, 0]]), Icoef[J].copy()]; beq = [n, -Iconst[J]]   # count, zero-mean S (periodic S_1)
    c = np.concatenate([-u, [0, 0]])
    r = linprog(c, A_ub=np.array(A), b_ub=np.array(bnd), A_eq=np.array(Aeq), b_eq=np.array(beq),
                bounds=[(0, None)]*J + ([(0, 0)]*2 if anchor else [(None, None)]*2), method='highs')   # anchor: S = S_1 = 0 at the period start
    return r
def single(a, wf, D, Smax, S1max, dmax=0.37, J=240):
    ap = a*(1 - wf); P = math.pi/ap; n = D*P
    phi = 2*ap*(np.arange(J)*P/J) + math.pi             # theta = pi: phase 0 of the period = most negative K
    ch = math.cosh(2*ap*dmax)
    u = np.where(np.cos(phi) < 0, -ch*np.cos(phi), -np.cos(phi))   # gain per zero (off-line where cos < 0)
    r = lp(u, n, D, P, Smax, S1max, J)
    return (1 - (-r.fun)/n) if r.status == 0 else float('nan')
def single_anchored(a, wf, D, Smax, S1max, dmax=0.37, J=240, npsi=32):
    # block starts and ends where S = S_1 = 0 (explore30's ends), phase of the kernel at the block start searched over npsi values
    ap = a*(1 - wf); P = math.pi/ap; n = D*P; ch = math.cosh(2*ap*dmax); best = float('nan')
    for psi in np.arange(npsi)*2*math.pi/npsi:
        phi = 2*ap*(np.arange(J)*P/J) + psi
        u = np.where(np.cos(phi) < 0, -ch*np.cos(phi), -np.cos(phi))
        r = lp(u, n, D, P, Smax, S1max, J, anchor=True)
        if r.status == 0:
            w_ = 1 - (-r.fun)/n
            if not (w_ >= best): best = w_
    return best
def bounds(t):
    lt = math.log(t); D = math.log(t/(2*math.pi))/(2*math.pi)
    Smax = min(0.112*lt + 0.278*math.log(lt) + 2.510, 0.0946*lt + 5.746) - 1.0; S1max = 0.059*lt + 2.067 - 1.0
    return D, Smax, S1max
if __name__ == '__main__':
    D, Sm, S1m = bounds(3e12 + 1e4)
    print(f"T=3e12: D={D:.4f} Smax={Sm:.3f} S1max={S1m:.3f}", flush=True)
    for wf in (0.03, 0.0):
        print(f"  single harmonic, LP placement, w = {wf}a: " + " ".join(f"a={a:.2f}(c={math.exp(2*a):.1f}):{single(a, wf, D, Sm, S1m):+.3f}" for a in (1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7)), flush=True)
    # asymptotic (t -> infinity): scale by D: Smax/D -> 0.0946*2pi, S1max/D -> 0.059*2pi
    print(f"t->inf: Smax/D={0.0946*2*math.pi:.4f} S1max/D={0.059*2*math.pi:.4f}", flush=True)
    for wf in (0.03, 0.0):
        print(f"  single harmonic, LP placement, w = {wf}a: " + " ".join(f"a={a:.2f}(c={math.exp(2*a):.1f}):{single(a, wf, 1.0, 0.0946*2*math.pi, 0.059*2*math.pi):+.3f}" for a in (1.3, 1.4, 1.5, 1.6, 1.7, 1.78, 1.925)), flush=True)
    for t in (1e20, 1e100):
        D, Sm, S1m = bounds(t)
        print(f"  t={t:.0e} (D={D:.2f}, Smax={Sm:.2f}, S1max={S1m:.2f}), w=0.03a: " + " ".join(f"a={a:.2f}:{single(a, 0.03, D, Sm, S1m):+.3f}" for a in (1.3, 1.4, 1.5, 1.6, 1.7, 1.8)), flush=True)
