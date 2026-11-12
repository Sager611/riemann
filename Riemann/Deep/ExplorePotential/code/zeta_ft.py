# zeta_ft.py (explore35): Fourier coefficients of {alpha*gamma} mod 1 for zeta zeros vs Landau's prime-power prediction,
# and the Erdos-Turan bound vs the measured discrepancy.  Input: zeros (one per line or .npy), T = max ordinate used.
import sys, json, math, numpy as np
def Lam(n):
    # von Mangoldt of an integer-valued x (x within 1e-6 of an integer, x <= 1e9); 0 otherwise
    if n < 2 or n > 1e9 or abs(n - round(n)) > 1e-6: return 0.0
    n = int(round(n))
    for p in range(2, n + 1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return math.log(p) if m == 1 else 0.0
    return 0.0
def load(path):
    if path.endswith('.npy'): return np.load(path)
    return np.array([float(l.split()[0]) for l in open(path) if l.strip() and not l.startswith('#')])
def ft(g, x):          # sum_gamma x^{i gamma}
    ph = np.mod(g*math.log(x), 2*math.pi)
    return complex(np.cos(ph).sum(), np.sin(ph).sum())
def disc(u):           # interval discrepancy (counts) of points u in [0,1)
    u = np.sort(np.mod(u, 1.0)); N = len(u); i = np.arange(1, N + 1)
    return float(np.max(i - N*u) + np.max(N*u - (i - 1)))
if __name__ == '__main__':
    g = load(sys.argv[1]); T0 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
    g = g[g > T0]; N = len(g); Tlo, Thi = T0, float(g.max())
    H = Thi - Tlo
    out = {'N': N, 'Tlo': Tlo, 'Thi': Thi}
    rows = []
    for x in list(range(2, 41)) + [math.e, math.pi, 1.5, 2.5, 6.5]:
        S = ft(g, x); P = -(H/(2*math.pi))*Lam(x)/math.sqrt(x)   # Landau main term (RH: x^{i gamma} = x^{rho - 1/2})
        rows.append((x, S.real/N, S.imag/N, P/N, (S.real - P)/math.sqrt(N), abs(S - P)))
    out['ft'] = rows
    print(f'N = {N}, T in ({Tlo:.1f}, {Thi:.1f}], H/2pi = {H/2/math.pi:.1f}, 1/sqrt(N) = {1/math.sqrt(N):.4f}')
    print('| x | Re S/N | Im S/N | Landau -H Lam(x)/(2 pi sqrt(x) N) | (Re S - P)/sqrt(N) | abs(S - P) |')
    for r in rows:
        xs = ('%.4g' % r[0]) if not float(r[0]).is_integer() else str(int(r[0]))
        print('| %s | %+.5f | %+.5f | %+.5f | %+.2f | %.1f |' % ((xs,) + tuple(r[1:])))
    # Erdos-Turan vs measured discrepancy
    et = {}
    for name, alpha in [('log2/2pi', math.log(2)/(2*math.pi)), ('log6/2pi', math.log(6)/(2*math.pi)), ('sqrt2/10', math.sqrt(2)/10), ('1', 1.0)]:
        u = alpha*g; D = disc(u)
        best = None; ck = []
        for K in range(1, 200):
            x = math.exp(2*math.pi*K*alpha)
            if x > 1e300: break
            ck.append(abs(ft(g, x)))
            b = N/(K + 1) + 3*sum(c/(k + 1) for k, c in enumerate(ck))
            if best is None or b < best[0]: best = (b, K)
        # first-order prediction of the discrepancy from Landau's coefficients alone (RH): Delta(u) = sum_k c_k (e(ku)-1)/(2 pi i k) + cc
        uu = np.linspace(0, 1, 4001); Dl = np.zeros_like(uu)
        for k in range(1, 60):
            x = math.exp(2*math.pi*k*alpha); c = -(H/(2*math.pi))*Lam(x)/math.sqrt(x)
            if c != 0: Dl += 2*(c*(np.sin(2*math.pi*k*uu))/(2*math.pi*k))
        Dpred = float(Dl.max() - Dl.min())
        et[name] = {'alpha': alpha, 'D_measured': D, 'D/N': D/N, 'ET_bound': best[0], 'ET_K': best[1], 'D_Landau_first_order': Dpred}
        print(f'alpha = {name}: D = {D:.1f} (D/N = {D/N:.5f}); Erdos-Turan bound {best[0]:.1f} at K = {best[1]}; Landau first-order D = {Dpred:.1f}; sqrt(N) = {math.sqrt(N):.1f}')
    out['et'] = et
    if len(sys.argv) > 3: json.dump(out, open(sys.argv[3], 'w'), indent=1)
