# t1_analyze.py — explore6: exponent table for mu(n) f(n) (Test A) and the digit-family transfer T1.
# Input: t1_sieve_out.json (from t1_sieve_twists.py). Output: printed tables (saved as t1_out.txt).
import json, math, numpy as np
d = json.load(open('t1_sieve_out.json'))
N = d['N']; LOG = int(round(math.log2(N)))
print(f"# Test A: S_f(y) = sum_(n<=y) mu(n) f(n), y < N = 2^{LOG} = {N}")
print("# m(K) = max_{2^K <= y < 2^(K+1)} |S_f(y)|; slope = least-squares d log2 m / dK over K = 12..%d" % (LOG-1))
print("%-16s %7s %7s | %s" % ('f', 'slope', 'slope2', 'm(K)/2^(K/2) for K = 14, 17, 20, 23, 25, 26'))
Ks = np.arange(12, LOG)
rows = {}
for k in d['names']:
    m = np.array(d['dmax'][k])
    y = np.log2(m[Ks])
    sl = np.polyfit(Ks, y, 1)[0]
    K2 = np.arange(18, LOG); sl2 = np.polyfit(K2, np.log2(m[K2]), 1)[0]
    norm = [m[K] / 2 ** (K / 2) for K in (14, 17, 20, 23, 25, 26)]
    rows[k] = (sl, sl2, norm)
    print("%-16s %7.3f %7.3f | %s" % (k, sl, sl2, ' '.join('%6.3f' % v for v in norm)))
# ---------------- T1: digit family F(alpha) = sum_{n<2^L} mu(n) e(alpha s_2(n)) = sum_k c_k z^k --------------
print()
print("# T1: digit-sum family; c_k(L) = sum_{n<2^L, s_2(n)=k} mu(n); F(alpha) = sum_k c_k e(k alpha); M = F(0)")
H = {int(L): np.array(v) for L, v in d['hist'].items()}
def binomw(L):
    return np.array([math.comb(L, k) for k in range(len(H[L]))], float) / 2.0**L
# Exact identity: mu(2m) = -mu(m) for odd m and s_2(2m) = s_2(m), so c_k(L) = sum over ODD m in [2^(L-1), 2^L) with s_2(m) = k.
# Random model: independent signs on odd squarefree m (density 8/pi^2), class size binom(L-2, k-2).
def oddtop(L, n):
    return np.array([math.comb(L - 2, k - 2) if k >= 2 else 0 for k in range(n)], float)
print("L   M(2^L)   rms z_k naive (6/pi^2*#class)   rms z_k exact-model (8/pi^2*#odd top-half class)   |M|/sqrt(2^L)")
for L in range(14, LOG + 1):
    c = H[L]; M = c.sum(); b = binomw(L); cnt = b * 2.0**L
    r = c - M * b; sel = cnt >= 1000
    z = r[sel] / np.sqrt(6 / math.pi**2 * cnt[sel])
    ot = oddtop(L, len(c)); sel2 = ot >= 1000
    z2 = r[sel2] / np.sqrt(8 / math.pi**2 * ot[sel2])
    print("%2d %9d   %.3f (%d classes)              %.3f (%d classes)                                   %.4f" % (L, M, math.sqrt((z**2).mean()), sel.sum(), math.sqrt((z2**2).mean()), sel2.sum(), abs(M) / 2**(L/2)))
# Bernstein-Walsh / two-constants on the arc A = {e^{i th}: |th| >= 2 pi delta}
def arc_nodes(delta, Q=400000):
    a0 = 2 * math.pi * delta; a = math.tan(a0 / 2)
    t = (np.arange(Q) + 0.5) / Q * math.pi - math.pi / 2       # harmonic measure at z=1 is dt/pi, t in (-pi/2, pi/2)
    zeta = -1.0 / np.cos(2 * t)
    th = 2 * np.arctan(a * zeta)
    G = math.log(1.0 / math.tan((math.pi - a0) / 4))            # Green function of C\A at z = 1
    return th, G
def logabs_poly(c, th):
    z = np.exp(1j * th)
    v = np.zeros_like(z)
    for ck in c[::-1]:
        v = v * z + ck
    return np.log(np.abs(v) + 1e-300)
L = LOG; c = H[L]; M = c.sum(); lx = L * math.log(2)
print()
print(f"# x = 2^{L}: M = {int(M)}, Theta_true = log|M|/log x = {math.log(abs(M))/lx:.4f}")
print("# check: P = (1+z)^L must give equality in the two-constants bound (harmonic measure); slack of the sup form = log2(1+sin pi d)")
print("delta   G/ln2   c(d)=-log2cos   slack  | (1+z)^L: lhs/L  2const/L  sup/L | mu: Theta_2const Theta_sup | M-part only | random model Theta_true Theta_2c")
rng = np.random.default_rng(7)
b = binomw(L); cnt = b * 2.0**L
crand = rng.normal(0, 1, len(c)) * np.sqrt(8 / math.pi**2 * oddtop(L, len(c)))
Mr = crand.sum()
for delta in (0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3):
    th, G = arc_nodes(delta)
    g = G / math.log(2); cd = -math.log2(math.cos(math.pi * delta)); sl = math.log2(1 + math.sin(math.pi * delta))
    lp = L * np.log(np.abs(1 + np.exp(1j * th)) + 1e-300)
    two_c = (L * G + lp.mean()) / L / math.log(2); sup_c = (L * G + lp.max()) / L / math.log(2)
    lf = logabs_poly(c, th)
    th2 = (L * G + lf.mean()) / lx; ths = (L * G + lf.max()) / lx
    lm = math.log(abs(M)) + L * np.log(np.abs((1 + np.exp(1j * th)) / 2) + 1e-300)   # closed form, no rounding
    thm = (L * G + lm.mean()) / lx
    lr = logabs_poly(crand, th); thr = (L * G + lr.mean()) / lx
    print("%.3f  %6.4f  %6.4f   %6.4f  |  %.4f  %.4f  %.4f |  %.4f  %.4f  |  %.4f  |  %.4f  %.4f" % (
        delta, g, cd, sl, 1.0, two_c, sup_c, th2, ths, thm, math.log(abs(Mr)) / lx, thr))
# needed extra saving to beat B_2 via the sup form, and the harmonic-measure criterion
B2 = (777 - math.sqrt(1081)) / 852
print()
print(f"# To beat B_2 = {B2:.10f}: sup form needs sigma(d) > (1-B_2) + c(d) + log2(1+sin pi d); two-constants form needs")
print("#   int_A (sigma - c) d omega_1 > 1 - B_2 = %.5f, where sigma(alpha) is the proven saving exponent and c(alpha) = -log2|cos pi alpha|." % (1 - B2))
print("delta  sup-form needed sigma(d)   frequency-0 value c(d)   proven-method cap (= c(d), Type I term m = 1)")
for delta in (0.01, 0.05, 0.1, 0.2, 0.3, 0.4):
    cd = -math.log2(math.cos(math.pi * delta)); sl = math.log2(1 + math.sin(math.pi * delta))
    print("%.2f   %.4f                     %.4f                   %.4f" % (delta, (1 - B2) + cd + sl, cd, cd))
