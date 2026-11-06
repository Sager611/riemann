# glitch.py (explore23, T2): acoustic-glitch analysis of zeta zero ordinates.
# Direct glitches (Landau): G(lam) = sum_n w_n exp(i lam g_n) / sum_n w_n  ~ -Lambda(e^lam) e^{-lam/2} / Lbar at lam = log n.
# Dual glitches (Gonek's lemma, functional equation): C_m = sum_n w_n exp(2i theta(g_n)) m^{-i g_n} / sum w_n ~ -m^{-1/2} / Lbar for EVERY m.
import sys, json, numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from zeros import theta
PI = np.pi
g = np.load(sys.argv[1] if len(sys.argv) > 1 else 'zeta_zeros.npy')
out = {}
def vM(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return np.log(p) if m == 1 else 0.0
    return 0.0
def win(gs, T1, T2):
    x = (gs - T1)/(T2 - T1)
    return np.where((x > 0) & (x < 1), np.sin(PI*x)**2, 0.0)
def spec(gs, w, lams):
    res = np.empty(len(lams), dtype=complex)
    for i in range(0, len(lams), 64):
        L = lams[i:i+64]
        res[i:i+64] = (w[None, :]*np.exp(1j*L[:, None]*gs[None, :])).sum(1)
    return res/w.sum()
T1, T2 = 2000.0, 20000.0
w = win(g, T1, T2)
Lbar = (w*np.log(g/(2*PI))).sum()/w.sum()
out['window'] = [T1, T2, int((w > 0).sum()), float(Lbar)]
# direct peaks at log n, n = 2..40; also refine the max near each log n (peak width ~ 2pi/(T2-T1))
rows = []
for n in range(2, 41):
    lam0 = np.log(n)
    loc = lam0 + np.linspace(-2e-3, 2e-3, 81)
    s = spec(g, w, loc)
    k = np.argmax(np.abs(s))
    pred = -vM(n)/np.sqrt(n)/Lbar
    rows.append([n, round(float(s[40].real), 5), round(float(s[40].imag), 5), round(float(pred), 5), round(float(loc[k] - lam0), 6), round(float(abs(s[k])), 5)])
out['direct'] = rows
# background (generic lam away from all log n): rms
rng = np.random.default_rng(1)
lr = rng.uniform(0.3, 3.7, 400)
lr = lr[np.min(np.abs(lr[:, None] - np.log(np.arange(2, 60))[None, :]), axis=1) > 0.01]
bg = spec(g, w, lr)
out['background_rms'] = float(np.sqrt(np.mean(np.abs(bg)**2)))
# dual peaks: C_m, m = 1..16
th = theta(g)
dual = []
for m in range(1, 17):
    c = (w*np.exp(1j*(2*th - g*np.log(m)))).sum()/w.sum()
    dual.append([m, round(float(c.real), 5), round(float(c.imag), 5), round(float(-1/np.sqrt(m)/Lbar), 5), round(float(-vM(m)/np.sqrt(m)/Lbar) if m > 1 else 0.0, 5)])
out['dual'] = dual
# dual peak seen in the plain ordinate spectrum: narrow window so the moving frequency log(T/2pi m) is sharp
dw = []
for (a, b) in [(9000.0, 11000.0), (17000.0, 19000.0)]:
    ww = win(g, a, b)
    Tm = 0.5*(a + b)
    for m in [1, 2, 3, 4, 6, 10]:
        lam = np.log(Tm/(2*PI*m))
        loc = lam + np.linspace(-0.08, 0.08, 321)
        s = spec(g, ww, loc)
        k = np.argmax(np.abs(s))
        dw.append([a, b, m, round(float(lam), 4), round(float(loc[k]), 4), round(float(abs(s[k])), 5)])
    bgs = spec(g, ww, np.linspace(1.0, 4.0, 300))
    dw.append([a, b, 'median_abs_bg', round(float(np.median(np.abs(bgs))), 5)])
out['dual_in_plain_spectrum'] = dw
json.dump(out, open('glitch_zeta.json', 'w'), indent=1)
print(json.dumps(out)[:6000])
