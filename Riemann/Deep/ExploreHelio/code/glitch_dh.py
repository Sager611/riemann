# glitch_dh.py (explore23, T2/T1 adversarial): glitch spectrum and Duvall/Turing count of Davenport-Heilbronn ordinates,
# and the moved-pair test for zeta. Ordinate-only statistics vs explicit-formula residuals.
import json, numpy as np
PI = np.pi
D = json.load(open('dh_zeros.json'))
on = np.array(D['online']); off = np.array(D['offline'])   # off: [beta, gamma] with beta > 1/2
k = (np.sqrt(10 - 2*np.sqrt(5)) - 2)/(np.sqrt(5) - 1)
amod = {0: 0.0, 1: 1.0, 2: k, 3: -k, 4: -1.0}
NB = 100000
a = np.array([amod[n % 5] for n in range(NB + 1)]); a[0] = 0
b = np.zeros(NB + 1); S = np.zeros(NB + 1)
lg = np.log(np.maximum(np.arange(NB + 1), 1))
for n in range(1, NB + 1):
    b[n] = a[n]*lg[n] - S[n]
    if b[n] != 0.0:
        m = np.arange(2, NB//n + 1)
        S[m*n] += a[m]*b[n]
out = {'b_small': [[n, round(float(b[n]), 6)] for n in range(1, 31)]}
# glitch amplitude growth: psi_DH(x) = sum_{n<=x} b(n) (no pole) vs zeta's psi(x) - x
psi = np.cumsum(b)
xs = np.unique(np.logspace(1, 5, 400).astype(int))
r = np.abs(psi[xs])/np.sqrt(xs)
out['psiDH_over_sqrtx_max_by_decade'] = [[int(10**e), round(float(np.max(r[(xs >= 10**e) & (xs < 10**(e + 1))])), 3)] for e in range(1, 5)]
# zeta: psi(x) - x
lam = np.zeros(NB + 1)
isp = np.ones(NB + 1, bool); isp[:2] = False
for p in range(2, int(NB**0.5) + 1):
    if isp[p]: isp[p*p::p] = False
for p in np.nonzero(isp)[0]:
    q = p
    while q <= NB: lam[q] = np.log(p); q *= p
psiz = np.cumsum(lam)
rz = np.abs(psiz[xs] - xs)/np.sqrt(xs)
out['psiZeta_minus_x_over_sqrtx_max_by_decade'] = [[int(10**e), round(float(np.max(rz[(xs >= 10**e) & (xs < 10**(e + 1))])), 3)] for e in range(1, 5)]
# fit growth exponent of max|psi_DH| on [10^3, 10^5]
blk = np.unique(np.logspace(3, 5, 9).astype(int))
mx = [np.max(np.abs(psi[blk[i]:blk[i + 1]])) for i in range(len(blk) - 1)]
mid = np.sqrt(blk[:-1]*blk[1:])
out['psiDH_growth_exponent_1e3_1e5'] = round(float(np.polyfit(np.log(mid), np.log(mx), 1)[0]), 3)
# --- ordinate lists
gfull_c = np.concatenate([on.astype(complex), off[:, 1] - 1j*(off[:, 0] - 0.5), off[:, 1] + 1j*(off[:, 0] - 0.5)])  # complex gamma_rho
gord = np.concatenate([on, off[:, 1], off[:, 1]])        # ordinates, off-line pair = double ordinate
gon = on.copy()                                            # sign changes only (Turing deficit)
T1, T2 = 20.0, 420.0
def w(z):
    x = (z - T1)/(T2 - T1)
    return np.where((x.real > 0) & (x.real < 1), np.sin(PI*x)**2, 0.0)
def G(gs, lams, norm):
    return np.array([(w(gs)*np.exp(1j*l*gs)).sum() for l in lams])/norm
Nw = w(gord).real.sum()
Lb = (w(gord).real*np.log(5*gord/(2*PI))).sum()/Nw
out['dh_window'] = [T1, T2, int(len(on)), int(len(off)), round(float(Nw), 2), round(float(Lb), 4)]
rows = []
for n in [2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14, 16, 17, 18, 19, 21, 22, 23, 24]:
    l = np.log(n)
    rows.append([n, round(float(-b[n]/np.sqrt(n)/Lb), 4), round(float(G(gfull_c, [l], Nw)[0].real), 4), round(float(G(gord, [l], Nw)[0].real), 4), round(float(G(gon, [l], Nw)[0].real), 4)])
out['dh_direct_peaks_[n,pred,complex,ordinates,online]'] = rows
# generic-lambda residual (away from log n for b(n) != 0 large): explicit formula says complex list ~ prediction; ordinates deviate by 2(cosh(lam delta)-1)
lr = np.linspace(0.35, 3.3, 600)
ln_ = np.log(np.arange(2, 200)); big = np.abs(b[2:200]/np.sqrt(np.arange(2, 200))) > 0.02
lr = lr[np.min(np.abs(lr[:, None] - ln_[big][None, :]), axis=1) > 0.03]
def pred(lams):
    # prime side: -(1/2pi) sum_n b(n) n^{-1/2} [W(lam - log n)] / Nw, W = FT of the window (numerical), conj term negligible
    tt = np.linspace(T1, T2, 40001); wt = np.sin(PI*(tt - T1)/(T2 - T1))**2; dt = tt[1] - tt[0]
    res = []
    for l in lams:
        s = 0j
        for n in range(2, 2000):
            if b[n] != 0 and abs(l - np.log(n)) < 0.5:
                s += -b[n]/np.sqrt(n)*(wt*np.exp(1j*(l - np.log(n))*tt)).sum()*dt/(2*PI)
        res.append(s)
    return np.array(res)/Nw
P = pred(lr)
for name, gs in [('complex', gfull_c), ('ordinates', gord), ('online_only', gon)]:
    g_ = G(gs, lr, Nw)
    out['dh_resid_rms_' + name] = float(np.sqrt(np.mean(np.abs(g_ - P)**2)))
# predicted ordinate residual from the off-line pairs alone: -(1/Nw) sum 2(cosh(lam d) - 1) e^{i lam g} w(g)
dd = off[:, 0] - 0.5
Rth = -np.array([(2*(np.cosh(l*dd) - 1)*np.exp(1j*l*off[:, 1])*w(off[:, 1])).sum() for l in lr])/Nw
g_ = G(gord, lr, Nw)
out['dh_ordinate_resid_vs_cosh_formula_rms'] = [float(np.sqrt(np.mean(np.abs(g_ - P)**2))), float(np.sqrt(np.mean(np.abs(Rth)**2))), float(np.sqrt(np.mean(np.abs(g_ - P - Rth)**2)))]
# Duvall/Turing diagram: pi(n - 1/2)/g_n - F(5 g_n) residual; F from the gamma factor (vartheta/g)
json.dump(out, open('glitch_dh.json', 'w'), indent=1)
print(json.dumps(out)[:5000])
