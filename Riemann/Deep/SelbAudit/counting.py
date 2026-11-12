#!/usr/bin/env python3
"""selbaudit check 4: (a) counting lemma m >= meas(G)/H vs Dousselin's pairing on toy random signals;
(b) independent mpmath check of the special-function constants of cert.log (float, not certified).
uv run --no-project --with numpy --with mpmath python3 counting.py"""
import math, numpy as np, mpmath as mp
rng = np.random.default_rng(1)
print("== (a) counting: m (sign changes in (T,2T+H)) vs meas(G)/H vs pairing floor(T/2H) - n1 ==")
for trial in range(6):
    Tlen = 200.0; H = [0.3, 0.7, 1.5, 3.0, 0.5, 1.0][trial]; h = 0.002
    t = np.arange(0, Tlen + H + h, h)
    K = 60; fr = rng.uniform(0.2, 4.0, K); ph = rng.uniform(0, 2*np.pi, K); am = rng.standard_normal(K)/np.sqrt(fr)
    sig = sum(a*np.cos(f*t + p) for a, f, p in zip(am, fr, ph))
    if trial >= 4: sig = sig + 2.5*np.std(sig)*np.sign(np.sin(0.05*t))   # long same-sign stretches (many no-change windows)
    s = np.sign(sig); ch = np.nonzero(s[1:] != s[:-1])[0]   # sign change between t[i], t[i+1]
    m = len(ch); w = int(round(H/h)); nT = int(round(Tlen/h))
    c = np.concatenate([[0], np.cumsum(s[1:] != s[:-1])])   # c[i] = #changes in (t0, t_i]
    G = np.array([c[min(i + w, len(c) - 1)] - c[i] > 0 for i in range(nT)])
    measG = G.sum()*h
    nI = int(Tlen//(2*H)); n1 = 0; pair_z = 0
    for k in range(nI):
        a = int(round(2*k*H/h)); b = a + w; e = a + 2*w
        if c[e] - c[a] > 0: pair_z += 1
        if not G[a:b].any(): n1 += 1
    print("H=%.2f: m=%d  meas(G)/H=%.1f  [m >= meas(G)/H: %s]  pairing count floor(T/2H)-n1=%d (pairs with a change: %d)  T/H=%.0f meas(E)/H=%.1f" % (H, m, measG/H, m >= measG/H - 1e-9, nI - n1, pair_z, Tlen/H, (Tlen - measG)/H))
print("== (b) special functions (mpmath, 30 digits; cert.log: R0max 2.338139 DA 4.348005 Emax 4.334941 Psimax 2.382235 Ipsi[0,1] 2.119367) ==")
mp.mp.dps = 30
psi = lambda Y: mp.sqrt(mp.pi)*mp.exp(-1j*mp.pi/4)*mp.erf(mp.exp(1j*mp.pi/4)*mp.sqrt(Y))
D2 = lambda Y: -2*mp.sqrt(Y)*mp.exp(-1j*Y) - (1 + 2j*Y)*psi(Y)
E = lambda Y: D2(Y)/mp.sqrt(Y) if Y > 0 else mp.mpf(-4)
c1 = 2*mp.sqrt(mp.pi)*mp.exp(-3j*mp.pi/4)
Ys = [i/400 for i in range(1, 401)]
print("Emax ~ %.6f (sup |E| on [0,1]; |E(0)| = 4)" % max([4.0] + [float(abs(E(y))) for y in Ys]))
Ys2 = [i/50 for i in range(1, 2500)]
print("R0max ~ %.6f (sup |Delta_2(Y) - c1 Y| on (0,50])" % max(float(abs(D2(y) - c1*y)) for y in Ys2))
print("DA ~ %.6f (sup |Delta_2(Y)|/max(1,Y) on (0,50])" % max(float(abs(D2(y))/max(1, y)) for y in Ys2))
print("Psimax ~ %.6f (sup |psi| on (0,50])" % max(float(abs(psi(y))) for y in Ys2))
W = lambda l: l*mp.quad(lambda Y: abs(psi(Y))/Y**2, [l, 1, 5, 20, 100, mp.inf])
print("Ipsi[0,1] ~ %.6f (sup_l l*int_l^inf |psi|/Y^2, l in (0,1]; sampled)" % max(float(W(l)) for l in [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.85, 1.0]))
th = 0.15
i0 = 2*th + 2*(1 - th)/3; i1 = 12*th*th + 8*th + 2*(1 - th)*(16/3 - 3*(1 - th))
print("Phi_BB(t) = min(1,(1-t/2)/(1-theta)) exactly: i0 = %.5f, i1 = %.5f at theta=0.15 (cert step-function upper values 0.87381, 6.29385)" % (i0, i1))
