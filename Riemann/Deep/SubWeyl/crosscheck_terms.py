# Float cross-check: the per-block terms of certify.py (F1, F2) against the terms of Patel-Yang Lemma 3.4 as stated
# (before its square roots are split).  For random (a, t, theta, eta, h) with q = theta a^{36/41} t^{-11/41} treated
# as exact, each certify term times a t^{27/82} must equal the corresponding Lemma 3.4 term (squared constant / alpha^2).
# uv run --no-project --with python-flint --with numpy python3 crosscheck_terms.py
import random, math
from certify import *
random.seed(1)
worst = 0.0
for trial in range(200):
    L = random.uniform(80, 2000); lt = L
    beta = random.uniform(0.42, 0.5); la = beta*L
    th = random.choice(THETAS); eta = random.choice(ETA1)
    H = 1 + random.uniform(0.005, 0.06)
    rc = floats_of(region_consts(arb(H)))
    d = rc['aba'][eta]
    E6f = float(E6.mid())
    # certify terms (q = Q exactly, g = 1), log space
    T = {
     'i+ii': math.log((H - 1)/th + d['c2']*th**(11/30)) + 5/41*la - 5/82*lt,
     'iii': math.log(d['c3']*th**(61/120)) - 17/164*la + 7/164*lt,
     'iv': math.log(rc['c4']/math.sqrt(th)) + 87/82*la - 57/82*lt,
     'v': math.log(th**(1/3)*(rc['c5a']*math.exp(-la/3) + rc['c5b'])) - 5/123*la - 7/82*lt,
     'vi': math.log(E6f) - 27/82*lt}
    # Patel-Yang Lemma 3.4 squared-sum terms divided by alpha^2 (TeX l. 844-858), with theta_2 := a t^{-7/17} (so the
    # a^{-1/3} in E_5 is exact), q_0 -> infinity, A_5 at h_5 = HC5 h^9 kappa, and the same kappa^{1/8}:
    A5, B5 = [float(x.mid()) for x in AkBk(5, arb(eta), HC5*arb(H)**9*((arb(2336)/(2*arb(H)**2*(arb(H)-1)))/((arb(2336)/(2*arb(H)**2*(arb(H)-1))) - 1))**fr(9, 2))]
    kap = float((((arb(2336)/(2*arb(H)**2*(arb(H)-1)))/((arb(2336)/(2*arb(H)**2*(arb(H)-1))) - 1))**fr(9, 2)).mid())
    K3f, K4f = float(K3.mid()), float(K4.mid())
    P = {
     'i+ii': math.log((H - 1)/th + K3f*th**(11/30)*A5*kap**0.125*H**(21/8)*(H - 1)) + 46/41*la + 11/41*lt,
     'iii': math.log(K4f*th**(61/120)*B5*H**1.5*(H - 1)**(7/8)) + 147/164*la + 61/164*lt,
     'iv': math.log(8/3*4.686*math.sqrt(math.pi)*H**1.5*th**(-0.5)) + 169/82*la - 15/41*lt,
     'v': math.log(9/14*th**(1/3)*(0.893*5/math.pi*((H - 1)/math.pi)**(1/3)*math.exp(-la/3) + 6/math.pi*H**3*(H - 1))) + 118/123*la + 10/41*lt,
     'vi': math.log(E6f) + la}
    for k in T:
        diff = abs((T[k] + la + 27/82*lt) - P[k])
        worst = max(worst, diff)
print('max |log ratio| certify-term vs Lemma 3.4 term over 200 random cases: %.3e' % worst)
assert worst < 1e-9
print('TERMS OK')
