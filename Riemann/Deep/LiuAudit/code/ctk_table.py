#!/usr/bin/env python3
"""Local correlation table C_{t,k}(J) of Liu 2610.26835v2 §13.3 (number_fields_v4.tex l.8578-8590), brute force over
F = Q(zeta12) at the split prime pi of norm q = 13 (and 37 for t+k <= 3).
Reconstruction (our reading of the outer-variable Poisson in I_eta, l.8031-8035, with modulus s*c*n^3):
  F(t,k,h) = sum_{x mod pi^(t+k)} g_{pi^k}(-x) chi_pi(x)^t e(h x / pi^(t+k)),  g_s(y) = sum_{d mod s} chi_s(d) e(y d/s),
and the claim tested is F(t,k,h) = q^k * C_{t,k}(J) * rho^(k-t) with J = v_pi(h), rho = chi_pi(h/pi^J) (or its conjugate),
C_{0,0}=1, C_{0,k}=1_{J=0}, C_{t,0}=1_{J>=t-1} q^{t-1} T_t(p^{J-t+1}), C_{t,1}=1_{J>=t} tau_1 q^{t-1} T_{t-1}(p^{J-t}), C_{t,k>=2}=0,
T_v(1) = tau_v (tau_0 = -1), T_v(p^{>0}) = (q-1) 1_{6|v}; tau_v = sum_{x mod pi} chi^v(x) e(x/pi)."""
import sys, os, math, cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from gauss_q12 import Prime, roots_mod_p, mul, pw, tr, ev, ONE

def run(P, tk_max, big=()):
    p = P.p
    tau = [complex(np.sum(P.tab[np.arange(1, p)]**v * np.exp(2j*np.pi*((np.arange(1, p)*tr(P.pistar)) % p)/p))) for v in range(6)]
    tau[0] = -1+0j
    out = []
    for t in range(0, 8):
        for k in range(0, 3):
            if t+k > tk_max and (t, k) not in big: continue
            Nm = p**(t+k); Nk = p**k
            # g_{pi^k}(y) for y mod p^k
            if k == 0:
                gk = np.ones(1, dtype=complex)
            else:
                d = np.arange(Nk, dtype=np.int64)
                chik = P.tab[d % p]**k
                Tk = tr(pw(P.pistar, k)) % Nk
                Gf = np.fft.ifft(chik)*Nk          # Gf[t'] = sum_d chi^k(d) exp(2 pi i d t'/Nk)
                y = np.arange(Nk, dtype=np.int64)
                gk = Gf[(y*Tk) % Nk]
            Tm = pw(P.pistar, t+k)
            for J in range(0, t+k+2):
                for h0 in ((1, 0, 0, 0), (2, 1, 0, 0), (0, 1, 1, 0)):
                    if ev(h0, P.r, p) == 0: continue
                    h = mul(pw(P.pi, J), h0)
                    th = tr(mul(h, Tm)) % Nm
                    tot = 0j; CH = 1 << 20
                    for s0 in range(0, Nm, CH):
                        x = np.arange(s0, min(Nm, s0+CH), dtype=np.int64)
                        gx = gk[(-x) % Nk] if k else 1.0
                        tot += np.sum(gx * (P.tab[x % p]**t if t else (np.ones(len(x)))) * np.exp(2j*np.pi*((x*th) % Nm)/Nm))
                    rho = P.chi(h0)
                    def T(v, e):
                        if e == 0: return tau[v % 6] if v % 6 else -1+0j
                        return (p-1) if v % 6 == 0 else 0
                    if t == 0 and k == 0: C = 1
                    elif t == 0: C = 1 if J == 0 else 0
                    elif k == 0: C = p**(t-1)*T(t, J-t+1) if J >= t-1 else 0
                    elif k == 1: C = tau[1]*p**(t-1)*T(t-1, J-t) if J >= t else 0
                    else: C = 0
                    pred = p**k*C*rho**(k-t); predc = p**k*C*np.conj(rho)**(k-t)
                    sc = max(1, abs(pred))
                    ok = abs(tot-pred) < 1e-6*sc; okc = abs(tot-predc) < 1e-6*sc
                    out.append((t, k, J, h0, ok, okc, tot, pred))
    return out

def main():
    allres = []
    for P, tkm, big in ((Prime(13, roots_mod_p(13)[0]), 6, ((6, 1), (7, 0))), (Prime(37, roots_mod_p(37)[0]), 3, ())):
        res = run(P, tkm, big)
        for t, k, J, h0, ok, okc, tot, pred in res:
            print(('PASS' if ok else ('PASS(conj)' if okc else 'FAIL')), f'q={P.p} t={t} k={k} J={J} h0={h0}: F={tot:.6g} pred={pred:.6g}', flush=True)
        allres += res
    n = len(allres); a = sum(r[4] for r in allres); b = sum(r[5] for r in allres); c = sum(r[4] or r[5] for r in allres)
    print('SUMMARY_CTK', c, 'pass of', n, f'(rho^(k-t) convention {a}, conjugate convention {b})')

if __name__ == '__main__':
    main()

