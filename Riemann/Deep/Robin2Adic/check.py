#!/usr/bin/env python3
# check.py -- independent float cross-check of cert.py (mpmath; no Arb Gamma, no sigma-cells, no running min).  ~1-2 min, < 100 MB
#   uv run --no-project --with python-flint --with mpmath python3 Riemann/Deep/Robin2Adic/check.py
# (1) W(s) = 4 int_{H0}^oo Nbar(s,t) t^-3 dt by quadrature in u = log t (exact log kt), against cert.W_kln / W_cs (Arb closed form);
# (2) zero sum as the continuous Stieltjes integral x^{-1/2} tau + int_{1/2}^{B} (L/2)(x^{s-1} - x^{-s}) 2 min(tau/2, W(s)) ds,
#     W from mpmath gammainc (closed form, independent code), composite Simpson, step 2^-16 on [0.8, B];
# (3) full majorant at L1 = 108.57 against cert.py (cert must be >= check, ratio close to 1).
import os, sys, json, importlib.util, mpmath as mp
mp.mp.dps = 20
sys.dont_write_bytecode = True                             # no __pycache__ next to the imported files
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('r2a_cert', os.path.join(HERE, 'cert.py'))
C = importlib.util.module_from_spec(spec); spec.loader.exec_module(C)          # library import (main guarded)
from flint import arb
def P(*a): print(*a, flush=True)
H0 = mp.mpf(3e12); h = mp.log(H0); B = mp.mpf('0.87338196661320752'); tau = mp.mpf(C.TAU.str(20, radius=False))
I2 = (h*h/2 + h/2 + mp.mpf(1)/4)/H0**2; I1 = (h/2 + mp.mpf(1)/4)/H0**2
def af(x): return mp.mpf(x.str(22, radius=False))
KP = C.I4.build_kln_H0('3e12')
KR = [(af(p.c0), af(p.c1), af(p.q0), af(p.q1), af(p.k), af(p.smin)) for p in KP]
CSR = [C.I4.CSRow(row) for row in C.I4.CS_ROWS]
CS = [(af(r.s1), af(r.s2), mp.mpf(row[2]), af(r.B1), mp.mpf(row[3]), mp.mpf(row[4])) for r, row in zip(CSR, C.I4.CS_ROWS)]
def Wk_quad(s):
    best = mp.inf
    for c0, c1, q0, q1, k, smin in KR:
        if s > smin:
            A = mp.exp(c0 + c1*s); Bc = max(q0 - q1*s, 0); a = 8*(1 - s)/3
            g = lambda u: (A*(u + mp.log(k))**(2*s)*u**(5 - 4*s)*mp.e**(a*u) + Bc*u**2)*mp.e**(-2*u)
            best = min(best, 4*mp.quad(g, [h, h + 5, h + 20, h + 80, mp.inf]))
    return best
def Wc_quad(s):
    best = mp.inf
    for s1, s2, b1r, b1c, b2, b3 in CS:
        if s1 <= s <= s2:
            a = 3*(1 - s)/(2 - s); e = (7 - 5*s)/(2 - s)
            g = lambda u: (b1c*mp.e**(a*u)*u**e + b2*u**2 + b3*u)*mp.e**(-2*u)
            best = min(best, 4*mp.quad(g, [h, h + 5, h + 20, h + 80, mp.inf]))
    return best
def Ig(q, e): return mp.gammainc(e + 1, q*h)/q**(e + 1)
def Wk_cf(s):
    best = mp.inf
    for c0, c1, q0, q1, k, smin in KR:
        if s > smin:
            best = min(best, 4*(mp.exp(c0 + c1*s)*Ig(2 - 8*(1 - s)/3, 5 - 2*s) + max(q0 - q1*s, 0)*I2))
    return best
def Wc_cf(s):
    best = mp.inf
    for s1, s2, b1r, b1c, b2, b3 in CS:
        if s1 <= s <= s2:
            a = 3*(1 - s)/(2 - s); e = (7 - 5*s)/(2 - s)
            best = min(best, 4*(b1c*Ig(2 - a, e) + b2*I2 + b3*I1))
    return best
P('(1) W(s): Arb closed form (cert) / quadrature; KLN quad uses the exact (log kt)^{2s}, so ratio >= 1:')
C.KLN_P = KP; C.CS_R = [(r.s1, r.s2, arb(row[2]), r.B1, arb(row[3]), arb(row[4])) for r, row in zip(CSR, C.I4.CS_ROWS)]
for s in ['0.7', '0.8', '0.85', '0.86', '0.87', '0.873']:
    j = int(float(s)*2**16); sa = arb(j)/2**16; sm = mp.mpf(j)/2**16
    wk, wc = Wk_quad(sm), Wc_quad(sm)
    P('  s=%s  KLN %s ratio %s | CS %s ratio %s | mpmath-cf/quad KLN %s CS %s' % (s, mp.nstr(wk, 6), mp.nstr(af(C.W_kln(sa))/wk, 9),
      mp.nstr(wc, 6), mp.nstr(af(C.W_cs(sa))/wc, 9), mp.nstr(Wk_cf(sm)/wk, 9), mp.nstr(Wc_cf(sm)/wc, 9)))
def simpson(F, a, b, n):
    n += n % 2; hh = (b - a)/n
    return hh/3*(F(a) + F(b) + 4*sum(F(a + (2*i - 1)*hh) for i in range(1, n//2 + 1)) + 2*sum(F(a + 2*i*hh) for i in range(1, n//2)))
def zero_sum(L, Wf):
    x = mp.e**L
    F = lambda s: (L/2)*(x**(s - 1) - x**(-s))*2*min(tau/2, Wf(s))
    edges = [mp.mpf('0.52')] + [mp.mpf(j)/32 for j in range(17, 28)] + [B]      # CS rows change at multiples of 1/32
    tot = x**(-mp.mpf(1)/2)*tau + (x**(mp.mpf('0.52') - 1) + x**(-mp.mpf('0.52')))/2*tau   # zeros with s(rho) <= 0.52: g <= g_{0.52}
    for a, b in zip(edges[:-1], edges[1:]):
        n = int((b - a)*(2**16 if a >= mp.mpf('0.78') else 2**10)) + 2
        tot += simpson(F, a + mp.mpf('1e-15'), b - mp.mpf('1e-15'), n)
    return tot
L = mp.mpf('108.57'); x = mp.e**L; f = 1/L + 1/L**2
cr = json.load(open(os.path.join(HERE, 'logs', 'cert_result.json')))['res']['variants']
for name, Wf in [('CS_b4_corrected', Wc_cf), ('KLN_RH3e12', Wk_cf)]:
    Z = zero_sum(L, Wf)
    cpsi = 1 + mp.mpf('1.1199e-9') + mp.mpf('1.01624')*mp.e**(-L/6)*L/mp.log(2)
    low = f*mp.e**(-L/2)*(2*cpsi + (2 + mp.euler - mp.log(4*mp.pi))*(1 + 2*(L + 2)/(mp.mpf('14.13')*L*(L + 1))))
    d = mp.mpf('2.4530e-12'); Gt = d**2/2*(L + 1 + d)/((1 - d)**2*(L - 2*d)**2)
    tot = f*(1 + 2*(L + 2)/(L*(L + 1)*H0))*Z + low + mp.log(2*mp.pi)/(x*L) + Gt
    P('(2,3) %s: continuous zero term %s; majorant %s = 2^-%s; cert %s; cert/check %s' % (
        name, mp.nstr(f*Z, 7), mp.nstr(tot, 8), mp.nstr(-mp.log(tot, 2), 6), cr[name]['E'], mp.nstr(mp.mpf(cr[name]['E'])/tot, 7)))
