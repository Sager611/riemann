#!/usr/bin/env python3
"""lcaudit falsification sanity check, block (1000, 2000] only (copy of zeros_F.py with these parameters).
 (setting only; the theorem is a liminf statement).
q = 5: chi_a = (5,2) and chi_b = (5,3) (complex conjugate quartic, both odd): F = b1 e^{i th_a} L_a + b2 e^{i th_b} L_b, b = (1, +-1).
Counts sign changes of the real function Z_F(t) = e^{i vartheta(t)} F(1/2+it) in dyadic blocks against N_F's main term;
L_a alone as calibration; opposite parity chi_a (odd) + chi_q = (5,4) (even): no real normalisation exists.
Run: uv run --no-project --with python-flint python3 zeros_F.py
"""
import math, cmath, time
from flint import acb, arb, dirichlet_char, ctx
ctx.prec = 64
q = 5
def gauss(chi):
    return sum(complex(chi(n)) * cmath.exp(2j*math.pi*n/q) for n in range(1, q + 1))
def phase_theta(chi):          # e^{i theta_chi} = eps^{-1/2}, eps = tau(chi)/(i^a sqrt q)
    a = chi.parity()
    eps = gauss(chi)/((1j)**a*math.sqrt(q))
    return cmath.exp(-0.5j*cmath.phase(eps)), eps, a
def vartheta(t, a):            # phase of (q/pi)^{(s+a)/2} Gamma((s+a)/2), s = 1/2 + it
    lg = acb((0.5 + a)/2, t/2).lgamma()
    return float(lg.imag.mid()) + (t/2)*math.log(q/math.pi)
def L(chi, s):
    v = acb.dirichlet_l(acb(s.real, s.imag), chi)
    return complex(float(v.real.mid()), float(v.imag.mid()))
out = []
def say(m): print(m, flush=True); out.append(m)
ca, cb, cq = dirichlet_char(5, 2), dirichlet_char(5, 3), dirichlet_char(5, 4)
ea, epsa, pa = phase_theta(ca); eb, epsb, pb = phase_theta(cb); eq, epsq, pq = phase_theta(cq)
say("chi_a(2) = %s, parity %d, eps = %s; chi_b(2) = %s, parity %d, eps = %s; chi_q parity %d" % (ca(2), pa, epsa, cb(2), pb, epsb, pq))
say("a_F(1) for b = (1, 1): %s ; for b = (1, -1): %s" % (ea + eb, ea - eb))
# Davenport-Heilbronn check: known off-line zero 0.808517 + 85.699348 i
for sgn in (1, -1):
    s0 = complex(0.808517, 85.699348)
    Fv = ea*L(ca, s0) + sgn*eb*L(cb, s0)
    say("b2 = %+d: |F(0.808517 + 85.699348 i)| = %.2e  (|L_a| there = %.3f)" % (sgn, abs(Fv), abs(L(ca, s0))))
def Nmain(T): return (T/(2*math.pi))*math.log(q*T/(2*math.pi*math.e))
h = 0.04
blocks = [(1000.0, 2000.0)]
t = 1000.0
prev = None; maxim = 0.0
cnt = {k: [0, 0, 0] for k in range(len(blocks))}     # sign changes: L_a, F(b2=+1), F(b2=-1)
minopp = (1e9, None); dphi = []
t0 = time.time()
while t <= 2000.0:
    s = complex(0.5, t)
    va, vb = L(ca, s), L(cb, s)
    ph = cmath.exp(1j*vartheta(t, pa))
    Za, Zb = ph*ea*va, ph*eb*vb
    maxim = max(maxim, abs(Za.imag)/(abs(Za) + 1e-300), abs(Zb.imag)/(abs(Zb) + 1e-300))
    cur = (Za.real, Za.real + Zb.real, Za.real - Zb.real)
    if prev is not None:
        for k, (lo, hi) in enumerate(blocks):
            if lo < t <= hi:
                for j in range(3):
                    if prev[j]*cur[j] < 0: cnt[k][j] += 1
    prev = cur
    if t <= 500.0 and abs(round(t/h) % 5) == 0:          # opposite parity on a coarser grid
        vq = L(cq, s)
        Fo = ea*va + eq*vq
        r = abs(Fo)/(abs(va) + abs(vq))
        if r < minopp[0]: minopp = (r, t)
        if t > 100 and len(dphi) < 5 and abs(t - round(t)) < 1e-9 and round(t) % 100 == 0:
            dphi.append((t, (vartheta(t, pa) - vartheta(t, pq)) % math.pi))
    t = round(t + h, 10)
say("max |Im Z_j|/|Z_j| on the grid (reality of e^{i vartheta} e^{i theta_j} L_j): %.1e" % maxim)
say("block       N_main  #sgn L_a   frac  #sgn F(1,1)  frac  #sgn F(1,-1)  frac")
for k, (lo, hi) in enumerate(blocks):
    N = Nmain(hi) - Nmain(lo)
    a_, b_, c_ = cnt[k]
    say("(%5.1f,%6.1f] %7.1f  %5d  %6.3f  %8d  %6.3f  %9d  %6.3f" % (lo, hi, N, a_, a_/N, b_, b_/N, c_, c_/N))
if minopp[1] is not None: say("opposite parity (chi_a odd + chi_q even, b = (1,1)), t in [10,500]: min |F|/(|L_a|+|L_q|) = %.3e at t = %.2f" % minopp)
if dphi: say("  vartheta_odd - vartheta_even mod pi at t = 200..500: %s (pi/4 = %.4f)" % (", ".join("%.4f" % d for _, d in dphi), math.pi/4))
say("elapsed %.0f s, step h = %.2f" % (time.time() - t0, h))
open('logs/zeros_F2.log', 'w').write('\n'.join(out) + '\n')
