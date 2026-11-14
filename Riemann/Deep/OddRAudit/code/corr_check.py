#!/usr/bin/env python3
# corr_check.py -- oddraudit: Lemma O2 against the exact alpha_2 (floor formula, mpmath 60 digits) and the normalised odd
# correction (sqrt x log x)[2/x + A/(L_N log L_N)] of transfer.py at sample x; also RobinCH margins at the same boxes.
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mpmath import mp, mpf, log, floor, sqrt
mp.dps = 60
for xs in ['41', '1e3', '5e8', '1e9', '1e12', '1e20', '1e30', '1e40', '6.46e46', '1.43e47']:
    x = mpf(xs); eps = log(1 + 1/x)/log(x)
    a = int(floor(log((2**(1 + eps) - 1)/(2**eps - 1))/log(2))) - 1
    k = int(floor(log(x)/log(2)))
    A = log(((x + 1)*log(x)/log(2) + 2)/2)
    LN = x*(1 - 1/log(x)) - log(2)
    corr = sqrt(x)*log(x)*(2/x + A/(LN*log(LN)))
    print(f'x={xs:>8}: alpha_2 = {a}, floor(log2 x) = {k}, a >= k: {a >= k}, 2^-a < 2/x: {mpf(2)**-a < 2/x}, a log2 <= A: {a*log(2) <= A} (A/log2 = {mp.nstr(A/log(2), 6)}), normalised corr = {mp.nstr(corr, 4)}')
from flint import arb
import transfer as T
import certify3 as Q
e1 = Q.eta1_ch(); b0 = Q.parse_b0('1')
for Ls in ['20.03', '20.5', '30', '60', '100', '107.5']:
    La = arb(Ls); Lb = La + arb('0.01')
    m, arg = Q.margins4('R', La, Lb, b0, e1)
    print(f'box log x in [{Ls}, +0.01]: RobinCH m >= {float(m.lower()):.5g}, odd corr <= {float(arb(T.corr(La, Lb)).upper()):.4g}')
