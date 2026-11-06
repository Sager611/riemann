# Rigorous wrap-up of rhow.c (k = 6): turn the computed max ratio into a proved upper bound for ||B||.
# Per level-6 cell the numerator is a recursive double sum of <= T positive products fl(WUP*h) (round-to-nearest),
# WUP >= w exactly (weights.py), h = h_5 the stored test function; then one division.
# Higham, Accuracy and Stability, Lemma 3.1/eq. (4.4): s_true <= s_hat/((1-u)(1-gamma_{T-1})), ratio <= fl(s/h)/(1-u).
from fractions import Fraction as F
import re, sys
log = open(sys.argv[1]).read()
m = re.search(r'RESULT maxr_hex (\S+) minh_hex (\S+) maxterms (\d+) flags (\d+)', log)
maxr = F(float.fromhex(m.group(1))); minh = float.fromhex(m.group(2)); T = int(m.group(3)); flags = int(m.group(4))
assert flags == 0 and minh > 0
u = F(1, 2**53)
gam = (T-1)*u/(1-(T-1)*u)
bound = maxr/((1-u)**2*(1-gam))
used = F(11101037276778, 10**12)
print('computed max ratio  %.17g' % float(maxr))
print('proved bound ||B|| <= %.17g  (inflation %.3e, T = %d terms)' % (float(bound), float(bound/maxr-1), T))
print('used in the certificate rho_W^up = %.15g ; slack %.3e' % (float(used), float(used-bound)))
print('min h_5 = %.17g > 0' % minh)
print('VERDICT', 'rho_W^up VALID' if bound < used else 'rho_W^up NOT PROVED')
