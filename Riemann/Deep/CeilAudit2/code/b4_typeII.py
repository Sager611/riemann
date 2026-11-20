# B4 test: BeyondB2's barrier "the one provable signed tool, a large sieve in the Mobius variable (LS1), loses Z^0.0333 at the contact".
# Freedom: a different inequality in the same step: split mu (Vaughan / Heath-Brown) into bilinear Type II pieces sum alpha_a beta_b chi_u(ab)
# and bound sum_{u in B} b_u A_u B_u by Cauchy-Schwarz and a family mean value with ANY large-sieve constant Delta(N) >= max(U, N) (even the
# conjectural diagonal sextic one). Exponents in Z (U = Z^h, #B = U^R, |F_u| <= Z^a, |b_u| <= 1):
#   absolute: a + hR;  LS1: 1 + hR/2;  Type II: >= 1/2 + min_alpha [max(h, alpha) + max(h, 1 - alpha)]/2 = 1/2 + h (as h > 1/2).
# Subset large sieve (Delta = N + #B U^(1/2)): 1/2 + [max(alpha, hR + h/2) + max(1 - alpha, hR + h/2)]/2.
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'B2Below'))
import numpy as np
from system import Rd_native, D2, R2
h2 = 0.8096863494; d2 = float(D2); R2 = float(R2); a2 = (1 + d2) / 2
print("contact: delta=%.6f a=%.6f h=%.6f R=%.6f" % (d2, a2, h2, R2))
absr = a2 + h2 * R2; ls1 = 1 + h2 * R2 / 2; t2 = 0.5 + h2
al = np.linspace(0, 1, 100001); sub = 0.5 + (np.maximum(al, h2 * R2 + h2 / 2) + np.maximum(1 - al, h2 * R2 + h2 / 2)) / 2
print("absolute %.6f  LS1 %.6f (loss %.7f)  TypeII %.6f (loss %.6f)  subset-LS TypeII %.6f (loss %.6f)" % (
    absr, ls1, ls1 - absr, t2, t2 - absr, sub.min(), sub.min() - absr))
d = np.linspace(1 / 8, 5 / 6, 1401); x = np.linspace(0, 0.5, 201); dd, xx = np.meshgrid(d, x)
R = Rd_native(dd, xx); a = (1 + dd) / 2
loss2 = 0.5 + h2 - a - h2 * R; loss1 = 1 - a - h2 * R / 2
i = np.unravel_index(np.argmin(loss2), loss2.shape)
print("class domain [1/8,5/6]x[0,1/2], native count: min LS1 loss %.5f; min TypeII loss %.5f at (delta, x) = (%.3f, %.3f), R = %.4f" % (
    loss1.min(), loss2.min(), dd[i], xx[i], R[i]))
print("TypeII beats absolute on %.2f%% of the grid (classes with h(1-R) < delta/2)" % (100 * (loss2 < 0).mean()))
