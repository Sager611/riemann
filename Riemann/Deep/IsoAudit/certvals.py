# certvals.py: exact certified lower endpoints at the stated thresholds, via the unmodified IsoZero cert.py (imported read-only).
# Run: cd Riemann/Deep/IsoAudit && uv run --no-project --with python-flint --with mpmath python3 certvals.py
import sys, time
sys.path.insert(0, '../IsoZero')
import cert
from flint import arb
t0 = time.time()
for c, tgt in (('0.00476', 0.6724), ('0.01106', 0.672), ('0.02485', 0.67), ('0.05495', 0.66), ('0.06427', 0.6552), ('0.001', 0.672487)):
    p, rho, ar, Bv, lam, r0 = cert.hinge_p(c, back=(0.02,))
    print('c = %s: certified p >= %.10f (target %s, margin %.2e); rho %.4f rho_max %.4f alpha <= %.9f B <= %.6e lam %.6e  [%.0f s]' % (c, p, tgt, p - tgt, rho, r0, ar, Bv, lam, time.time() - t0)); sys.stdout.flush()
w = arb('3.75').sqrt()
p, rho, ar, Bv, lam, r0 = cert.hinge_p('0.1666666666666666666666666667', w, back=(0.02,))
print('c = 1/6, w^2 = 3.75: certified p >= %.10f; rho %.4f alpha <= %.6e B <= %.3e lam %.6f  [%.0f s]' % (p, rho, ar, Bv, lam, time.time() - t0))
