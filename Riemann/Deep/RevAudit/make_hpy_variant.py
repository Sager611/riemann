# Builds /tmp/revaudit_hpy/certify_hpy.py: Explicit611's certify.py with Revers's Lemmas 2.7/2.9 (constant 1 - 6/pi, beta with
# (6/pi-1)/(c W^{2/3}) + (1-6/pi)/W^{1/3}) replaced by Hiary-Patel-Yang's refereed Lemmas 2.5/2.6 (arXiv 2207.02366: constant 2 - 4/pi,
# beta = 32/(3 sqrt(pi eta)) + (2-4/pi)/W^{1/3}); everything else unchanged. Every replacement must match exactly once.
import os, sys
src = open('/Users/adriansager/playground/riemann/Riemann/Deep/Explicit611/certify.py').read()
R = [
 # blk_arb (relaxed per-block bound)
 ("    be = 32 / (3 * (PI * cc).sqrt()) + (6 / PI - 1) / (cc * W3 ** 2) + (1 - 6 / PI) / W3\n    Q = (L + cc * W3)",
  "    be = 32 / (3 * (PI * cc).sqrt()) + (2 - 4 / PI) / W3   # HPY Lemma 2.6 beta\n    Q = (L + cc * W3)"),
 # blk_f float twin
 ("    be = 32 / (3 * np.sqrt(math.pi * c)) + (6 / math.pi - 1) / (c * W3 ** 2) + (1 - 6 / math.pi) / W3\n",
  "    be = 32 / (3 * np.sqrt(math.pi * c)) + (2 - 4 / math.pi) / W3\n"),
 # blk_exact_f: per-m bound X_m (HPY Lemma 2.5)
 ("+ 4 * math.sqrt(W / math.pi) / np.sqrt(m) + 1 - 6 / math.pi\n",
  "+ 4 * math.sqrt(W / math.pi) / np.sqrt(m) + 2 - 4 / math.pi\n"),
 # blk_exact_arb
 ("+ 4 * (W / PI).sqrt() * wmh + (1 - 6 / PI) * wz)", "+ 4 * (W / PI).sqrt() * wmh + (2 - 4 / PI) * wz)"),
 # tail beta_hat: W^{1/3} >= pi^{1/3}(r+1) (K >= t^{1/3})
 ("    return 32 / (3 * (PI * cc).sqrt()) + (6 / PI - 1) / (cc * PI.root(3) ** 2 * arb(r + 1) ** 2) + (1 - 6 / PI) / (PI.root(3) * (r + 1) * kappa)",
  "    return 32 / (3 * (PI * cc).sqrt()) + (2 - 4 / PI) / (PI.root(3) * (r + 1))"),
 # tail exact-Weyl S0c
 ("    S0c = 8 * (r32 * kiv ** 3).sqrt() * wmh + 2 * (1 - 6 / PI) * wz", "    S0c = 8 * (r32 * kiv ** 3).sqrt() * wmh + 2 * (2 - 4 / PI) * wz"),
 # tail r >= r1: Bhat
 ("    Bhat = 32 / (3 * (PI * cc).sqrt()) + (6 / PI - 1) / (cc * PI.root(3) ** 2 * arb(r1 + 1) ** 2)",
  "    Bhat = 32 / (3 * (PI * cc).sqrt()) + (2 - 4 / PI) / (PI.root(3) * (r1 + 1))"),
 # Part D uniform block (not run in ABCEF; replaced for completeness)
 ("    be = 32 / (3 * (PI * cc).sqrt()) + (6 / PI - 1) * Wm3 ** 2 / cc + (1 - 6 / PI) * Wm3",
  "    be = 32 / (3 * (PI * cc).sqrt()) + (2 - 4 / PI) * Wm3"),
 # float twins in tail_opt / tail_optM
 ("        return 32 / (3 * np.sqrt(math.pi * c)) + (6 / math.pi - 1) / (c * p3 ** 2 * (r + 1) ** 2) + (1 - 6 / math.pi) / (p3 * (r + 1) * kappa)",
  "        return 32 / (3 * np.sqrt(math.pi * c)) + (2 - 4 / math.pi) / (p3 * (r + 1))"),
 ("            S0c = 8 * r32 * np.sum(w / np.sqrt(m)) + 2 * (1 - 6 / math.pi) * np.sum(w)",
  "            S0c = 8 * r32 * np.sum(w / np.sqrt(m)) + 2 * (2 - 4 / math.pi) * np.sum(w)"),
]
for a, b in R:
    n = src.count(a)
    if n != 1: sys.exit("pattern count %d for: %s" % (n, a[:80]))
    src = src.replace(a, b)
assert '1 - 6 / PI' not in src and '1 - 6 / math.pi' not in src and '6 / PI - 1' not in src and '6 / math.pi - 1' not in src, 'leftover Revers constant'
os.makedirs('/tmp/revaudit_hpy', exist_ok=True)
open('/tmp/revaudit_hpy/certify_hpy.py', 'w').write("# HPY-only variant built by Riemann/Deep/RevAudit/make_hpy_variant.py\n" + src)
print("written /tmp/revaudit_hpy/certify_hpy.py with %d replacements" % len(R))
