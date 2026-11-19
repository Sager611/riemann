"""padhi2 side check: can Padhi's certified scan reach past log x = 3159 (Axler 2022 needs the faulty Dudek-Platt case
below about e^3160.3)? q_arb on [3159, 3161] and [3161, 3163] with his top kernels (c = 70, f = 0.40), published N3."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import common as CM
S = CM.S
for L1, L2, c, f in [(3150, 3159, 70, 0.40), (3159, 3161, 70, 0.40), (3161, 3163, 70, 0.40)]:
    q, parts = S.q_arb(L1, L2, 3.0e12, c, f)
    print(L1, L2, 'q_upper = %.9f' % S.fup_safe(q), parts, flush=True)
