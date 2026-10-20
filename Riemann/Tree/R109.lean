/-!
# Round 109 — ideas to break the 13/15 floor of the OpenAI/Liu positivity system

Floor (Riemann/Deep/QRHImprove/REPORT.md): the upper-bound side forces B ≥ 13/15 ≈ 0.8667 via the
reflected-energy term (5ℓ − 1)₊/4 once the detector length ℓ > 1/5.

| idea | verdict |
|------|---------|
| mollify the detector (Levinson/Conrey) | speculative: cross terms may restore the floor |
| two-piece detector (short + long, ℓ ≤ 1/5 on the reflected part) | speculative: cross terms |
| higher-degree CM field instead of Q(√−3) | speculative: no reason the structural term disappears |
| GL(2) / symmetric-square auxiliary family | speculative |
| Selberg-sieve weights in the kernel | speculative: constants only |
| bootstrap the 0.8749 half-plane back in | dead: floor is structural |
| Guth–Maynard density only in (13/15, 7/8) | dead: floor is pointwise |
| reflect to σ < 1/2 | dead: destroys positivity |

Most promising to test numerically: the two-piece detector (cheapest to add to refined_system.py).
-/
