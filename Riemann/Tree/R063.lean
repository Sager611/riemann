/-!
# Round 63 — zeta23 reuse map (Anthropic formal-math, Mathlib v4.33, Apache-2.0)

| subtree | key theorems | unlocks for us | port effort |
|---------|--------------|----------------|-------------|
| WeilEF | rectangle_identity, full_line_identity, EF_lit_zeta | explicit formula ⇒ Weil positivity ⇒ RH, Li | ~4.5k lines, high drift risk |
| ExplicitFormula | paperFT_weilTest, literatureRHS_eq_integral_nu | normalisation / integrability | ~0.9k |
| RvM | rectangleIntegral'_logDeriv_completedZeta_eq_Ncount, N_eq_halfContour_completedZeta | N(T), argument principle | ~2k |
| FromPNTPlus | ResidueTheoremInRectangle | residues on rectangles (Speiser, counting) | core files of ~10.8k |
| GammaFacts | StirlingVert.mu_stirling, Complete.gammaFacts | Stirling in strips | ~1.25k |
| XiPrime | hardyW_ne_zero_of_one_le_re, riemannVonMangoldt_W | Hardy Z, zero-free strip templates | ~2k |
| Unconditional / Statement | two_thirds_on_critical_line (no hypotheses, standard axioms) | — | — |
| ZeroSide / Analytic / Poisson | weighted argument principle, Poisson trace bounds | zero counting | ~1.5k |

Plan: port ResidueTheoremInRectangle + WeilEF to v4.35, then formalize "Weil positivity on all windows ⇒ RH"
on top of Riemann/Deep/Weil/Defs.lean (queued as an Orca deep dive).
-/
