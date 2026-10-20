# Riesz criterion, sufficiency — `Riemann/Deep/Riesz/Suff.lean`

**Status:** sorry-free, no extra hypotheses. `lake build Riemann.Deep.Riesz.Suff` is green.
`#print axioms Fayed.Riesz.RH_of_rieszBound` → `[propext, Classical.choice, Quot.sound]`.

**Main theorem:** `Fayed.Riesz.RH_of_rieszBound (h : ∀ ε : ℝ, 0 < ε → ∃ C, ∀ x ≥ 1, |rieszM x| ≤ C * x^(1/4+ε)) : RiemannHypothesis`,
where `rieszM x = ∑' n, μ n * (x/n²) * exp(-(x/n²))`.

**Proof steps (lemma names):**
1. `F_eq`: rieszM t = t · G(t), with G(t) = Σ (μ(n)/n²) e^{-t/n²}.
2. `mellin_eq`: M(s) := mellin F (-s) = Γ(1-s) · L(μ, 2s) for 1/2 < Re s < 1. Uses Mathlib's `hasSum_mellin` and `mellin_cpow_smul`; the summability side condition is `summable_coeff`.
3. `differentiableAt_M`: M is holomorphic on 1/4 < Re s < 1, from `mellin_differentiableAt_of_isBigO_rpow`. It needs three inputs: the bound gives O(x^{1/4+ε}) at ∞; `abs_rieszM_le` gives |rieszM x| ≤ ζ(2)·x near 0; `continuousOn_rieszM` gives local integrability.
4. `eq_on_convex`: the identity theorem gives ζ(2s) M(s) = Γ(1-s). It is applied on two convex strips that avoid s = 1/2: {Im − Re > −1/2} and {Re + Im < 1/2}. This avoids proving that a punctured strip is connected.
5. `zeta_two_mul_ne_zero`: Γ ≠ 0 implies ζ(w) ≠ 0 for 1/2 < Re w < 1. The final step uses `Fayed.R7.RH_iff_QuasiRH_half` (the functional equation).

**Not done:** the link to `Fayed.R13.rieszSeries` (Σ (−x)^{k+1}/(k! ζ(2k+2)) = rieszM x) appears only in the module docstring and is not proved.
