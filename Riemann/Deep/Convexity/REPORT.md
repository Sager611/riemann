# Convexity bound — Fayed `Inputs.convexity` discharged (sorry-free)

## Result
Ladder rung (2) done: the exact field `Inputs.convexity` is a theorem from Mathlib alone.
Sharp exponent η/2 + ε on Re s = 1-η, so rung (3) is covered too.
Stretch done: `fayed_RH_noConvexity`, i.e. Fayed's argument with convexity removed from the inputs.

## Lean names (Riemann/Deep/Convexity/)
- `Convexity.lean`: `Fayed.Convexity.convexity_holds` (same statement as the field; this is checked
  when `inputs_of : InputsNoConvexity → Inputs` is built), `fayed_RH_noConvexity`.
- `Partial.lean`: `zeta_eq_on_U`: ζ(s) = s/(s-1) - s∫₁^∞{t}t^{-s-1}dt on Re>0, Im>0.
  Comes from Mathlib `LSeries_eq_mul_integral'` plus the identity theorem; the integral is
  written as a Mellin transform to get holomorphy. `norm_zeta_le_on_U`: ‖ζ‖ ≤ ‖s‖/‖s-1‖ + ‖s‖/σ.
- `FE.lean`: `left_edge`: ‖ζ(it)‖ ≤ √(4t/25)‖ζ(1+it)‖. Uses the exact value |Γ(1-it)|² = πt/sinh πt
  (from the reflection formula) and `riemannZeta_one_sub`.
  `apriori`: ‖ζ(σ+it)‖ ≤ 24e^{3t} for 0≤σ≤2, t≥1. `norm_Gamma_le`: |Γ(w)| ≤ Γ(Re w).
- `Bootstrap.lean`, `MaxMod.lean`: g = ζ(s)(-is)^{(s-1-δ)/2}e^{εs²} on the half-strip
  0≤σ≤1+δ, t≥1. Edges: left |g(it)| ≤ (4/5)|g(1+it)|; right bounded by ζ(1+δ);
  the Gaussian factor makes the top edge small.
  `step`: max modulus on rectangles [0,1+δ]×[1,T] (`Complex.norm_le_of_forall_mem_frontier_norm_le`).
  `iterate_bound` then gives |g| ≤ 5K with K independent of ε (`gf_bound`).
  `convexity_pos`: take ε = 1/t² and δ = min(1, 2ε); t<0 follows from `riemannZeta_conj`.

## Method notes
- Phragmén–Lindelöf is not used. The pole at s=1 is avoided by working only on Im s ≥ 1.
- No bound on Re s = 1 is needed: the left edge is tied to the interior line Re s = 1 and closed by
  iteration with factor 4/5 < 1.
- No complex Stirling formula: only exact |Γ| values on Re = 1, and |Γ(w)| ≤ Γ(Re w).

## Axioms
`#print axioms` for convexity_holds and fayed_RH_noConvexity: [propext, Classical.choice, Quot.sound].
No sorry, no new axiom (`rg sorry Riemann/Deep/Convexity` finds nothing).

## Gaps
None for this field. Rung (1), polynomial growth on any compact σ-range, is not stated as its own
lemma. It follows from the same machinery: on σ ≥ σ₀ > 0 from `norm_zeta_le_on_U`, and with
convexity in the strip.

## Build
`lake build Riemann.Deep.Convexity.Convexity` is green. Not yet imported by `Riemann.lean`.
