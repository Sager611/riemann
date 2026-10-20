# Mertens bound ⇒ zero-free half-plane — REPORT

Status: DONE. Sorry-free. Axioms: [propext, Classical.choice, Quot.sound] (both files).

## Proved
- `Fayed.R2.mertensZeroFree : Fayed.R2.MertensZeroFree` (Riemann/Deep/Mertens/ZeroFree.lean;
  internal proof `Fayed.Mertens.mertensZeroFree`). Target statement unchanged, for every real θ
  (θ < 0 needs no special case: the same argument works, the hypothesis is just never true).
- Companion (Riemann/Deep/Mertens/Converge.lean): `Fayed.Mertens.tendsto_partialSum` and
  `mertens_LSeries_converges`: MertensBound θ, 0 ≤ θ < Re s ⇒ Σ_{n≤N} μ(n)/n^s → g(s).

## Proof route
1. Mstep t := M(⌊t⌋) (complex). |M(⌊t⌋)| = O(t^θ) at ∞ (any sign of θ, via
   `isEquivalent_nat_floor` + `IsEquivalent.rpow`); Mstep = 0 on (0,1); locally integrable.
2. g(s) := s · mellin Mstep (−s) = s ∫_1^∞ M(t) t^{−s−1} dt; holomorphic on Re s > θ by
   `mellin_differentiableAt_of_isBigO_rpow` (no hand-rolled dominated differentiation needed).
3. Re s > max θ 1: `LSeries_eq_mul_integral` (r = max θ 0) + `LSeriesSummable_moebius_iff` give
   L(μ,s) = g(s); `LSeries_zeta_mul_Lseries_moebius`, `LSeries_zeta_eq_riemannZeta` give ζ·g = 1.
4. U = {Re s > θ} \ {1} is open and preconnected (hand proof: union of 4 convex pieces via
   `isPreconnected_of_forall`); ζ·g analytic on U (`differentiableAt_riemannZeta`,
   `DifferentiableOn.analyticOnNhd`); identity theorem `eqOn_of_preconnected_of_eventuallyEq`
   ⇒ ζ·g = 1 on U ⇒ ζ ≠ 0 there. s = 1: `riemannZeta_one_ne_zero`.

## Other Mathlib lemmas
`setIntegral_eq_of_subset_of_ae_sdiff_eq_zero`, `integral_Ici_eq_integral_Ioi`,
`Nat.measurable_floor`, `LocallyIntegrableOn.mono`, `abs_moebius_le_one`, `Complex.norm_intCast`;
companion copies the body of Mathlib's private `LSeries_eq_mul_integral_aux`
(`tendsto_sum_mul_atTop_nhds_one_sub_integral₀`) minus the summability step.

## Gaps
None for the target. Companion needs 0 ≤ θ (for θ < 0 the bound is false anyway, since
|M(N)| < 1 eventually forces μ(p) = 0 at large primes; not formalized).
Not wired into Riemann.lean / TREE.md (outside my allowed paths) — orchestrator to add imports.

## Build
`lake build Riemann.Deep.Mertens.ZeroFree Riemann.Deep.Mertens.Converge`

## Progress log
- ZeroFree.lean compiled green, axioms clean (~50 s first compile, ~16 s rebuild).
- Converge.lean (companion) compiled green, axioms clean.
