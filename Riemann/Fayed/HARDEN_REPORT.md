# HARDEN report

## Milestones
- Split Paper.lean: Defs.lean (definitions, Inputs, Setting, algebraic identities, G lemmas). Built green.
- Asymptotics.lean: J at a zero, identity (a), bounds for F, F' and the boundary terms. Built green (138 s).
- Paper.lean: eq39, eq42, eq50, eq54 proved; fayed_RH and fayed_gap unchanged. Built green.
- StripReflect.lean: Inputs.strip and Inputs.reflect proved from Mathlib (strip_holds, reflect_holds).
- lake build green; Audit: every theorem depends only on [propext, Classical.choice, Quot.sound].

## Final summary
All four HARDEN sorries are gone. Names and statements of eq39, eq42, eq50, eq54, fayed_RH,
fayed_gap and every Inputs field are unchanged. No statement was false.

How each step is proved (all from the Inputs hypotheses):
- (34) at a zero (J_at_zero): perron33 (ε = 1/2) + em17 give
  ‖u^s J(u,s) + u^s EMcorr(u,s)‖ ≤ C/u on D_γ, for s₀ and for 1 - s̄₀ (a zero by Inputs.reflect).
- Identity (a) (Ffun_eq): F(u) = u^{s₀}J(u,s₀) - u^{1-s̄₀}J(u,1-s̄₀). Both sides are integrals of
  the same function Kf(y)/A(y) - Kf(y)/B(y); the integrands are continuous, so the integral splits.
- eq39 / eq42: the B̄₁ terms cancel, so ‖F‖ ≤ K/u and ‖F/G‖ ≤ K u^{p-1} → 0 (FG_bound);
  F = (F/G)·G → 0.
- eq50: Fder_eq writes F' = u⁻¹(s₀u^{s₀}J(s₀) - s₁u^{s₁}J(s₁)) + boundary (eq44_algebra under the
  integral). X_bound: the first part is (1-2σ₀)B̄₁(u) + O(1/u). bdry_bound: the boundary is
  O(u^{-2}) (1/A - 1/B = (2σ₀-1)/(AB), |AB| ≥ u^16, convexity with ε = 1/4). On D_γ, B̄₁(u) > γ,
  so ‖F'/G'‖ ≥ c·u^p → ∞.
- eq54: for a_M ≥ U, the bound of FG_bound on the wider domain D_{γ/2} (which also contains the
  right endpoints b_M) gives δ_M = 2K a_M^{p-1} → 0; for the finitely many small M, F/G is
  continuous on [a_M, b_M] (Inputs.leibniz44), hence bounded.
- Bounds |B̄₂|, |B̄₃| ≤ c: Bernoulli polynomials are continuous on [0,1] (compactness).
  B̄₁ = {x} - 1/2 from Polynomial.bernoulli_one.

Optional part: strip and reflect are now theorems of Mathlib (functional equation, conjugation,
no zeros on Re s ≥ 1). The fields stay in Inputs because the brief fixes that interface.
So Inputs reduces to four real assumptions: em17, perron33, leibniz44, convexity.

Consequence: given Inputs, the paper's argument is checked step by step except (55), and
fayed_gap pins the failure on Theorem 1 (refuted in Theorem1.lean).

Files: Riemann/Fayed/{Defs,Asymptotics,Paper,StripReflect,Audit}.lean. Riemann.lean imports
StripReflect (Defs and Asymptotics come in through Paper).
