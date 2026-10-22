# Lagarias / Salem — sufficiency directions (no sorry)

**Status 2026-10-22:** fixed. LagariasPositivity (Riemann/Tree/R009.lean) now assumes xi s ≠ 0, so the junk value ξ'/0 = 0 no longer enters, and the proof uses the pole of ξ'/ξ. Both modules are imported by Riemann.lean.

- `Suff.lean` (ns `Fayed.Lagarias`): `xi_eq` (ξ = s(s-1)/2·Λ, s≠0,1), `xi_zero_iff` (Re s>0, s≠1: ξ=0 ↔ ζ=0),
  `RH_of_LagariasPositivity : Fayed.R9.LagariasPositivity → RiemannHypothesis` (via `Fayed.R7.RH_iff_QuasiRH_half`).
  Lean proof uses only ξ'/0 = 0 at a zero; docstring notes the real content (pole of ξ'/ξ with residue m>0) is stronger.
- `Salem.lean` (ns `Fayed.Lagarias.Salem`): `salem_solution`: ζ(ρ)=0, 1/2<Re ρ<1 ⇒ φ(y)=exp(-i·Im ρ·y), |φ|=1, and
  ∀x, y ↦ e^{-σy}φ(y)/(exp(e^{x-y})+1) integrable with integral 0 (σ = Re ρ).
  Pieces: `mellin_K_eq` (∫₀^∞ t^{s-1}/(e^t+1) dt = Γ(s)η(s) on Re s>0, s≠1; η = (1-2^{1-s})ζ, `eta_eq`),
  via geometric series + Mathlib `hasSum_mellin` on Re s>1, even/odd split of ζ's series, identity theorem on
  {Re>0}\{1} (reuses `Fayed.Mertens.isPreconnected_U`), substitution `integral_comp_exp` + `integral_sub_left_eq_self`.
- Verify: `lake build Riemann.Deep.Lagarias.Suff`, `lake build Riemann.Deep.Lagarias.Salem` green;
  #print axioms of all main theorems = [propext, Classical.choice, Quot.sound].
- Caveat: these are the easy directions (hypothesis ⇒ RH); both hypotheses are open (≡ RH). Not wired into Riemann.lean.
