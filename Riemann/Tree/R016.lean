import Mathlib

/-!
# Round 16 — personas: best unexplored direction per school

Every persona produced an infrastructure milestone; none carries an RH bridge.

| id | school | direction / first milestone | status |
|----|--------|-----------------------------|--------|
| R16.1 | analytic NT | hybrid large sieve → uniform zero detection | INFRA: density ≠ exclusion |
| R16.2 | NCG | adelic scaling representation + trace formula | INFRA: positivity remains RH |
| R16.3 | arithmetic geometry | polarized H¹(Spec ℤ) with Frobenius flow | INFRA: cohomology missing (polarized_eigen) |
| R16.4 | homogeneous dynamics | Weil form as Hecke-packet variance | INFRA: uniformity is the gap |
| R16.5 | RMT / probability | negative fractional moments of ζ off the line | INFRA: not zero-excluding |
| R16.6 | de Branges | canonical system for Ξ, Hamiltonian positivity | CIRCULAR: positivity ≡ RH |
| R16.7 | physics | PT-symmetric scattering determinant | INFRA: arithmetic realization missing |
| R16.8 | formal/AI | exact SOS certificates for Weil kernels | INFRA: finite only (gram_nonneg) |
| R16.9 | real-rootedness | Lorentzian arrays for Jensen polynomials | DEAD: degree-2 milestone already known (Csordas–Norfolk–Varga 1986) |
| R16.10 | model theory | uniform o-minimal complexity on strips | DEAD |
-/

namespace Fayed.R16

/-- R16.3: in a polarized setting (‖F v‖ = c‖v‖) every eigenvalue has modulus c — the purity
mechanism of Deligne/Weil; for ζ the space and F are what is missing. -/
theorem polarized_eigen {V : Type*} [NormedAddCommGroup V] [NormedSpace ℂ V] (F : V →ₗ[ℂ] V)
    (c : ℝ) (hF : ∀ v, ‖F v‖ = c * ‖v‖) {μ : ℂ} {v : V} (hv : v ≠ 0) (h : F v = μ • v) :
    ‖μ‖ = c := by
  have h1 := hF v
  rw [h, norm_smul] at h1
  exact mul_right_cancel₀ (norm_pos_iff.2 hv).ne' h1

/-- R16.8: a Gram (sum-of-squares) certificate gives a nonnegative form — finite positivity is free;
the analytic limit is the hard part. -/
theorem gram_nonneg {n m : ℕ} (B : Fin n → Fin m → ℝ) (x : Fin m → ℝ) :
    0 ≤ ∑ i, (∑ j, B i j * x j) ^ 2 :=
  Finset.sum_nonneg (fun i _ => sq_nonneg _)

end Fayed.R16
