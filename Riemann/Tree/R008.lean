import Mathlib

/-!
# Round 8 — invention: new mechanisms built to pass both filters

Filters: (A) σ-blindness, (B) Davenport–Heilbronn (FE, off-line zeros, no Euler product).

| id | invented mechanism | status |
|----|--------------------|--------|
| R8.1 | entropy deficit of the prime-log law | CIRCULAR: zero-weighted lower bound encodes RH |
| R8.2 | Gaussian multiplicative chaos of the Euler product | DEAD: moments see prime variance, not zeros |
| R8.3 | Euler–Choi complete positivity | CIRCULAR: Weil positivity in Choi form |
| R8.4 | Euler–Schur GCD kernel positivity | DEAD: PSD for every σ, fails filter A (gcd_kernel_psd) |
| R8.5 | prime heat energy monotonicity | DEAD: decay automatic for Λ ≥ 0 (prime_heat_antitone) |
| R8.6 | Morse flux of log ζ | CIRCULAR: vanishing defect ≡ no off-line zeros |
| R8.7 | additive energy of log p | CIRCULAR: zero–energy bridge is RH |
| R8.8 | Möbius joinings ⇒ uniform ζ'/ζ bound | DEAD: bound false (universality: ζ'/ζ unbounded on σ ∈ (1/2,1)) |
| R8.9 | prime-constrained Weil SDP | CIRCULAR: strict certificate ≡ Weil positivity |
| R8.10 | prime-phase energy (mean values) | DEAD: Montgomery–Vaughan mean values are unconditional, zero-blind |
-/

namespace Fayed.R8

/-- R8.4 no-go: a GCD kernel with nonnegative weights is positive semidefinite for any
weights — in particular for w(d) = Λ(d) d^{-σ} at every σ — so it carries no σ-information. -/
theorem gcd_kernel_psd (N : ℕ) (D : Finset ℕ) (w : ℕ → ℝ) (hw : ∀ d ∈ D, 0 ≤ w d)
    (x : ℕ → ℝ) :
    0 ≤ ∑ i ∈ Finset.range N, ∑ j ∈ Finset.range N,
      x i * x j * ∑ d ∈ D, w d * (if d ∣ i ∧ d ∣ j then (1 : ℝ) else 0) := by
  set a : ℕ → ℕ → ℝ := fun d i => if d ∣ i then x i else 0 with ha
  have h1 : ∀ i j, x i * x j * ∑ d ∈ D, w d * (if d ∣ i ∧ d ∣ j then (1 : ℝ) else 0)
      = ∑ d ∈ D, w d * (a d i * a d j) := by
    intro i j
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl (fun d _ => ?_)
    simp only [ha]
    by_cases hi : d ∣ i <;> by_cases hj : d ∣ j <;> simp [hi, hj] <;> ring
  have h2 : ∀ d, w d * (∑ i ∈ Finset.range N, a d i) ^ 2
      = ∑ i ∈ Finset.range N, ∑ j ∈ Finset.range N, w d * (a d i * a d j) := by
    intro d
    rw [sq, Finset.sum_mul_sum, Finset.mul_sum]
    refine Finset.sum_congr rfl (fun i _ => ?_)
    rw [Finset.mul_sum]
  calc (0 : ℝ) ≤ ∑ d ∈ D, w d * (∑ i ∈ Finset.range N, a d i) ^ 2 :=
        Finset.sum_nonneg (fun d hd => mul_nonneg (hw d hd) (sq_nonneg _))
    _ = ∑ d ∈ D, ∑ i ∈ Finset.range N, ∑ j ∈ Finset.range N, w d * (a d i * a d j) :=
        Finset.sum_congr rfl (fun d _ => h2 d)
    _ = ∑ i ∈ Finset.range N, ∑ j ∈ Finset.range N, ∑ d ∈ D, w d * (a d i * a d j) := by
        rw [Finset.sum_comm]
        refine Finset.sum_congr rfl (fun i _ => ?_)
        rw [Finset.sum_comm]
    _ = _ := Finset.sum_congr rfl (fun i _ => Finset.sum_congr rfl (fun j _ => (h1 i j).symm))

/-- R8.5 no-go: heat-type energies of nonnegative coefficients decrease in t automatically,
so monotonicity carries no zero information. -/
theorem prime_heat_antitone (N : ℕ) (a : ℕ → ℝ) (ha : ∀ n, 0 ≤ a n) (σ : ℝ) {t u : ℝ}
    (htu : t ≤ u) :
    ∑ n ∈ Finset.range N, a n * Real.exp (-u * (Real.log n) ^ 2 - σ * Real.log n)
      ≤ ∑ n ∈ Finset.range N, a n * Real.exp (-t * (Real.log n) ^ 2 - σ * Real.log n) := by
  refine Finset.sum_le_sum (fun n _ => mul_le_mul_of_nonneg_left ?_ (ha n))
  apply Real.exp_le_exp.2
  have := sq_nonneg (Real.log n)
  nlinarith

end Fayed.R8
