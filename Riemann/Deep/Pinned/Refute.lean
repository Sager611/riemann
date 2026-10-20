import Mathlib

/-!
# Refutations of pinned false statements (R075 follow-up)

* **Robin-threshold, Zenodo 21869761 (M. Bitto), §5 eqs. (17)–(19).** The paper bounds
  `q₁ > ∏_{p ≤ pₙ} (1 - p⁻³)⁻¹`, takes the limit (19) `= ζ(3)`, and then asserts
  "it can be stated that q₁ > ζ(3)", i.e. it uses `∏_{p ≤ pₙ} (1 - p⁻³)⁻¹ > ζ(3)` for every prime
  `pₙ`. Refuted at `pₙ = 13`: the product is `1.20133…` while `ζ(3) ≥ ∑_{n ≤ 40} n⁻³ = 1.20171…`.

Skipped (no cheap explicit counterexample):
* Q(A), Zenodo 18923112, Thm 4.8 / eq. (38): the identity itself is correct (summation by parts of
  `∑ (S(a)² - S(a-1)²)/a²`); the error is interpretive (annihilating the diagonal does not bound
  the intermediate band), not a false finite statement.
* Vélez, Zenodo 21173721, Lemma 7.1: a statement about tempered-distribution limits of semilocal
  spectral measures; given Prop. 6.1 (`q_X ≥ 0`) the positivity of the limit is actually valid, so the
  real gap is upstream and RH-strength — no finite witness.
* Li-criterion, Zenodo 16900110: the deposited PDF has one page (definitions (1)–(4) only); no
  false claim is stated in the available text.
-/

namespace Fayed.Pinned

open Finset

/-- Bitto, Zenodo 21869761, §5 ((17)–(19) and the sentence after): for every prime `P`,
the partial Euler product `∏_{p ≤ P} (1 - p⁻³)⁻¹` exceeds `ζ(3)`. -/
def RobinEq19 : Prop :=
  ∀ P : ℕ, P.Prime →
    (riemannZeta 3).re < ∏ p ∈ (range (P + 1)).filter Nat.Prime, (1 - ((p : ℝ) ^ 3)⁻¹)⁻¹

/-- `ζ(3)` as a real series. -/
theorem zeta_three_re : (riemannZeta 3).re = ∑' n : ℕ, 1 / (n : ℝ) ^ 3 := by
  have h := zeta_nat_eq_tsum_of_gt_one (k := 3) (by norm_num)
  push_cast at h
  rw [h]
  have h2 : (∑' n : ℕ, 1 / (n : ℂ) ^ 3) = ((∑' n : ℕ, 1 / (n : ℝ) ^ 3 : ℝ) : ℂ) := by
    rw [Complex.ofReal_tsum]; push_cast; rfl
  rw [h2, Complex.ofReal_re]

/-- `ζ(3) ≥ ∑_{n ≤ 40} n⁻³`. -/
theorem partial_le_zeta_three :
    ∑ n ∈ range 41, 1 / (n : ℝ) ^ 3 ≤ (riemannZeta 3).re := by
  rw [zeta_three_re]
  exact Summable.sum_le_tsum _ (fun _ _ => by positivity)
    (Real.summable_one_div_nat_pow.mpr (by norm_num))

theorem primes_le_13 : (range (13 + 1)).filter Nat.Prime = {2, 3, 5, 7, 11, 13} := by
  ext n
  simp only [mem_filter, mem_range, mem_insert, mem_singleton]
  constructor
  · rintro ⟨h1, h2⟩
    interval_cases n <;> first | (norm_num at h2; done) | norm_num
  · rintro (rfl | rfl | rfl | rfl | rfl | rfl) <;> norm_num

/-- The partial Euler product at `P = 13` is `< ζ(3)` (`1.2013314… < 1.2017…≤ ζ(3)`). -/
theorem euler13_lt_zeta_three :
    ∏ p ∈ (range (13 + 1)).filter Nat.Prime, (1 - ((p : ℝ) ^ 3)⁻¹)⁻¹ < (riemannZeta 3).re := by
  refine lt_of_lt_of_le ?_ partial_le_zeta_three
  rw [primes_le_13]
  simp only [sum_range_succ, sum_range_zero]
  norm_num [prod_insert, prod_singleton]

/-- Bitto's eq. (19) step "q₁'s lower bound exceeds ζ(3)" is false (witness `P = 13`). -/
theorem not_robinEq19 : ¬ RobinEq19 := fun h =>
  lt_asymm (h 13 (by norm_num)) euler13_lt_zeta_three

end Fayed.Pinned

#print axioms Fayed.Pinned.not_robinEq19
