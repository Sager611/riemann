module

public import Riemann.Deep.Combos.Lean.ZetaA

@[expose] public section

/-! # ζ(1 + r) ≤ 1 + 1/r: the telescoped sum (NOTE l.89) -/

open Finset

namespace RingClass

/-- Σ_{k=1}^{M} k^{−(1+r)} ≤ 1 + (1 − M^{−r})/r -/
lemma sum_Icc_le (r : ℝ) (hr : 0 < r) (M : ℕ) (hM : 1 ≤ M) :
    ∑ k ∈ Icc 1 M, (((k : ℕ) : ℝ) ^ (1 + r))⁻¹ ≤ 1 + (1 - (((M : ℕ) : ℝ) ^ r)⁻¹) / r := by
  induction M, hM using Nat.le_induction with
  | base => simp
  | succ n hn ih =>
    rw [Finset.sum_Icc_succ_top (by omega)]
    have hn1 : (1 : ℝ) ≤ (n : ℝ) := by exact_mod_cast hn
    have ht := tele_step r (n : ℝ) hr hn1
    have e : ((n + 1 : ℕ) : ℝ) = (n : ℝ) + 1 := by push_cast; ring
    rw [e]
    have e2 : 1 + (1 - ((n : ℝ) ^ r)⁻¹) / r + (((n : ℝ) ^ r)⁻¹ - (((n : ℝ) + 1) ^ r)⁻¹) / r
        = 1 + (1 - (((n : ℝ) + 1) ^ r)⁻¹) / r := by ring
    linarith

/-- **ζ(1 + r) ≤ 1 + 1/r**: Σ_{k∈K} k^{−(1+r)} ≤ 1 + 1/r for every finite K ⊂ ℕ_{>0}. -/
theorem sum_rpow_neg_le (r : ℝ) (hr : 0 < r) (K : Finset ℕ) (hK : ∀ k ∈ K, 0 < k) :
    ∑ k ∈ K, (((k : ℕ) : ℝ) ^ (1 + r))⁻¹ ≤ 1 + 1 / r := by
  set M := K.sup id with hMdef
  have hsub : K ⊆ Icc 1 (max M 1) := by
    intro k hk
    rw [Finset.mem_Icc]
    refine ⟨hK k hk, le_trans ?_ (le_max_left _ _)⟩
    exact Finset.le_sup (f := id) hk
  have h1 : ∑ k ∈ K, (((k : ℕ) : ℝ) ^ (1 + r))⁻¹ ≤ ∑ k ∈ Icc 1 (max M 1), (((k : ℕ) : ℝ) ^ (1 + r))⁻¹ :=
    Finset.sum_le_sum_of_subset_of_nonneg hsub (fun k _ _ => by positivity)
  have h2 := sum_Icc_le r hr (max M 1) (le_max_right _ _)
  have h3 : 0 ≤ ((((max M 1 : ℕ)) : ℝ) ^ r)⁻¹ := by positivity
  have h4 : (1 - ((((max M 1 : ℕ)) : ℝ) ^ r)⁻¹) / r ≤ 1 / r := by
    apply div_le_div_of_nonneg_right _ hr.le; linarith
  linarith

end RingClass

