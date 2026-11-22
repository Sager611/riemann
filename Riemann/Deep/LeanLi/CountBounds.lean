module

public import Riemann.Deep.LeanLi.ZerosCount
public import Mathlib.Analysis.Real.Pi.Bounds

/-!
# LeanLi: Riemann–von Mangoldt counting bounds from the literature inputs
-/

@[expose] public section

open Complex Filter Topology Set ComplexConjugate

namespace LeanLi

/-- The main term `(T/2π) log(T/2πe)` of `N(T)`. -/
noncomputable def mainN (T : ℝ) : ℝ := T / (2 * Real.pi) * Real.log (T / (2 * Real.pi * Real.exp 1))

lemma mainN_eq {T : ℝ} (hT : 0 < T) :
    mainN T = T / (2 * Real.pi) * (Real.log (T / (2 * Real.pi)) - 1) := by
  unfold mainN
  rw [show T / (2 * Real.pi * Real.exp 1) = T / (2 * Real.pi) / Real.exp 1 by ring,
    Real.log_div (by positivity) (Real.exp_pos 1).ne', Real.log_exp]

/-- Increments of the main term. -/
lemma mainN_sub_le {s t : ℝ} (hs : 0 < s) (hst : s ≤ t) :
    mainN t - mainN s ≤ (t - s) / (2 * Real.pi) * Real.log (t / (2 * Real.pi)) := by
  have ht : 0 < t := hs.trans_le hst
  rw [mainN_eq ht, mainN_eq hs]
  have hπ : 0 < 2 * Real.pi := by positivity
  have hlog : Real.log (s / (2 * Real.pi)) = Real.log (t / (2 * Real.pi)) - Real.log (t / s) := by
    rw [Real.log_div (by positivity) hs.ne', Real.log_div ht.ne' hπ.ne', Real.log_div hs.ne' hπ.ne']
    ring
  have hts : Real.log (t / s) ≤ t / s - 1 := Real.log_le_sub_one_of_pos (by positivity)
  rw [hlog]
  have : s * Real.log (t / s) ≤ t - s := by
    have := mul_le_mul_of_nonneg_left hts hs.le
    rw [mul_sub, mul_div_cancel₀ _ hs.ne', mul_one] at this
    exact this
  rw [div_mul_eq_mul_div, div_mul_eq_mul_div, div_mul_eq_mul_div, ← sub_div, div_le_div_iff_of_pos_right hπ]
  nlinarith [this]

variable (L : LitCore)
include L

/-- Titchmarsh's averaged count: `N(T) = (T/2π) log(T/2πe) + O(log T)`. -/
lemma avg_asymp : ∃ C T₀ : ℝ, 3 ≤ T₀ ∧ ∀ T, T₀ ≤ T → |zeroCountAvg T - mainN T| ≤ C * Real.log T := by
  obtain ⟨C₁, T₁, h₁⟩ := L.rvm
  obtain ⟨C₂, t₂, h₂⟩ := L.rs_theta
  refine ⟨|C₁| + |C₂| + 1, max (max T₁ t₂) 3, le_max_right _ _, fun T hT => ?_⟩
  have hT1 : T₁ ≤ T := le_trans (le_trans (le_max_left _ _) (le_max_left _ _)) hT
  have hT2 : t₂ ≤ T := le_trans (le_trans (le_max_right _ _) (le_max_left _ _)) hT
  have hT3 : 3 ≤ T := le_trans (le_max_right _ _) hT
  have hTpos : 0 < T := by linarith
  have hlog1 : 1 ≤ Real.log T := by
    rw [← Real.log_exp 1]; exact Real.log_le_log (Real.exp_pos 1) (by linarith [Real.exp_one_lt_d9])
  have hS := h₁ T hT1
  have hθ := h₂ T hT2
  have hdef : zeroCountAvg T = zetaS T + rsTheta T / Real.pi + 1 := by unfold zetaS; field_simp; ring
  rw [hdef, mainN_eq hTpos]
  have hπ : 0 < Real.pi := Real.pi_pos
  -- rsTheta T / π = main + (-1/8) + err/π
  set E := rsTheta T - (T / 2 * (Real.log (T / (2 * Real.pi)) - 1) - Real.pi / 8) with hE
  have hθ' : rsTheta T / Real.pi = T / (2 * Real.pi) * (Real.log (T / (2 * Real.pi)) - 1) - 1 / 8 + E / Real.pi := by
    rw [hE]; field_simp; ring
  rw [hθ']
  have hEb : |E / Real.pi| ≤ |C₂| := by
    rw [abs_div, abs_of_pos hπ, div_le_iff₀ hπ]
    have : |E| ≤ C₂ / T := hθ
    have h3 : C₂ / T ≤ |C₂| / 3 := by
      calc C₂ / T ≤ |C₂| / T := div_le_div_of_nonneg_right (le_abs_self _) hTpos.le
        _ ≤ |C₂| / 3 := div_le_div_of_nonneg_left (abs_nonneg _) (by norm_num) hT3
    nlinarith [Real.pi_gt_three, abs_nonneg C₂]
  have key : zetaS T + (T / (2 * Real.pi) * (Real.log (T / (2 * Real.pi)) - 1) - 1 / 8 + E / Real.pi) + 1 -
      T / (2 * Real.pi) * (Real.log (T / (2 * Real.pi)) - 1) = zetaS T + 7 / 8 + E / Real.pi := by ring
  rw [key]
  calc |zetaS T + 7 / 8 + E / Real.pi| ≤ |zetaS T| + 7 / 8 + |E / Real.pi| := by
        have := abs_add_le (zetaS T + 7 / 8) (E / Real.pi)
        have := abs_add_le (zetaS T) (7 / 8)
        rw [abs_of_pos (by norm_num : (0:ℝ) < 7 / 8)] at this
        linarith
    _ ≤ |C₁| * Real.log T + 7 / 8 + |C₂| := by
        have : |zetaS T| ≤ |C₁| * Real.log T :=
          hS.trans (mul_le_mul_of_nonneg_right (le_abs_self _) (by linarith))
        linarith
    _ ≤ (|C₁| + |C₂| + 1) * Real.log T := by nlinarith [abs_nonneg C₂, abs_nonneg C₁]

/-- The jump of `N` at `T` (the multiplicity of zeros at height `T`) is `O(log T)`. -/
lemma jump_le : ∃ C T₀ : ℝ, 4 ≤ T₀ ∧ ∀ T, T₀ ≤ T →
    (zeroCount T : ℝ) - ∑ ρ ∈ zsetBelow T, (mult ρ : ℝ) ≤ C * Real.log T := by
  obtain ⟨C, T₀, hT₀, h⟩ := avg_asymp L
  refine ⟨3 * |C| + 1, T₀ + 1, by linarith, fun T hT => ?_⟩
  obtain ⟨t₀, ht₀, hn₀⟩ := exists_not_ordinate (a := T - 1) (b := T) (by linarith) (by linarith)
  obtain ⟨t₁, ht₁, hn₁⟩ := exists_not_ordinate (a := T) (b := T + 1) (by linarith) (by linarith)
  have hlow : (zeroCount t₀ : ℝ) ≤ ∑ ρ ∈ zsetBelow T, (mult ρ : ℝ) := by
    rw [zeroCount_eq]
    refine Finset.sum_le_sum_of_subset_of_nonneg (fun ρ hρ => ?_) (fun _ _ _ => Nat.cast_nonneg _)
    rw [mem_zset] at hρ; rw [mem_zsetBelow]
    exact ⟨hρ.1, hρ.2.1, lt_of_le_of_lt hρ.2.2 ht₀.2⟩
  have hup : (zeroCount T : ℝ) ≤ zeroCount t₁ := zeroCount_mono ht₁.1.le
  have e0 := avg_eq_of_not_ordinate hn₀
  have e1 := avg_eq_of_not_ordinate hn₁
  have h0 := h t₀ (by linarith [ht₀.1])
  have h1 := h t₁ (by linarith [ht₁.1])
  rw [e0] at h0; rw [e1] at h1
  have hm := mainN_sub_le (s := t₀) (t := t₁) (by linarith [ht₀.1]) (by linarith [ht₀.2, ht₁.1])
  have hTpos : 0 < T := by linarith
  have hlogT : 1 ≤ Real.log T := by
    rw [← Real.log_exp 1]; exact Real.log_le_log (Real.exp_pos 1) (by linarith [Real.exp_one_lt_d9])
  have hl0 : Real.log t₀ ≤ Real.log T := Real.log_le_log (by linarith [ht₀.1]) ht₀.2.le
  have hl1 : Real.log t₁ ≤ Real.log T + 1 := by
    have : Real.log t₁ ≤ Real.log (2 * T) := Real.log_le_log (by linarith [ht₁.1]) (by linarith [ht₁.2])
    rw [Real.log_mul (by norm_num) hTpos.ne'] at this
    linarith [Real.log_two_lt_d9]
  have hlogm : Real.log (t₁ / (2 * Real.pi)) ≤ Real.log T := by
    have h2π : 2 * Real.pi ≥ 2 := by linarith [Real.pi_gt_three]
    calc Real.log (t₁ / (2 * Real.pi)) ≤ Real.log ((T + 1) / 2) :=
          Real.log_le_log (div_pos (by linarith [ht₁.1]) (by positivity)) (div_le_div₀ (by linarith) ht₁.2.le (by norm_num) h2π)
      _ ≤ Real.log T := Real.log_le_log (by linarith) (by linarith)
  have hm' : mainN t₁ - mainN t₀ ≤ Real.log T := by
    have hd : (t₁ - t₀) / (2 * Real.pi) ≤ 1 := by
      rw [div_le_one (by positivity)]; linarith [ht₀.1, ht₁.2, Real.pi_gt_three]
    have hlogm0 : 0 ≤ Real.log (t₁ / (2 * Real.pi)) ∨ Real.log (t₁ / (2 * Real.pi)) < 0 := le_or_gt _ _
    rcases hlogm0 with hpos | hneg
    · calc mainN t₁ - mainN t₀ ≤ (t₁ - t₀) / (2 * Real.pi) * Real.log (t₁ / (2 * Real.pi)) := hm
        _ ≤ 1 * Real.log T := mul_le_mul hd hlogm hpos zero_le_one
        _ = Real.log T := one_mul _
    · have : (t₁ - t₀) / (2 * Real.pi) * Real.log (t₁ / (2 * Real.pi)) ≤ 0 :=
        mul_nonpos_of_nonneg_of_nonpos (div_nonneg (by linarith [ht₀.2, ht₁.1]) (by positivity)) hneg.le
      linarith
  have a0 := (abs_le.1 h0).1
  have a1 := (abs_le.1 h1).2
  have hC0 : C * Real.log t₀ ≤ |C| * Real.log T :=
    (mul_le_mul_of_nonneg_left hl0 (abs_nonneg C)) |>.trans' (by
      have : 0 ≤ Real.log t₀ := Real.log_nonneg (by linarith [ht₀.1])
      exact mul_le_mul_of_nonneg_right (le_abs_self C) this)
  have hC1 : C * Real.log t₁ ≤ |C| * (Real.log T + 1) :=
    (mul_le_mul_of_nonneg_left hl1 (abs_nonneg C)) |>.trans' (by
      have : 0 ≤ Real.log t₁ := Real.log_nonneg (by linarith [ht₁.1])
      exact mul_le_mul_of_nonneg_right (le_abs_self C) this)
  have hCl : |C| ≤ |C| * Real.log T := le_mul_of_one_le_right (abs_nonneg C) hlogT
  have e2 : (zeroCount t₁ : ℝ) - zeroCount t₀ ≤ mainN t₁ - mainN t₀ + C * Real.log t₁ + C * Real.log t₀ := by
    linarith [a0, a1]
  nlinarith [abs_nonneg C]

/-- `N(T) = (T/2π) log(T/2πe) + O(log T)` (right-continuous count). -/
lemma count_asymp : ∃ C T₀ : ℝ, 4 ≤ T₀ ∧ ∀ T, T₀ ≤ T → |(zeroCount T : ℝ) - mainN T| ≤ C * Real.log T := by
  obtain ⟨C, T₀, hT₀, h⟩ := avg_asymp L
  obtain ⟨C', T₁, hT₁, h'⟩ := jump_le L
  refine ⟨|C| + |C'|, max T₀ T₁, le_trans hT₁ (le_max_right _ _), fun T hT => ?_⟩
  have hA := h T (le_trans (le_max_left _ _) hT)
  have hJ := h' T (le_trans (le_max_right _ _) hT)
  have hlog : 0 ≤ Real.log T := Real.log_nonneg (by linarith [le_max_right T₀ T₁])
  have havg : zeroCountAvg T = ((zeroCount T : ℝ) + ∑ ρ ∈ zsetBelow T, (mult ρ : ℝ)) / 2 := by
    unfold zeroCountAvg; rw [zeroCountBelow_eq]; ring
  have hle : ∑ ρ ∈ zsetBelow T, (mult ρ : ℝ) ≤ zeroCount T := by
    rw [zeroCount_eq]
    exact Finset.sum_le_sum_of_subset_of_nonneg (zsetBelow_subset T) (fun _ _ _ => Nat.cast_nonneg _)
  rw [havg] at hA
  have hCl : C * Real.log T ≤ |C| * Real.log T := mul_le_mul_of_nonneg_right (le_abs_self C) hlog
  have hC'l : C' * Real.log T ≤ |C'| * Real.log T := mul_le_mul_of_nonneg_right (le_abs_self C') hlog
  rw [abs_le] at hA ⊢
  constructor <;> nlinarith [hA.1, hA.2]

/-- `N(T) ≤ T log T` for large `T`. -/
lemma count_le : ∃ T₀ : ℝ, 4 ≤ T₀ ∧ ∀ T, T₀ ≤ T → (zeroCount T : ℝ) ≤ T * Real.log T := by
  obtain ⟨C, T₀, hT₀, h⟩ := count_asymp L
  refine ⟨max T₀ (2 * |C| + 4), le_trans hT₀ (le_max_left _ _), fun T hT => ?_⟩
  have hT0 : T₀ ≤ T := le_trans (le_max_left _ _) hT
  have hTC : 2 * |C| + 4 ≤ T := le_trans (le_max_right _ _) hT
  have hTpos : 0 < T := by linarith [abs_nonneg C]
  have hlog : 0 ≤ Real.log T := Real.log_nonneg (by linarith [abs_nonneg C])
  have h1 := (abs_le.1 (h T hT0)).2
  rw [mainN_eq hTpos] at h1
  have hπ : (1 : ℝ) ≤ 2 * Real.pi := by linarith [Real.pi_gt_three]
  have hl : Real.log (T / (2 * Real.pi)) - 1 ≤ Real.log T := by
    have : Real.log (T / (2 * Real.pi)) ≤ Real.log T :=
      Real.log_le_log (by positivity) (div_le_self hTpos.le hπ)
    linarith
  have hm : T / (2 * Real.pi) * (Real.log (T / (2 * Real.pi)) - 1) ≤ T / 2 * Real.log T := by
    rcases le_total 0 (Real.log (T / (2 * Real.pi)) - 1) with h0 | h0
    · calc T / (2 * Real.pi) * (Real.log (T / (2 * Real.pi)) - 1) ≤ T / 2 * (Real.log (T / (2 * Real.pi)) - 1) :=
            mul_le_mul_of_nonneg_right (div_le_div_of_nonneg_left hTpos.le (by norm_num) (by linarith [Real.pi_gt_three])) h0
        _ ≤ T / 2 * Real.log T := mul_le_mul_of_nonneg_left hl (by linarith)
    · have : T / (2 * Real.pi) * (Real.log (T / (2 * Real.pi)) - 1) ≤ 0 :=
        mul_nonpos_of_nonneg_of_nonpos (by positivity) h0
      nlinarith
  have hC : C * Real.log T ≤ T / 2 * Real.log T :=
    mul_le_mul_of_nonneg_right (by linarith [le_abs_self C]) hlog
  linarith

end LeanLi
