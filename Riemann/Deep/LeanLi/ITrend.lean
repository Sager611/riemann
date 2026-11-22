module

public import Riemann.Deep.LeanLi.IClosed
public import Mathlib.NumberTheory.Harmonic.EulerMascheroni

/-!
# LeanLi: `I(n) = t(n) + O(1)`, `t(n) = (n/2)(log n - 1 - log 2π + γ)`
-/

@[expose] public section

open Real Filter Topology Set

namespace LeanLi

/-- `t(n)` (same expression as `liTrend` of `Defs.lean`). -/
noncomputable def trendT (n : ℕ) : ℝ :=
  (n / 2 : ℝ) * (log n - 1 - log (2 * π) + eulerMascheroniConstant)

/-- `Σ_{k odd < n} 1/k`. -/
noncomputable def oddRecip (n : ℕ) : ℝ := ∑ k ∈ (Finset.range n).filter Odd, (1 : ℝ) / k

lemma Iclosed_eq (n : ℕ) :
    Iclosed n = n * oddRecip n - ((Finset.range n).filter Odd).card - (n / 2) * log (4 * π) := by
  unfold Iclosed oddRecip
  rw [Finset.mul_sum, Finset.card_eq_sum_ones, Nat.cast_sum, ← Finset.sum_sub_distrib]
  congr 1
  refine Finset.sum_congr rfl fun k hk => ?_
  have hk0 : (k : ℝ) ≠ 0 := by
    have := (Finset.mem_filter.1 hk).2
    rcases this with ⟨r, rfl⟩; push_cast; positivity
  field_simp; push_cast; ring

lemma filter_odd_succ (n : ℕ) : (Finset.range (n + 1)).filter Odd =
    if Odd n then insert n ((Finset.range n).filter Odd) else (Finset.range n).filter Odd := by
  rw [Finset.range_add_one, Finset.filter_insert]

lemma oddRecip_succ (n : ℕ) : oddRecip (n + 1) = oddRecip n + if Odd n then (1 : ℝ) / n else 0 := by
  unfold oddRecip
  rw [filter_odd_succ]
  split_ifs with h
  · rw [Finset.sum_insert (by simp)]; ring
  · simp

lemma card_odd_succ (n : ℕ) : (((Finset.range (n + 1)).filter Odd).card : ℝ) =
    ((Finset.range n).filter Odd).card + if Odd n then 1 else 0 := by
  rw [filter_odd_succ]
  split_ifs with h
  · rw [Finset.card_insert_of_notMem (by simp)]; push_cast; ring
  · simp

lemma oddRecip_two_mul (m : ℕ) :
    oddRecip (2 * m) = (harmonic (2 * m) : ℝ) - (harmonic m : ℝ) / 2 ∧
      oddRecip (2 * m + 1) = oddRecip (2 * m) := by
  induction m with
  | zero =>
    refine ⟨by simp [oddRecip], ?_⟩
    simp [oddRecip, Finset.filter_singleton]
  | succ m ih =>
    have h1 : oddRecip (2 * m + 1) = oddRecip (2 * m) := ih.2
    have hodd : Odd (2 * m + 1) := ⟨m, rfl⟩
    have heven : ¬ Odd (2 * (m + 1)) := by rw [Nat.not_odd_iff_even]; exact ⟨m + 1, by ring⟩
    have h2 : oddRecip (2 * (m + 1)) = oddRecip (2 * m) + 1 / (2 * m + 1) := by
      rw [show 2 * (m + 1) = (2 * m + 1) + 1 by ring, oddRecip_succ, if_pos hodd, h1]
      push_cast; ring
    have h3 : oddRecip (2 * (m + 1) + 1) = oddRecip (2 * (m + 1)) := by
      rw [oddRecip_succ, if_neg heven, add_zero]
    refine ⟨?_, h3⟩
    rw [h2, ih.1, show 2 * (m + 1) = 2 * m + 1 + 1 by ring, harmonic_succ, harmonic_succ, harmonic_succ]
    push_cast
    field_simp
    ring

lemma card_odd_two_mul (m : ℕ) :
    (((Finset.range (2 * m)).filter Odd).card : ℝ) = m ∧
      (((Finset.range (2 * m + 1)).filter Odd).card : ℝ) = m := by
  induction m with
  | zero => simp
  | succ m ih =>
    have hodd : Odd (2 * m + 1) := ⟨m, rfl⟩
    have heven : ¬ Odd (2 * (m + 1)) := by rw [Nat.not_odd_iff_even]; exact ⟨m + 1, by ring⟩
    have h2 : (((Finset.range (2 * (m + 1))).filter Odd).card : ℝ) = m + 1 := by
      rw [show 2 * (m + 1) = (2 * m + 1) + 1 by ring, card_odd_succ, if_pos hodd, ih.2]
    refine ⟨by rw [h2]; push_cast; ring, ?_⟩
    rw [card_odd_succ, if_neg heven, h2]; push_cast; ring

/-- `0 ≤ H_n - log n - γ ≤ 1/n`. -/
lemma harmonic_err (n : ℕ) (hn : 1 ≤ n) :
    0 ≤ (harmonic n : ℝ) - log n - eulerMascheroniConstant ∧
      (harmonic n : ℝ) - log n - eulerMascheroniConstant ≤ 1 / n := by
  have h1 := eulerMascheroniSeq_lt_eulerMascheroniConstant n
  have h2 := eulerMascheroniConstant_lt_eulerMascheroniSeq' n
  have hn0 : n ≠ 0 := by omega
  simp only [eulerMascheroniSeq, eulerMascheroniSeq', hn0, if_false] at h1 h2
  have hnpos : (0 : ℝ) < n := by exact_mod_cast (by omega : 0 < n)
  have hlog : log (n + 1) - log n ≤ 1 / n := by
    rw [← log_div (by positivity) hnpos.ne']
    have := log_le_sub_one_of_pos (show (0 : ℝ) < (n + 1) / n by positivity)
    rw [show ((n : ℝ) + 1) / n - 1 = 1 / n by field_simp; ring] at this
    exact this
  constructor <;> linarith

/-- `|I(n) - t(n)| ≤ 3`. -/
theorem Iclosed_sub_trend_le (n : ℕ) : |Iclosed n - trendT n| ≤ 3 := by
  have hγ0 := one_half_lt_eulerMascheroniConstant
  have hγ1 := eulerMascheroniConstant_lt_two_thirds
  have hπ : 0 < π := pi_pos
  have hlog4 : log (4 * π) = log 2 + log (2 * π) := by
    rw [show (4 : ℝ) * π = 2 * (2 * π) by ring, log_mul (by norm_num) (by positivity)]
  rcases Nat.even_or_odd' n with ⟨m, rfl | rfl⟩
  · -- n = 2m
    rcases Nat.eq_zero_or_pos m with h0 | hm
    · subst h0; simp [Iclosed, trendT]
    have hm1 : 1 ≤ m := hm
    obtain ⟨hO, _⟩ := oddRecip_two_mul m
    obtain ⟨hP, _⟩ := card_odd_two_mul m
    rw [Iclosed_eq, hO, hP]
    obtain ⟨e1, e2⟩ := harmonic_err (2 * m) (by omega)
    obtain ⟨f1, f2⟩ := harmonic_err m hm1
    have hmpos : (0 : ℝ) < m := by exact_mod_cast hm
    have hlog2m : log ((2 * m : ℕ) : ℝ) = log 2 + log m := by
      push_cast; rw [log_mul (by norm_num) hmpos.ne']
    unfold trendT
    rw [hlog2m] at e1 e2 ⊢
    push_cast at e1 e2 ⊢
    rw [hlog4]
    have hb1 : (2 * m : ℝ) * (1 / (2 * m)) = 1 := by field_simp
    have hb2 : (m : ℝ) * (1 / m) = 1 := by field_simp
    rw [abs_le]
    constructor <;> nlinarith [mul_le_mul_of_nonneg_left e2 (by positivity : (0:ℝ) ≤ 2 * m),
      mul_le_mul_of_nonneg_left f2 (by positivity : (0:ℝ) ≤ m),
      mul_nonneg (by positivity : (0:ℝ) ≤ 2 * m) e1, mul_nonneg hmpos.le f1]
  · -- n = 2m + 1
    rcases Nat.eq_zero_or_pos m with h0 | hm
    · subst h0
      simp only [Iclosed, trendT, mul_zero, zero_add, Finset.range_one]
      have : (Finset.filter Odd {0}) = (∅ : Finset ℕ) := by decide
      rw [this, Finset.sum_empty]
      push_cast
      rw [log_one, hlog4]
      have hl2 := log_two_lt_d9
      have hl2' := log_two_gt_d9
      rw [abs_le]; constructor <;> nlinarith
    have hm1 : 1 ≤ m := hm
    obtain ⟨hO, hO'⟩ := oddRecip_two_mul m
    obtain ⟨_, hP'⟩ := card_odd_two_mul m
    rw [Iclosed_eq, hO', hO, hP']
    obtain ⟨e1, e2⟩ := harmonic_err (2 * m) (by omega)
    obtain ⟨f1, f2⟩ := harmonic_err m hm1
    have hmpos : (0 : ℝ) < m := by exact_mod_cast hm
    have hlog2m : log ((2 * m : ℕ) : ℝ) = log 2 + log m := by
      push_cast; rw [log_mul (by norm_num) hmpos.ne']
    have hlog2m1 : log (((2 * m + 1 : ℕ)) : ℝ) = log 2 + log m + log (1 + 1 / (2 * m)) := by
      push_cast
      rw [show (2 * (m : ℝ) + 1) = 2 * m * (1 + 1 / (2 * m)) by field_simp,
        log_mul (by positivity) (by positivity), log_mul (by norm_num) hmpos.ne']
    have hl0 : 0 ≤ log (1 + 1 / (2 * m)) := log_nonneg (by have : 0 ≤ 1 / (2 * (m : ℝ)) := by positivity
                                                           linarith)
    have hl1 : log (1 + 1 / (2 * m)) ≤ 1 / (2 * m) := by
      have := log_le_sub_one_of_pos (show (0 : ℝ) < 1 + 1 / (2 * m) by positivity)
      linarith
    unfold trendT
    rw [hlog2m] at e1 e2
    rw [hlog2m1]
    push_cast at e1 e2 ⊢
    rw [hlog4]
    have hb1 : (2 * m + 1 : ℝ) * (1 / (2 * m)) ≤ 3 / 2 := by
      rw [show (2 * m + 1 : ℝ) * (1 / (2 * m)) = 1 + 1 / (2 * m) by field_simp]
      have hm1' : (1 : ℝ) ≤ m := by exact_mod_cast hm1
      have : 1 / (2 * (m : ℝ)) ≤ 1 / 2 := by
        rw [div_le_div_iff₀ (by positivity) (by norm_num)]; linarith
      linarith
    have hb2 : (2 * m + 1 : ℝ) * (1 / m) ≤ 3 := by
      rw [show (2 * m + 1 : ℝ) * (1 / m) = 2 + 1 / m by field_simp]
      have : 1 / (m : ℝ) ≤ 1 := by rw [div_le_one hmpos]; exact_mod_cast hm1
      linarith
    rw [abs_le]
    constructor <;> nlinarith [mul_le_mul_of_nonneg_left e2 (by positivity : (0:ℝ) ≤ 2 * m + 1),
      mul_le_mul_of_nonneg_left f2 (by positivity : (0:ℝ) ≤ 2 * m + 1),
      mul_nonneg (by positivity : (0:ℝ) ≤ 2 * m + 1) e1,
      mul_nonneg (by positivity : (0:ℝ) ≤ 2 * m + 1) f1,
      mul_le_mul_of_nonneg_left hl1 (by positivity : (0:ℝ) ≤ 2 * m + 1),
      mul_nonneg (by positivity : (0:ℝ) ≤ 2 * m + 1) hl0]

end LeanLi
