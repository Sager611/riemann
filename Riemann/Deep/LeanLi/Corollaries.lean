module

public import Riemann.Deep.LeanLi.Thm1
public import Riemann.Deep.LeanLi.Lemma0
public import Riemann.Deep.LeanLi.LR.Cor3
public import Riemann.Deep.LeanLi.Statement

/-!
# LeanLi: Corollaries 1, 2, 4 and Lemma 0 for `R(n) = λ_n - t(n)` (Claim parts 2, 5; parts of 1, 3)
-/

@[expose] public section

open Real Finset Filter

namespace LeanLi

lemma c₁_eq : c₁ = LR.c₁ := rfl

lemma c₁_num : |c₁ - 0.0558| < 0.00005 := by
  have ha : (0.77459 : ℝ) < √(3 / 5) := by rw [Real.lt_sqrt (by norm_num)]; norm_num
  have ha' : √(3 / 5) < (0.77460 : ℝ) := by rw [Real.sqrt_lt' (by norm_num)]; norm_num
  have hb : (1.41421 : ℝ) < √2 := by rw [Real.lt_sqrt (by norm_num)]; norm_num
  have hb' : √2 < (1.41422 : ℝ) := by rw [Real.sqrt_lt' (by norm_num)]; norm_num
  have hp := Real.pi_gt_d4
  have hp' := Real.pi_lt_d4
  have hb0 : 0 < √(2 : ℝ) := by positivity
  have hd1 : π * √2 < 3.1416 * 1.41422 := mul_lt_mul'' hp' hb' Real.pi_pos.le hb0.le
  have hd2 : 3.1415 * 1.41421 < π * √2 := mul_lt_mul'' hp hb (by norm_num) (by norm_num)
  unfold c₁
  rw [abs_lt]
  constructor
  · rw [lt_sub_iff_add_lt, lt_div_iff₀ (by positivity)]
    nlinarith
  · rw [sub_lt_iff_lt_add, div_lt_iff₀ (by positivity)]
    nlinarith

lemma abs_liTrend_le (n : ℕ) : |liTrend n| ≤ 5 * (n : ℝ) ^ 2 := by
  rcases Nat.eq_zero_or_pos n with h0 | hn
  · subst h0; simp [liTrend]
  have hn1 : (1 : ℝ) ≤ n := by exact_mod_cast hn
  have hl0 : 0 ≤ Real.log n := Real.log_nonneg hn1
  have hl1 : Real.log n ≤ n := log_le_self' (by linarith)
  have hp : Real.log (2 * π) ≤ 6 := by
    have := Real.log_le_sub_one_of_pos (show (0 : ℝ) < 2 * π by positivity)
    have := Real.pi_lt_d2; linarith
  have hp0 : 0 ≤ Real.log (2 * π) := Real.log_nonneg (by linarith [Real.pi_gt_three])
  have hγ0 := Real.one_half_lt_eulerMascheroniConstant
  have hγ1 := Real.eulerMascheroniConstant_lt_two_thirds
  unfold liTrend
  rw [abs_mul, abs_of_nonneg (by positivity : (0 : ℝ) ≤ (n : ℝ) / 2)]
  have h : |Real.log n - 1 - Real.log (2 * π) + Real.eulerMascheroniConstant| ≤ (n : ℝ) + 8 := by
    rw [abs_le]; constructor <;> linarith
  calc (n : ℝ) / 2 * |Real.log n - 1 - Real.log (2 * π) + Real.eulerMascheroniConstant|
      ≤ (n : ℝ) / 2 * ((n : ℝ) + 8) := mul_le_mul_of_nonneg_left h (by positivity)
    _ ≤ 5 * (n : ℝ) ^ 2 := by nlinarith

/-- **Lemma 0** for `R` (from `lam_exp_growth`). -/
theorem lemma0_R (L : LitCore) : ¬ RiemannHypothesis → LR.Lemma0 R := by
  intro hRH
  obtain ⟨c, hc, hf⟩ := lam_exp_growth L hRH
  refine ⟨c / 2, by positivity, ?_⟩
  have hev : ∀ᶠ n : ℕ in atTop, 80 / c ^ 3 ≤ (n : ℝ) := tendsto_natCast_atTop_atTop.eventually_ge_atTop _
  refine (hf.and_eventually (hev.and (eventually_ge_atTop 1))).mono fun n ⟨h1, h2, h3⟩ => ?_
  have hn1 : (1 : ℝ) ≤ n := by exact_mod_cast h3
  set x := c / 2 * n with hx
  have hx0 : 0 ≤ x := by positivity
  have hexp : exp (c * n) = exp x * exp x := by rw [← exp_add]; congr 1; rw [hx]; ring
  have hq := Real.quadratic_le_exp_of_nonneg hx0
  have ha := Real.add_one_le_exp x
  have hdiff : x ^ 3 / 2 ≤ exp (c * n) - exp x := by
    rw [hexp]
    have : exp x * exp x - exp x = exp x * (exp x - 1) := by ring
    rw [this]
    have h1 : x ^ 2 / 2 ≤ exp x := by nlinarith
    have h2 : x ≤ exp x - 1 := by linarith
    calc x ^ 3 / 2 = (x ^ 2 / 2) * x := by ring
      _ ≤ exp x * (exp x - 1) := mul_le_mul h1 h2 hx0 (exp_pos x).le
  have ht := abs_liTrend_le n
  have hx3 : 5 * (n : ℝ) ^ 2 ≤ x ^ 3 / 2 := by
    have hc3 : 80 ≤ c ^ 3 * n := by
      rw [div_le_iff₀ (by positivity)] at h2; linarith
    have e : x ^ 3 / 2 = c ^ 3 * n * (n : ℝ) ^ 2 / 16 := by rw [hx]; ring
    rw [e]
    have := mul_le_mul_of_nonneg_right hc3 (sq_nonneg (n : ℝ))
    nlinarith
  have hR : |lam n| - |liTrend n| ≤ |R n| := by
    unfold R
    have := abs_sub_abs_le_abs_sub (lam n) (liTrend n)
    linarith
  show exp (c / 2 * n) ≤ |R n|
  rw [← hx]
  linarith

lemma rpow_quarter_nlogn (n : ℕ) :
    ((n : ℝ) * Real.log n) ^ (1 / 4 : ℝ) = √(√((n : ℝ) * Real.log n)) := by
  rcases Nat.eq_zero_or_pos n with h0 | hn
  · subst h0; simp
  · have hn1 : (1 : ℝ) ≤ n := by exact_mod_cast hn
    exact LR.rpow_quarter (mul_nonneg (by linarith) (Real.log_nonneg hn1))

/-- **Claim part 2.** -/
theorem part2 (L : LitCore) : Part2 := by
  have hRq : ∀ n : ℕ, ∀ c : ℝ, c * √(√((n : ℝ) * Real.log n)) ≤ |R n| →
      c * ((n : ℝ) * Real.log n) ^ (1 / 4 : ℝ) ≤ |R n| := fun n c h => by
    rw [rpow_quarter_nlogn]; exact h
  refine ⟨?_, ?_, lemma0_R L⟩
  · obtain ⟨c, hc, hf⟩ := LR.omega_uncond R RiemannHypothesis (thm1_li L) (lemma0_R L)
    exact ⟨c, hc, hf.mono fun n h => hRq n c h⟩
  · intro hRH
    have h1 := LR.cor1_num R (thm1_li L hRH)
    have h1' : ∃ᶠ n : ℕ in atTop, (0.289 : ℝ) * ((n : ℝ) * Real.log n) ^ (1 / 4 : ℝ) ≤ |R n| :=
      h1.mono fun n h => hRq n _ h
    refine ⟨h1', ?_⟩
    refine Filter.le_limsup_of_frequently_le ?_ (Filter.isBoundedUnder_of ⟨⊤, fun _ => le_top⟩)
    refine (h1'.and_eventually (eventually_ge_atTop 3)).mono fun n ⟨h, hn⟩ => ?_
    have hn3 : (3 : ℝ) ≤ n := by exact_mod_cast hn
    have hl : 1 ≤ Real.log n := by
      rw [Real.le_log_iff_exp_le (by linarith)]; linarith [Real.exp_one_lt_d9]
    have hpos : 0 < ((n : ℝ) * Real.log n) ^ (1 / 4 : ℝ) :=
      Real.rpow_pos_of_pos (by positivity) _
    rw [EReal.coe_le_coe_iff, le_div_iff₀ hpos]
    exact h

/-- Maślanka's trend is within `O(1)` of `t(n)`. -/
lemma vorosTrend_bdd (L : LitCore) : ∃ C : ℝ, ∀ n, |vorosTrend n - liTrend n| ≤ C := by
  obtain ⟨C, n₀, hC⟩ := L.voros
  refine ⟨3 / 4 + |C| + ∑ n ∈ Finset.range n₀, |vorosTrend n - liTrend n|, fun n => ?_⟩
  have hs0 : 0 ≤ ∑ n ∈ Finset.range n₀, |vorosTrend n - liTrend n| :=
    Finset.sum_nonneg fun _ _ => abs_nonneg _
  by_cases hn : n < n₀
  · have := Finset.single_le_sum (f := fun n => |vorosTrend n - liTrend n|) (fun _ _ => abs_nonneg _)
      (Finset.mem_range.2 hn)
    have := abs_nonneg C
    linarith
  · push Not at hn
    have h := hC n hn
    have hn1 : (0 : ℝ) ≤ n := Nat.cast_nonneg n
    have hCn : |C / n| ≤ |C| := by
      rcases Nat.eq_zero_or_pos n with h0 | h0
      · subst h0; simp
      · rw [abs_div, abs_of_pos (by exact_mod_cast h0 : (0 : ℝ) < n)]
        exact div_le_self (abs_nonneg _) (by exact_mod_cast h0)
    have : |vorosTrend n - liTrend n| ≤ 3 / 4 + |C / n| := by
      have := abs_sub_abs_le_abs_sub (vorosTrend n - liTrend n) (3 / 4)
      have e : vorosTrend n - liTrend n - 3 / 4 = vorosTrend n - liTrend n - 3 / 4 := rfl
      rw [abs_of_pos (by norm_num : (0 : ℝ) < 3 / 4)] at this
      linarith [le_abs_self (C / n), h]
    linarith

/-- **Claim part 5.** -/
theorem part5 (L : LitCore) : Part5 := by
  have key : ∀ trend : ℕ → ℝ, (∃ C : ℝ, ∀ n, |trend n - liTrend n| ≤ C) →
      (∀ a : ℝ, ¬ ∀ᶠ n : ℕ in atTop, |lam n - trend n| ≤ a * Real.log n) ∧
      (∀ ε : ℝ, ε < 1 / 4 →
        ¬ (fun n : ℕ => lam n - trend n) =O[atTop] (fun n : ℕ => (n : ℝ) ^ ε)) := by
    intro trend ⟨C, hC⟩
    have hlt : ∃ C, ∀ n, |(lam n - trend n) - R n| ≤ C := by
      refine ⟨C, fun n => ?_⟩
      have e : lam n - trend n - R n = -(trend n - liTrend n) := by unfold R; ring
      rw [e, abs_neg]; exact hC n
    refine ⟨fun a => LR.cor4 R (fun n => lam n - trend n) RiemannHypothesis (thm1_li L) (lemma0_R L) hlt a,
      fun ε hε hO => ?_⟩
    obtain ⟨C', hC'⟩ := hO.bound
    apply LR.cor4_rpow R (fun n => lam n - trend n) RiemannHypothesis (thm1_li L) (lemma0_R L) hlt hε C'
    filter_upwards [hC'] with n hn
    rw [Real.norm_eq_abs, Real.norm_eq_abs, abs_of_nonneg (Real.rpow_nonneg (Nat.cast_nonneg n) _)] at hn
    exact hn
  obtain ⟨h1, h2⟩ := key vorosTrend (vorosTrend_bdd L)
  exact ⟨key, h1, h2⟩

lemma keiperY_eq {m : ℕ} (hm : 1 ≤ m) : keiperY m = R m / m := by
  have hm0 : (m : ℝ) ≠ 0 := by exact_mod_cast (by omega : m ≠ 0)
  unfold keiperY R liTrend
  field_simp
  ring

lemma summable_keiper_iff (δ : ℝ) :
    Summable (fun m : ℕ => (m : ℝ) ^ δ * keiperY m ^ 2) ↔
      Summable (fun m : ℕ => (m : ℝ) ^ δ * (R m / m) ^ 2) := by
  rw [← summable_nat_add_iff 1, ← summable_nat_add_iff 1 (f := fun m : ℕ => (m : ℝ) ^ δ * (R m / m) ^ 2)]
  have : (fun n : ℕ => ((n + 1 : ℕ) : ℝ) ^ δ * keiperY (n + 1) ^ 2) =
      fun n : ℕ => ((n + 1 : ℕ) : ℝ) ^ δ * (R (n + 1) / ((n + 1 : ℕ) : ℝ)) ^ 2 := by
    funext n; rw [keiperY_eq (by omega)]
  simp only [this]

/-- **Claim part 3, last clause.** -/
theorem part3_half (L : LitCore) :
    RiemannHypothesis → ¬ Summable (fun m : ℕ => (m : ℝ) ^ (1 / 2 : ℝ) * keiperY m ^ 2) := by
  intro hRH
  rw [summable_keiper_iff]
  exact LR.not_summable_of_half_le R (thm1_li L hRH) le_rfl

lemma rpow_three_halves (x : ℝ) (hx : 0 ≤ x) : x ^ (3 / 2 : ℝ) = x * √x := by
  rw [show (3 / 2 : ℝ) = 1 + 1 / 2 by norm_num, Real.rpow_add' hx (by norm_num), Real.rpow_one,
    Real.sqrt_eq_rpow]

/-- **Claim part 1, first clause** (Theorem 1). -/
theorem part1_lower (L : LitCore) (hRH : RiemannHypothesis) :
    ∀ ε > (0 : ℝ), ∀ᶠ N : ℕ in atTop,
      (c₁ - ε) * ((N : ℝ) ^ (3 / 2 : ℝ) * Real.log N ^ (1 / 2 : ℝ)) ≤ ∑ n ∈ Finset.Icc 1 N, R n ^ 2 := by
  intro ε hε
  filter_upwards [thm1_li L hRH ε hε, eventually_ge_atTop 1] with N h hN
  have hN1 : (1 : ℝ) ≤ N := by exact_mod_cast hN
  rw [rpow_three_halves _ (by linarith), ← Real.sqrt_eq_rpow, c₁_eq]
  convert h using 2

end LeanLi
