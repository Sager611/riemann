/-
E1Lean.Numeric — rational enclosures of E1.md §4: √2; sin, cos at 353553/10⁶ by three angle doublings from
Mathlib's Taylor bounds (sin_bound, cos_bound); cot(1/(2√2)) and c(1/2); e^{8.10242} ≥ c_hi; κ₁ ≥ 37379/442500;
h(1/2; κ₁) > 0.0055 and k*(1/2) < κ₁.
-/
import Riemann.Deep.SantiAudit.Lean.AlgebraC

namespace E1Lean

variable (F : Fns)

theorem sqrt2_bounds : (141421356237 / 100000000000 : ℝ) ≤ F.sqrt 2 ∧ F.sqrt 2 ≤ 141421356238 / 100000000000 := by
  obtain ⟨hs, hss, _⟩ := sqrt2_facts F
  refine ⟨le_sqrt_of_sq_le F (by norm_num), ?_⟩
  by_contra hc
  rw [not_le] at hc
  have := mul_self_lt_mul_self (by norm_num) hc
  rw [hss] at this
  norm_num at this

/-- Taylor enclosures at 0 < a ≤ 1 (Real.sin_bound, Real.cos_bound). -/
theorem base_bounds {a : ℝ} (ha : 0 < a) (ha1 : a ≤ 1) :
    (a - a ^ 3 / 6 - a ^ 5 / 100 ≤ F.sin a ∧ F.sin a ≤ a - a ^ 3 / 6 + a ^ 5 / 100) ∧
    (1 - a ^ 2 / 2 - a ^ 4 * (5 / 96) ≤ F.cos a ∧ F.cos a ≤ 1 - a ^ 2 / 2 + a ^ 4 * (5 / 96)) := by
  have habs : |a| ≤ 1 := by rw [abs_of_pos ha]; exact ha1
  have hs := F.sin_bound a habs
  have hc := F.cos_bound a habs
  rw [abs_of_pos ha] at hs hc
  rw [abs_le] at hs hc
  exact ⟨⟨by linarith [hs.1], by linarith [hs.2]⟩, ⟨by linarith [hc.1], by linarith [hc.2]⟩⟩

/-- One angle doubling of nonnegative two-sided bounds: sin 2a = 2 sin a cos a, cos 2a = 1 − 2 sin² a. -/
theorem dbl {a sl sh cl ch sl' sh' cl' ch' : ℝ}
    (hs : sl ≤ F.sin a ∧ F.sin a ≤ sh) (hc : cl ≤ F.cos a ∧ F.cos a ≤ ch) (h0 : 0 ≤ sl) (h1 : 0 ≤ cl)
    (e1 : sl' ≤ 2 * sl * cl) (e2 : 2 * sh * ch ≤ sh') (e3 : cl' ≤ 1 - 2 * sh ^ 2) (e4 : 1 - 2 * sl ^ 2 ≤ ch') :
    (sl' ≤ F.sin (2 * a) ∧ F.sin (2 * a) ≤ sh') ∧ (cl' ≤ F.cos (2 * a) ∧ F.cos (2 * a) ≤ ch') := by
  have hcos2 : F.cos (2 * a) = 1 - 2 * F.sin a ^ 2 := by
    rw [F.cos_two_mul]; have := F.sin_sq_add_cos_sq a; linarith
  have hsa : 0 ≤ F.sin a := h0.trans hs.1
  have hca : 0 ≤ F.cos a := h1.trans hc.1
  have hsh : 0 ≤ sh := hsa.trans hs.2
  rw [F.sin_two_mul, hcos2]
  refine ⟨⟨?_, ?_⟩, ⟨?_, ?_⟩⟩
  · nlinarith [mul_le_mul hs.1 hc.1 h1 hsa]
  · nlinarith [mul_le_mul hs.2 hc.2 hca hsh]
  · nlinarith [mul_le_mul hs.2 hs.2 hsa hsh]
  · nlinarith [mul_le_mul hs.1 hs.1 h0 hsa]

theorem exp_nat : ∀ n : ℕ, F.exp n = F.exp 1 ^ n
  | 0 => by simp [F.exp_zero]
  | n + 1 => by rw [Nat.cast_succ, F.exp_add, exp_nat n, pow_succ]

/-- sin, cos at 353553/10⁶ (≤ 1/(2√2)): Taylor at 353553/(8·10⁶), then three doublings. -/
theorem trig_xlo : (34623307866193455119 / 100000000000000000000 : ℝ) ≤ F.sin (353553 / 1000000) ∧ F.cos (353553 / 1000000) ≤ (46907426160882655471 / 50000000000000000000 : ℝ) := by
  have b0 := base_bounds F (a := 353553 / 8000000) (by norm_num) (by norm_num)
  have s0 : (4417973723752344383 / 100000000000000000000 : ℝ) ≤ F.sin (353553 / 8000000) ∧ F.sin (353553 / 8000000) ≤ (4417974060925269981 / 100000000000000000000 : ℝ) :=
    ⟨le_trans (by norm_num) b0.1.1, le_trans b0.1.2 (by norm_num)⟩
  have c0 : (99902324097647092277 / 100000000000000000000 : ℝ) ≤ F.cos (353553 / 8000000) ∧ F.cos (353553 / 8000000) ≤ (99902363833901345223 / 100000000000000000000 : ℝ) :=
    ⟨le_trans (by norm_num) b0.2.1, le_trans b0.2.2 (by norm_num)⟩
  have d1 := dbl F (sl' := (8827316856103909843 / 100000000000000000000 : ℝ)) (sh' := (4413660520432949501 / 50000000000000000000 : ℝ)) (cl' := (99609630103939829577 / 100000000000000000000 : ℝ)) (ch' := (99609630163524676877 / 100000000000000000000 : ℝ)) s0 c0
    (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num)
  rw [show (2 : ℝ) * (353553 / 8000000) = 353553 / 4000000 by norm_num] at d1
  have d2 := dbl F (sl' := (17585715336935670213 / 100000000000000000000 : ℝ)) (sh' := (17585723684307038003 / 100000000000000000000 : ℝ)) (cl' := (24610392016207430907 / 25000000000000000000 : ℝ)) (ch' := (98441569542438875701 / 100000000000000000000 : ℝ)) d1.1 d1.2
    (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num)
  rw [show (2 : ℝ) * (353553 / 4000000) = 353553 / 2000000 by norm_num] at d2
  have d3 := dbl F (sl' := (34623307866193455119 / 100000000000000000000 : ℝ)) (sh' := (692466496409130273 / 2000000000000000000 : ℝ)) (cl' := (93814846449984049943 / 100000000000000000000 : ℝ)) (ch' := (46907426160882655471 / 50000000000000000000 : ℝ)) d2.1 d2.2
    (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num)
  rw [show (2 : ℝ) * (353553 / 2000000) = 353553 / 1000000 by norm_num] at d3
  exact ⟨d3.1.1, d3.2.2⟩

/-- c(1/2) ≥ c_lo (E1.md §4: c(1/2) = −0.16596385…; c_lo = -0.165967085215). -/
theorem cW_half_ge : (-16596708521464472031 / 100000000000000000000 : ℝ) ≤ cW F (1 / 2) := by
  obtain ⟨hs, hss, hs1⟩ := sqrt2_facts F
  obtain ⟨hsl, hsh⟩ := sqrt2_bounds F
  have ht := trig_xlo F
  have hcot : cot F (353553 / 1000000) ≤ (1354793318482534033 / 500000000000000000 : ℝ) := by
    unfold cot
    calc F.cos (353553 / 1000000) / F.sin (353553 / 1000000) ≤ (46907426160882655471 / 50000000000000000000 : ℝ) / (34623307866193455119 / 100000000000000000000 : ℝ) :=
          div_le_div₀ (by norm_num) ht.2 (by norm_num) ht.1
      _ ≤ (1354793318482534033 / 500000000000000000 : ℝ) := by norm_num
  set x₀ := (1 / 2 : ℝ) / F.sqrt 2 with hx₀
  have hxlo : (353553 / 1000000 : ℝ) ≤ x₀ := by rw [hx₀, le_div_iff₀ hs]; linarith
  have hx1 : x₀ ≤ 1 := by rw [hx₀, div_le_iff₀ hs]; linarith
  have hcotx : cot F x₀ ≤ (1354793318482534033 / 500000000000000000 : ℝ) := by
    rcases hxlo.lt_or_eq with hlt | heq
    · have := cot_sub_gt F (by norm_num) hlt hx1; linarith
    · rw [← heq]; exact hcot
  unfold cW
  rw [← hx₀]
  have h1 : cot F x₀ / F.sqrt 2 ≤ (1354793318482534033 / 500000000000000000 : ℝ) / F.sqrt 2 := div_le_div_of_nonneg_right hcotx hs.le
  have h2 : (1354793318482534033 / 500000000000000000 : ℝ) / F.sqrt 2 ≤ (1354793318482534033 / 500000000000000000 : ℝ) / (141421356237 / 100000000000) :=
    div_le_div_of_nonneg_left (by norm_num) (by norm_num) hsl
  have h3 : (-16596708521464472031 / 100000000000000000000 : ℝ) ≤ 2 - 1 / 2 / 2 - (1354793318482534033 / 500000000000000000 : ℝ) / (141421356237 / 100000000000) := by norm_num
  linarith

end E1Lean
