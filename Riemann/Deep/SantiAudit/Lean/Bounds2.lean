/-
E1Lean.Bounds2 — upper ends of the E1.md §4 enclosures: c(1/2) ≤ −0.1659636, κ₁ ≤ 0.0844882, h(1/2; κ₁) ≤ 0.005543,
κ₁ − k*(1/2) ∈ (0.0078, 0.00787) (from H4; Real.exp_one_lt_d9, Real.exp_bound_div_one_sub_of_interval).
-/
import Riemann.Deep.SantiAudit.Lean.Arch

namespace E1Lean

variable (F : Fns)

/-- cos, sin at 3535534/10⁷ (≥ 1/(2√2)): Taylor at 3535534/(8·10⁷), then three doublings. -/
theorem trig_xhi : (93814832600643526557 / 100000000000000000000 : ℝ) ≤ F.cos (3535534 / 10000000) ∧ F.sin (3535534 / 10000000) ≤ (3462336234640537867 / 10000000000000000000 : ℝ) ∧ 0 < F.sin (3535534 / 10000000) := by
  have b0 := base_bounds F (a := 3535534 / 80000000) (by norm_num) (by norm_num)
  have s0 : (176719148754743339 / 4000000000000000000 : ℝ) ≤ F.sin (3535534 / 80000000) ∧ F.sin (3535534 / 80000000) ≤ (2208989528021708209 / 50000000000000000000 : ℝ) :=
    ⟨le_trans (by norm_num) b0.1.1, le_trans b0.1.2 (by norm_num)⟩
  have c0 : (19980464775317285829 / 20000000000000000000 : ℝ) ≤ F.cos (3535534 / 80000000) ∧ F.cos (3535534 / 80000000) ≤ (19980472722604101671 / 20000000000000000000 : ℝ) :=
    ⟨le_trans (by norm_num) b0.2.1, le_trans b0.2.2 (by norm_num)⟩
  have d1 := dbl F (sl' := (8827326817045512201 / 100000000000000000000 : ℝ)) (sh' := (275854093807224061 / 3125000000000000000 : ℝ)) (cl' := (3112800913162726077 / 3125000000000000000 : ℝ)) (ch' := (49804814640396243099 / 50000000000000000000 : ℝ)) s0 c0
    (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num)
  rw [show (2 : ℝ) * (3535534 / 80000000) = 3535534 / 40000000 by norm_num] at d1
  have d2 := dbl F (sl' := (17585735025206457987 / 100000000000000000000 : ℝ)) (sh' := (8792871686312493781 / 50000000000000000000 : ℝ)) (cl' := (24610391136920551567 / 25000000000000000000 : ℝ)) (ch' := (6152597876581336433 / 6250000000000000000 : ℝ)) d1.1 d1.2
    (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num)
  rw [show (2 : ℝ) * (3535534 / 40000000) = 3535534 / 20000000 by norm_num] at d2
  have d3 := dbl F (sl' := (34623345392045946091 / 100000000000000000000 : ℝ)) (sh' := (3462336234640537867 / 10000000000000000000 : ℝ)) (cl' := (93814832600643526557 / 100000000000000000000 : ℝ)) (ch' := (9381483847246453637 / 10000000000000000000 : ℝ)) d2.1 d2.2
    (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num)
  rw [show (2 : ℝ) * (3535534 / 20000000) = 3535534 / 10000000 by norm_num] at d3
  exact ⟨d3.2.1, d3.1.2, lt_of_lt_of_le (by norm_num) d3.1.1⟩

/-- c(1/2) ≤ c_hi = -0.165963667646 (E1.md §4: c(1/2) = −0.16596385…). -/
theorem cW_half_le : cW F (1 / 2) ≤ (-16596366764648141727 / 100000000000000000000 : ℝ) := by
  obtain ⟨hs, hss, hs1⟩ := sqrt2_facts F
  obtain ⟨hsl, hsh⟩ := sqrt2_bounds F
  obtain ⟨hc, hsn, hsp⟩ := trig_xhi F
  have hcot : (270958180381298083559 / 100000000000000000000 : ℝ) ≤ cot F (3535534 / 10000000) := by
    unfold cot
    calc (270958180381298083559 / 100000000000000000000 : ℝ) ≤ (93814832600643526557 / 100000000000000000000 : ℝ) / (3462336234640537867 / 10000000000000000000 : ℝ) := by norm_num
      _ ≤ F.cos (3535534 / 10000000) / F.sin (3535534 / 10000000) := div_le_div₀ (by linarith) hc hsp hsn
  set x₀ := (1 / 2 : ℝ) / F.sqrt 2 with hx₀
  have hx0 : 0 < x₀ := div_pos (by norm_num) hs
  have hxhi : x₀ ≤ (3535534 / 10000000 : ℝ) := by rw [hx₀, div_le_iff₀ hs]; nlinarith
  have hcotx : (270958180381298083559 / 100000000000000000000 : ℝ) ≤ cot F x₀ := by
    rcases hxhi.lt_or_eq with hlt | heq
    · have := cot_sub_gt F hx0 hlt (by norm_num); linarith
    · rw [heq]; exact hcot
  unfold cW
  rw [← hx₀]
  have h1 : (270958180381298083559 / 100000000000000000000 : ℝ) / F.sqrt 2 ≤ cot F x₀ / F.sqrt 2 := div_le_div_of_nonneg_right hcotx hs.le
  have h2 : (270958180381298083559 / 100000000000000000000 : ℝ) / (141421356238 / 100000000000) ≤ (270958180381298083559 / 100000000000000000000 : ℝ) / F.sqrt 2 :=
    div_le_div_of_nonneg_left (by norm_num) hs hsh
  have h3 : 2 - 1 / 2 / 2 - (270958180381298083559 / 100000000000000000000 : ℝ) / (141421356238 / 100000000000) ≤ (-16596366764648141727 / 100000000000000000000 : ℝ) := by norm_num
  linarith

theorem exp_mul_nat (y : ℝ) : ∀ n : ℕ, F.exp (n * y) = F.exp y ^ n
  | 0 => by simp [F.exp_zero]
  | n + 1 => by rw [Nat.cast_succ, add_mul, one_mul, F.exp_add, exp_mul_nat y n, pow_succ]

/-- exp L ≤ c ⟹ L ≤ log c (from log_le_iff_le_exp, exp_add, quadratic_le_exp_of_nonneg). -/
theorem le_log_of_exp_le {c L : ℝ} (hc : 0 < c) (hle : F.exp L ≤ c) : L ≤ F.log c := by
  by_contra hlt
  rw [not_le] at hlt
  have h1 : c ≤ F.exp (F.log c) := (F.log_le_iff_le_exp c (F.log c) hc).1 le_rfl
  have h2 := F.quadratic_le_exp_of_nonneg (L - F.log c) (by linarith)
  have h3 : F.exp L = F.exp (F.log c) * F.exp (L - F.log c) := by rw [← F.exp_add]; congr 1; ring
  have h4 := F.exp_pos (F.log c)
  have h6 : 1 < F.exp (L - F.log c) := by nlinarith [sq_nonneg (L - F.log c)]
  have h5 : F.exp (F.log c) < F.exp L := by rw [h3]; nlinarith [mul_lt_mul_of_pos_left h6 h4]
  linarith

/-- e^{8.10228} ≤ c_lo = 3302.286923204427565638 (e < 2.7182818286, e^x ≤ 1/(1 − x)). -/
theorem exp_Llo_le : F.exp (8 + 2557 / 25000) ≤ (3302286923204427565638 / 1000000000000000000 : ℝ) := by
  have h8 : F.exp 8 = F.exp 1 ^ 8 := by exact_mod_cast exp_nat F 8
  have hd : F.exp (2557 / 25000) = F.exp (2557 / 1600000) ^ 64 := by
    rw [← exp_mul_nat F (2557 / 1600000) 64]; norm_num
  have hb := F.exp_bound_div_one_sub_of_interval (2557 / 1600000) (by norm_num) (by norm_num)
  have he1 := pow_le_pow_left₀ (F.exp_pos 1).le F.exp_one_lt_d9.le 8
  have he2 := pow_le_pow_left₀ (F.exp_pos _).le hb 64
  rw [F.exp_add, h8, hd]
  calc F.exp 1 ^ 8 * F.exp (2557 / 1600000) ^ 64
        ≤ (2.7182818286 : ℝ) ^ 8 * (1 / (1 - 2557 / 1600000)) ^ 64 :=
          mul_le_mul he1 he2 (pow_nonneg (F.exp_pos _).le 64) (by norm_num)
    _ ≤ 3302286923204427565638 / 1000000000000000000 := by norm_num

/-- For c ≥ c_lo: κ = 1 − (177/20)⁻¹ log c ≤ 6231/73750 = 0.08448813…. -/
theorem kappa_upper {c : ℝ} (hlo : (3302286923204427565638 / 1000000000000000000 : ℝ) ≤ c) :
    kappa F c (177 / 20) ≤ 6231 / 73750 := by
  have hc0 : 0 < c := lt_of_lt_of_le (by norm_num) hlo
  have hl := le_log_of_exp_le F hc0 ((exp_Llo_le F).trans hlo)
  unfold kappa
  have h1 : (177 / 20 : ℝ)⁻¹ * (8 + 2557 / 25000) ≤ (177 / 20 : ℝ)⁻¹ * F.log c :=
    mul_le_mul_of_nonneg_left hl (by norm_num)
  have h2 : (6231 / 73750 : ℝ) = 1 - (177 / 20 : ℝ)⁻¹ * (8 + 2557 / 25000) := by norm_num
  linarith

/-- h(1/2; κ) ≤ 0.005543 for κ ≤ 6231/73750. -/
theorem h_half_le {κ : ℝ} (hk : κ ≤ 6231 / 73750) : h F (1 / 2) κ ≤ 0.005543 := by
  have hc1 := cW_lt_one F (by norm_num : (0 : ℝ) < 1 / 2) (by norm_num)
  have m1 := hc_mono_k F hc1 hk (by norm_num)
  have m2 := hc_mono_c F (by norm_num : (6231 / 73750 : ℝ) < 1) (cW_half_le F) (by norm_num)
  have top : hc F (-16596366764648141727 / 100000000000000000000 : ℝ) (6231 / 73750) ≤ 0.005543 := by
    unfold hc
    have := le_sqrt_of_sq_le F (a := (1 - 6231 / 73750) * (9 - 6231 / 73750 - 8 * (-16596366764648141727 / 100000000000000000000 : ℝ)))
      (B := 3 + 6231 / 73750 - 4 * 0.005543) (by norm_num)
    linarith
  unfold h
  linarith

/-- **E1.md §4 constants**, from H4: c(1/2), κ₁, h(1/2; κ₁), κ₁ − k*(1/2) enclosed (Arb: −0.1659638504,
  0.0844779080, 0.0055353042, 0.0078543529). -/
theorem e1_enclosures (I : Inputs) (h4 : H4 I) :
    (-16596708521464472031 / 100000000000000000000 : ℝ) ≤ cW F (1 / 2) ∧ cW F (1 / 2) ≤ (-16596366764648141727 / 100000000000000000000 : ℝ) ∧
    (37379 / 442500 : ℝ) ≤ κ₁ F I ∧ κ₁ F I ≤ 6231 / 73750 ∧
    (0.0055 : ℝ) < h F (1 / 2) (κ₁ F I) ∧ h F (1 / 2) (κ₁ F I) ≤ 0.005543 ∧
    (0.0078 : ℝ) < κ₁ F I - kstar F (1 / 2) ∧ κ₁ F I - kstar F (1 / 2) ≤ 0.00787 := by
  obtain ⟨hk, hk1⟩ := kappa_bounds F h4.1 h4.2
  have hku := kappa_upper F h4.1
  have hlo := cW_half_ge F
  have hhi := cW_half_le F
  have ks1 : kstar F (1 / 2) ≤ 1 - 2 / (2 - (-16596708521464472031 / 100000000000000000000 : ℝ)) := by
    unfold kstar
    have := div_le_div_of_nonneg_left (by norm_num : (0 : ℝ) ≤ 2) (by linarith) (by linarith : 2 - cW F (1 / 2) ≤ 2 - (-16596708521464472031 / 100000000000000000000 : ℝ))
    linarith
  have ks2 : 1 - 2 / (2 - (-16596366764648141727 / 100000000000000000000 : ℝ)) ≤ kstar F (1 / 2) := by
    unfold kstar
    have := div_le_div_of_nonneg_left (by norm_num : (0 : ℝ) ≤ 2) (by norm_num) (by linarith : 2 - (-16596366764648141727 / 100000000000000000000 : ℝ) ≤ 2 - cW F (1 / 2))
    linarith
  have g1 : (0.0078 : ℝ) < 37379 / 442500 - (1 - 2 / (2 - (-16596708521464472031 / 100000000000000000000 : ℝ))) := by norm_num
  have g2 : (6231 / 73750 : ℝ) - (1 - 2 / (2 - (-16596366764648141727 / 100000000000000000000 : ℝ))) ≤ 0.00787 := by norm_num
  unfold κ₁
  exact ⟨hlo, hhi, hk, hku, h_half_gt F hk hk1, h_half_le F hku, by linarith, by linarith⟩

end E1Lean
