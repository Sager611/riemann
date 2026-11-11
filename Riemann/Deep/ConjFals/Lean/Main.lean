import Riemann.Deep.ConjFals.Lean.Bounds
import Mathlib.Analysis.SpecialFunctions.Pow.Asymptotics
import Mathlib.Order.LiminfLimsup

/-!
# lim_{a→∞} s_a = 1 (conjecture of Johnston–Leong–Tudzi, arXiv 2408.04143)

`W a x = Σ_{1 ≤ n ≤ x} (−a)^{Ω(n)}`, `s a = limsup_{x→∞} |W a x| / x^{log₂ a}`.
Main theorem `tendsto_s`: `s a → 1` as `a → ∞`. For `a ≥ 64` the limsup is squeezed:
`L a ≤ s a ≤ U a` (`le_s`, `s_le`) with `U a = (a/(a−1))(1 + 9a^{−1/3})` and
`L a = a/(a+1) − (a/(a−1)) 9a^{−1/3}`, both tending to `1`. `explicit_bounds` is the
limsup-free form: for every `ε > 0` and all large `a`, `|W a x| ≤ (1+ε) x^α` for all `x ≥ 1`
and `|W a 2^k| ≥ (1−ε) (2^k)^α` for all large `k`.
-/

open Filter Topology Finset Real
open scoped ArithmeticFunction.Omega

namespace OmegaLean

/-- `W a x = Σ_{1 ≤ n ≤ x} (−a)^{Ω(n)}`. -/
noncomputable def W (a x : ℝ) : ℝ := ∑ n ∈ Icc 1 ⌊x⌋₊, (-a) ^ (ArithmeticFunction.cardFactors n)

/-- `s a = limsup_{x→∞} |W a x| / x^{log₂ a}`. -/
noncomputable def s (a : ℝ) : ℝ := limsup (fun x : ℝ => |W a x| / x ^ logb 2 a) atTop

lemma W_eq_WN (a x : ℝ) : W a x = WN a ⌊x⌋₊ := rfl

/-- Upper constant `U a = (a/(a−1))(1 + 9a^{−1/3})`. -/
noncomputable def U (a : ℝ) : ℝ := (a / (a - 1)) * (1 + 9 * a ^ (-(1/3 : ℝ)))

/-- Lower constant `L a = a/(a+1) − (a/(a−1)) 9a^{−1/3}`. -/
noncomputable def L (a : ℝ) : ℝ := a / (a + 1) - (a / (a - 1)) * (9 * a ^ (-(1/3 : ℝ)))

lemma logb_ge_six {a : ℝ} (ha : 64 ≤ a) : 6 ≤ logb 2 a := by
  have h64 : logb 2 64 = 6 := by
    rw [show (64:ℝ) = 2 ^ (6:ℕ) by norm_num, Real.logb_pow, Real.logb_self_eq_one (by norm_num)]
    norm_num
  rw [← h64]
  unfold Real.logb
  exact div_le_div_of_nonneg_right (Real.log_le_log (by norm_num) ha)
    (Real.log_pos (by norm_num)).le

/-- `|W a x| / x^α ≤ U a` for `a ≥ 64`, `x ≥ 1`. -/
lemma ratio_le {a : ℝ} (ha : 64 ≤ a) {x : ℝ} (hx : 1 ≤ x) :
    |W a x| / x ^ logb 2 a ≤ U a := by
  have ha1 : 1 < a := by linarith
  have hα := logb_ge_six ha
  have hx0 : 0 < x := by linarith
  have hxα : 0 < x ^ logb 2 a := Real.rpow_pos_of_pos hx0 _
  rw [div_le_iff₀ hxα, W_eq_WN]
  have hN : ((⌊x⌋₊ : ℕ) : ℝ) ^ logb 2 a ≤ x ^ logb 2 a :=
    Real.rpow_le_rpow (Nat.cast_nonneg _) (Nat.floor_le hx0.le) (by linarith)
  have hU : 0 ≤ U a := by
    unfold U
    apply mul_nonneg (div_nonneg (by linarith) (by linarith))
    have : 0 ≤ a ^ (-(1/3:ℝ)) := (Real.rpow_pos_of_pos (by linarith) _).le
    linarith
  calc |WN a ⌊x⌋₊| ≤ U a * ((⌊x⌋₊ : ℕ) : ℝ) ^ logb 2 a := abs_WN_le ha1 hα _
    _ ≤ U a * x ^ logb 2 a := mul_le_mul_of_nonneg_left hN hU

/-- `|W a 2^k| ≥ L a · a^k − 1` for `a ≥ 64`. -/
lemma W_two_pow_ge {a : ℝ} (ha : 64 ≤ a) (k : ℕ) :
    L a * a ^ k - 1 ≤ |W a ((2:ℝ) ^ k)| := by
  have ha1 : 1 < a := by linarith
  have hα := logb_ge_six ha
  have hfloor : ⌊(2:ℝ) ^ k⌋₊ = 2 ^ k := by
    rw [(by norm_num : ((2:ℝ) ^ k) = ((2 ^ k : ℕ) : ℝ)), Nat.floor_natCast]
  rw [W_eq_WN, hfloor]
  have hlow := abs_WN_two_pow_ge ha1 hα k
  have e1 : (a ^ (k + 1) - 1) / (a + 1) = a / (a + 1) * a ^ k - 1 / (a + 1) := by ring
  have h1a : 1 / (a + 1) ≤ 1 := by
    rw [div_le_one (by linarith)]; linarith
  rw [e1] at hlow
  unfold L
  linarith [hlow, h1a]

lemma ratio_nonneg (a x : ℝ) : 0 ≤ |W a x| / x ^ logb 2 a := by
  rcases lt_or_ge x 1 with hx | hx
  · have : W a x = 0 := by
      rw [W_eq_WN, Nat.floor_eq_zero.mpr hx, WN_zero]
    rw [this, abs_zero, zero_div]
  · exact div_nonneg (abs_nonneg _) (Real.rpow_nonneg (by linarith) _)

/-- Along `x = 2^k` the ratio is frequently at least `L a − ε`. -/
lemma frequently_ratio_ge {a : ℝ} (ha : 64 ≤ a) {ε : ℝ} (hε : 0 < ε) :
    ∃ᶠ x in atTop, L a - ε ≤ |W a x| / x ^ logb 2 a := by
  rw [frequently_atTop]
  intro x₀
  obtain ⟨k, hk⟩ := pow_unbounded_of_one_lt (max x₀ ε⁻¹) (by norm_num : (1:ℝ) < 2)
  refine ⟨(2:ℝ) ^ k, le_trans (le_max_left _ _) hk.le, ?_⟩
  have ha0 : 0 < a := by linarith
  have hak : 0 < a ^ k := pow_pos ha0 k
  have h2a : (2:ℝ) ^ k ≤ a ^ k := pow_le_pow_left₀ (by norm_num) (by linarith) k
  rw [two_pow_rpow_logb ha0 k, le_div_iff₀ hak]
  have hεk : 1 ≤ ε * a ^ k := by
    have h1 : ε⁻¹ < a ^ k := lt_of_le_of_lt (le_max_right _ _) (hk.trans_le h2a)
    have h2 := mul_lt_mul_of_pos_left h1 hε
    rw [mul_inv_cancel₀ hε.ne'] at h2
    exact h2.le
  have := W_two_pow_ge ha k
  linarith [this, hεk]

/-- `s a ≤ U a` for `a ≥ 64`. -/
theorem s_le {a : ℝ} (ha : 64 ≤ a) : s a ≤ U a := by
  apply limsup_le_of_le
  · exact isCoboundedUnder_le_of_le atTop (ratio_nonneg a)
  · filter_upwards [eventually_ge_atTop (1:ℝ)] with x hx
    exact ratio_le ha hx

/-- `L a ≤ s a` for `a ≥ 64` (the limsup is a genuine one: the ratio is bounded). -/
theorem le_s {a : ℝ} (ha : 64 ≤ a) : L a ≤ s a := by
  have hb : IsBoundedUnder (· ≤ ·) atTop (fun x : ℝ => |W a x| / x ^ logb 2 a) :=
    isBoundedUnder_of_eventually_le (a := U a)
      (by filter_upwards [eventually_ge_atTop (1:ℝ)] with x hx; exact ratio_le ha hx)
  apply le_of_forall_pos_lt_add
  intro ε hε
  have h := le_limsup_of_frequently_le (frequently_ratio_ge ha (half_pos hε)) hb
  unfold s
  linarith

lemma tendsto_div_sub_one : Tendsto (fun a : ℝ => a / (a - 1)) atTop (𝓝 1) := by
  have h1 : Tendsto (fun a : ℝ => a - 1) atTop atTop :=
    tendsto_atTop_atTop.mpr fun b => ⟨b + 1, fun a ha => by linarith⟩
  have h2 : Tendsto (fun a : ℝ => 1 + (a - 1)⁻¹) atTop (𝓝 (1 + 0)) :=
    tendsto_const_nhds.add (tendsto_inv_atTop_zero.comp h1)
  rw [add_zero] at h2
  apply h2.congr'
  filter_upwards [eventually_gt_atTop (1:ℝ)] with a ha
  have h : a - 1 ≠ 0 := by linarith
  rw [eq_div_iff h, add_mul, inv_mul_cancel₀ h]
  ring

lemma tendsto_div_add_one : Tendsto (fun a : ℝ => a / (a + 1)) atTop (𝓝 1) := by
  have h1 : Tendsto (fun a : ℝ => a + 1) atTop atTop :=
    tendsto_atTop_atTop.mpr fun b => ⟨b, fun a ha => by linarith⟩
  have h2 : Tendsto (fun a : ℝ => 1 - (a + 1)⁻¹) atTop (𝓝 (1 - 0)) :=
    tendsto_const_nhds.sub (tendsto_inv_atTop_zero.comp h1)
  rw [sub_zero] at h2
  apply h2.congr'
  filter_upwards [eventually_gt_atTop (0:ℝ)] with a ha
  have h : a + 1 ≠ 0 := by linarith
  rw [eq_div_iff h, sub_mul, inv_mul_cancel₀ h]
  ring

lemma tendsto_U : Tendsto U atTop (𝓝 1) := by
  have h3 : Tendsto (fun a : ℝ => a ^ (-(1/3 : ℝ))) atTop (𝓝 0) :=
    tendsto_rpow_neg_atTop (by norm_num)
  have h4 : Tendsto (fun a : ℝ => 1 + 9 * a ^ (-(1/3 : ℝ))) atTop (𝓝 (1 + 9 * 0)) :=
    tendsto_const_nhds.add (tendsto_const_nhds.mul h3)
  have h5 := tendsto_div_sub_one.mul h4
  rw [mul_zero, add_zero, mul_one] at h5
  exact h5

lemma tendsto_L : Tendsto L atTop (𝓝 1) := by
  have h3 : Tendsto (fun a : ℝ => a ^ (-(1/3 : ℝ))) atTop (𝓝 0) :=
    tendsto_rpow_neg_atTop (by norm_num)
  have h4 : Tendsto (fun a : ℝ => a / (a - 1) * (9 * a ^ (-(1/3 : ℝ)))) atTop (𝓝 (1 * (9 * 0))) :=
    tendsto_div_sub_one.mul (tendsto_const_nhds.mul h3)
  have h5 := tendsto_div_add_one.sub h4
  rw [mul_zero, mul_zero, sub_zero] at h5
  exact h5

/-- **Main theorem** (Johnston–Leong–Tudzi's conjecture): `lim_{a→∞} s_a = 1`. -/
theorem tendsto_s : Tendsto s atTop (𝓝 1) := by
  apply tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_L tendsto_U
  · filter_upwards [eventually_ge_atTop (64:ℝ)] with a ha using le_s ha
  · filter_upwards [eventually_ge_atTop (64:ℝ)] with a ha using s_le ha

/-- Limsup-free form: for every `ε > 0` there is `a₀` such that for all `a ≥ a₀`:
`|W a x| ≤ (1+ε) x^α` for every `x ≥ 1`, and `|W a 2^k| ≥ (1−ε) (2^k)^α` for all large `k`. -/
theorem explicit_bounds {ε : ℝ} (hε : 0 < ε) : ∃ a₀ : ℝ, ∀ a ≥ a₀,
    (∀ x ≥ (1:ℝ), |W a x| ≤ (1 + ε) * x ^ logb 2 a) ∧
    (∀ᶠ k : ℕ in atTop, (1 - ε) * ((2:ℝ) ^ k) ^ logb 2 a ≤ |W a ((2:ℝ) ^ k)|) := by
  have hU := (tendsto_order.1 tendsto_U).2 (1 + ε) (by linarith)
  have hL := (tendsto_order.1 tendsto_L).1 (1 - ε / 2) (by linarith)
  obtain ⟨a₀, ha₀⟩ := eventually_atTop.mp (hU.and (hL.and (eventually_ge_atTop (64:ℝ))))
  refine ⟨a₀, fun a ha => ?_⟩
  obtain ⟨hUa, hLa, h64⟩ := ha₀ a ha
  have ha0 : 0 < a := by linarith
  refine ⟨fun x hx => ?_, ?_⟩
  · have hxα : 0 < x ^ logb 2 a := Real.rpow_pos_of_pos (by linarith) _
    have h := ratio_le h64 hx
    rw [div_le_iff₀ hxα] at h
    linarith [mul_le_mul_of_nonneg_right hUa.le hxα.le]
  · have hk : ∀ᶠ k : ℕ in atTop, 2 / ε ≤ a ^ k :=
      (tendsto_pow_atTop_atTop_of_one_lt (by linarith : (1:ℝ) < a)).eventually_ge_atTop (2 / ε)
    filter_upwards [hk] with k hk
    rw [two_pow_rpow_logb ha0 k]
    have h := W_two_pow_ge h64 k
    have h2 : 2 ≤ ε * a ^ k := by
      rw [div_le_iff₀ hε] at hk; linarith
    linarith [mul_le_mul_of_nonneg_right hLa.le (pow_pos ha0 k).le]

end OmegaLean
