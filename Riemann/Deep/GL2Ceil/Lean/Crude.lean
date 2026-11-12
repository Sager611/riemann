import Riemann.Deep.GL2Ceil.Lean.InnerMin

/-!
# GL2Ceil (gl2lean): the derivative-free crude bound (DERIVATION §4.3)

Drop γ∫w'² ≥ 0; fold y ↦ 1-y and minimise e^{2Ry}q² + e^{2R(1-y)}(β-q)² pointwise:
W₀ ≥ β² G(R), G(R) = ∫₀^{1/2} dy / (e^{2R(y-1)} + e^{-2Ry}) (= e^R gd(R)/(4R), G_eq); then minimise
over β: J ≥ crude(R, αG(R)), crude(R, K) = 1/2 + e^{2R}K/(e^{2R} + 2K). crude increases in R and K,
G increases in R, C_λ increases in λ: a crude cell [R_a, R_b] × [λ_a, λ_b] is closed by the
single arithmetic inequality e^{(1-κ_up)R_b} ≤ crude(R_a, (C_{λ_a}/ξ) G(R_a)).
-/

open Real MeasureTheory Set

noncomputable section

namespace GL2Ceil

/-- G(R) = ∫₀^{1/2} 1/(e^{2R(y-1)} + e^{-2Ry}) dy. -/
def G (R : ℝ) : ℝ :=
  ∫ y in (0:ℝ)..(1 / 2), 1 / (Real.exp (2 * R * (y - 1)) + Real.exp (-(2 * R * y)))

/-- crude(R, K) = 1/2 + e^{2R}K/(e^{2R} + 2K). -/
def crude (R K : ℝ) : ℝ := 1 / 2 + Real.exp (2 * R) * K / (Real.exp (2 * R) + 2 * K)

lemma G_integrand_cont (R : ℝ) :
    Continuous fun y : ℝ => 1 / (Real.exp (2 * R * (y - 1)) + Real.exp (-(2 * R * y))) :=
  continuous_const.div₀ (by fun_prop) (fun y => by positivity)

lemma G_nonneg (R : ℝ) : 0 ≤ G R :=
  intervalIntegral.integral_nonneg (by norm_num) (fun y _ => by positivity)

lemma G_mono {R1 R2 : ℝ} (h12 : R1 ≤ R2) : G R1 ≤ G R2 := by
  unfold G
  apply intervalIntegral.integral_mono_on (by norm_num) ((G_integrand_cont R1).intervalIntegrable _ _)
    ((G_integrand_cont R2).intervalIntegrable _ _)
  intro y hy
  obtain ⟨hy0, hy1⟩ := hy
  apply one_div_le_one_div_of_le (by positivity)
  apply add_le_add
  · exact Real.exp_le_exp.2 (by nlinarith)
  · exact Real.exp_le_exp.2 (by nlinarith)

lemma crude_mono {R1 R2 K1 K2 : ℝ} (hR : R1 ≤ R2) (hK1 : 0 ≤ K1) (hK : K1 ≤ K2) :
    crude R1 K1 ≤ crude R2 K2 := by
  unfold crude
  have he1 := Real.exp_pos (2 * R1)
  have he2 := Real.exp_pos (2 * R2)
  have he : Real.exp (2 * R1) ≤ Real.exp (2 * R2) := Real.exp_le_exp.2 (by linarith)
  have hK2 : 0 ≤ K2 := hK1.trans hK
  have : Real.exp (2 * R1) * K1 / (Real.exp (2 * R1) + 2 * K1)
      ≤ Real.exp (2 * R2) * K2 / (Real.exp (2 * R2) + 2 * K2) := by
    rw [div_le_div_iff₀ (by positivity) (by positivity)]
    nlinarith [mul_nonneg (mul_nonneg he1.le he2.le) (sub_nonneg.2 hK),
      mul_nonneg (mul_nonneg hK1 hK2) (sub_nonneg.2 he)]
  linarith

/-- Pointwise minimum: β² A²B²/(A²+B²) ≤ A²q² + B²(β-q)². -/
lemma pointwise_min (A B q β : ℝ) (hA : 0 < A) (hB : 0 < B) :
    β ^ 2 * (1 / (1 / B ^ 2 + 1 / A ^ 2)) ≤ A ^ 2 * q ^ 2 + B ^ 2 * (β - q) ^ 2 := by
  have e : 1 / (1 / B ^ 2 + 1 / A ^ 2) = A ^ 2 * B ^ 2 / (A ^ 2 + B ^ 2) := by
    field_simp
  rw [e, mul_div_assoc', div_le_iff₀ (by positivity)]
  nlinarith [sq_nonneg (A ^ 2 * q - B ^ 2 * (β - q))]

/-- W₀ ≥ β² G(R). -/
lemma W0_ge_G {R : ℝ} {Q : ℝ → ℝ} (hQc : Continuous Q) {β : ℝ}
    (hβ : ∀ y ∈ Icc (0:ℝ) 1, Q y + Q (1 - y) = β) : β ^ 2 * G R ≤ W0 R Q := by
  have hwc := w_cont (R := R) hQc
  have hw2 : Continuous fun y => w R Q y ^ 2 := hwc.fun_pow 2
  have hw2' : Continuous fun y => w R Q (1 - y) ^ 2 := by
    have : Continuous fun y : ℝ => w R Q (1 - y) := hwc.comp (by fun_prop)
    exact this.fun_pow 2
  have split : W0 R Q = ∫ y in (0:ℝ)..(1 / 2), (w R Q y ^ 2 + w R Q (1 - y) ^ 2) := by
    unfold W0
    rw [← intervalIntegral.integral_add_adjacent_intervals (b := 1 / 2)
      (hw2.intervalIntegrable _ _) (hw2.intervalIntegrable _ _),
      intervalIntegral.integral_add (hw2.intervalIntegrable _ _) (hw2'.intervalIntegrable _ _)]
    congr 1
    have h := intervalIntegral.integral_comp_sub_left (fun y => w R Q y ^ 2) (a := 0) (b := 1 / 2) 1
    norm_num at h
    rw [h]
  rw [split]
  unfold G
  rw [← intervalIntegral.integral_const_mul]
  apply intervalIntegral.integral_mono_on (by norm_num)
    (((G_integrand_cont R).intervalIntegrable _ _).const_mul _)
    ((hw2.fun_add hw2').intervalIntegrable _ _)
  intro y hy
  obtain ⟨hy0, hy1⟩ := hy
  have hq : Q (1 - y) = β - Q y := by
    have := hβ y ⟨hy0, by linarith⟩
    linarith
  set A := Real.exp (R * y) with hA
  set B := Real.exp (R * (1 - y)) with hB
  have hApos : 0 < A := Real.exp_pos _
  have hBpos : 0 < B := Real.exp_pos _
  have hu : Real.exp (2 * R * (y - 1)) = 1 / B ^ 2 := by
    rw [hB, eq_div_iff (by positivity), ← Real.exp_nat_mul, ← Real.exp_add]
    norm_num
    ring_nf
  have hv : Real.exp (-(2 * R * y)) = 1 / A ^ 2 := by
    rw [hA, eq_div_iff (by positivity), ← Real.exp_nat_mul, ← Real.exp_add]
    norm_num
    ring_nf
  have hw1 : w R Q y ^ 2 = A ^ 2 * Q y ^ 2 := by simp only [w, hA]; ring
  have hw2y : w R Q (1 - y) ^ 2 = B ^ 2 * (β - Q y) ^ 2 := by simp only [w, hB, hq]; ring
  rw [hu, hv, hw1, hw2y]
  exact pointwise_min A B (Q y) β hApos hBpos

/-- §4.3: J ≥ crude(R, α G(R)) for every Q in the class. -/
theorem crude_le_Jag {R α γ : ℝ} (hα : 0 ≤ α) (hγ : 0 ≤ γ) {Q : ℝ → ℝ} (hQ : QAdm Q) :
    crude R (α * G R) ≤ Jag R α γ Q := by
  obtain ⟨hQd, _⟩ := contDiff_one_iff_deriv.mp hQ.smooth
  obtain ⟨β, hβ⟩ := hQ.symm
  have hW := W0_ge_G (R := R) hQd.continuous hβ
  have hQ1 : Q 1 = β - 1 := by
    have := hβ 0 ⟨le_refl 0, zero_le_one⟩
    rw [sub_zero, hQ.zero] at this
    linarith
  have hw1 : w R Q 1 ^ 2 = Real.exp (2 * R) * (β - 1) ^ 2 := by
    simp only [w, mul_one, hQ1]
    rw [show (2:ℝ) * R = R + R by ring, Real.exp_add]
    ring
  have hG := G_nonneg R
  have hK : 0 ≤ α * G R := mul_nonneg hα hG
  have hE := Real.exp_pos (2 * R)
  have h1 : α * (β ^ 2 * G R) ≤ α * W0 R Q := mul_le_mul_of_nonneg_left hW hα
  have h2 : 0 ≤ γ * W1 R Q := mul_nonneg hγ (W1_nonneg R Q)
  have h3 : Real.exp (2 * R) * (α * G R) / (Real.exp (2 * R) + 2 * (α * G R))
      ≤ Real.exp (2 * R) * (β - 1) ^ 2 / 2 + (α * G R) * β ^ 2 := by
    rw [div_le_iff₀ (by positivity)]
    nlinarith [sq_nonneg (Real.exp (2 * R) * (β - 1) / 2 + (α * G R) * β), hK, hE]
  unfold crude Jag
  rw [hw1]
  nlinarith [h1, h2, h3]

end GL2Ceil


