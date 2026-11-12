import Riemann.Deep.GL2Ceil.Lean.Tail
import Mathlib.Analysis.SpecialFunctions.Trigonometric.ArctanDeriv
import Mathlib.Analysis.Real.Pi.Bounds
import Mathlib.Analysis.Calculus.Deriv.MeanValue

/-!
# GL2Ceil (gl2lean): the R ≥ 28 tail proved (DERIVATION §4.3 (iii))

For R ≥ 28: e^R ≥ 4096, gd(R) ≥ gd(28) = 2 arctan tanh 14 ≥ 1.5687 (tanh 14 ≥ 0.999, arctan is
1-Lipschitz, π > 3.1415), e^{κR} ≥ e^{28κ} R/28 ≥ 7.0934 R/28 (28κ > 1); then
e^{(1-κ)R} ≤ crude(R, G(R)/ξ₀) by elementary algebra. So TailLargeR is a theorem.
-/

open Real MeasureTheory Set

noncomputable section

namespace GL2Ceil

lemma exp28_ge : 4096 ≤ Real.exp 28 := by
  have h7 : 8 ≤ Real.exp 7 := by linarith [Real.add_one_le_exp 7]
  have e : Real.exp 28 = Real.exp 7 ^ 4 := by
    rw [← Real.exp_nat_mul]
    norm_num
  rw [e]
  calc (4096:ℝ) = 8 ^ 4 := by norm_num
    _ ≤ Real.exp 7 ^ 4 := by gcongr

lemma tanh14_ge : (999:ℝ) / 1000 ≤ Real.tanh 14 := by
  rw [Real.tanh_eq_sinh_div_cosh, Real.sinh_eq, Real.cosh_eq, le_div_iff₀ (by positivity)]
  have h28 := exp28_ge
  have e : Real.exp 14 * Real.exp 14 = Real.exp 28 := by
    rw [← Real.exp_add]
    norm_num
  have hinv : Real.exp (-14) * Real.exp 14 = 1 := by
    rw [← Real.exp_add]
    norm_num
  have hpos := Real.exp_pos 14
  have hneg := Real.exp_pos (-14)
  nlinarith [hinv, e, h28, hpos, hneg]

lemma arctan_999_ge : Real.pi / 4 - 1 / 1000 ≤ Real.arctan (999 / 1000) := by
  have hd : ∀ x, deriv Real.arctan x ≤ 1 := by
    intro x
    rw [Real.deriv_arctan]
    show 1 / (1 + x ^ 2) ≤ 1
    rw [div_le_one (by positivity)]
    nlinarith [sq_nonneg x]
  have h := image_sub_le_mul_sub_of_deriv_le Real.differentiable_arctan hd
    (x := (999 / 1000 : ℝ)) (y := 1) (by norm_num)
  rw [Real.arctan_one] at h
  linarith

lemma gd28_ge : 15687 / 10000 ≤ gd 28 := by
  rw [gd_eq]
  have h1 := tanh14_ge
  have h2 : Real.arctan (999 / 1000) ≤ Real.arctan (Real.tanh (28 / 2)) := by
    rw [show (28:ℝ) / 2 = 14 by norm_num]
    exact Real.arctan_strictMono.monotone h1
  have h3 := arctan_999_ge
  have hpi := Real.pi_gt_d4
  linarith

lemma gd_mono28 {R : ℝ} (hR : 28 ≤ R) : gd 28 ≤ gd R := by
  unfold gd
  have hc : Continuous fun t : ℝ => 1 / Real.cosh t :=
    continuous_const.div₀ Real.continuous_cosh (fun t => (Real.cosh_pos t).ne')
  rw [← intervalIntegral.integral_add_adjacent_intervals (b := 28) (hc.intervalIntegrable 0 28)
    (hc.intervalIntegrable 28 R)]
  have : 0 ≤ ∫ t in (28:ℝ)..R, 1 / Real.cosh t :=
    intervalIntegral.integral_nonneg hR (fun t _ => by have := Real.cosh_pos t; positivity)
  linarith

lemma expk_ge {R : ℝ} (hR : 28 ≤ R) : 70934 / 10000 / 28 * R ≤ Real.exp (kup * R) := by
  have hE1 := Real.exp_one_gt_d9
  have h2 : Real.exp 2 = Real.exp 1 ^ 2 := by
    rw [← Real.exp_nat_mul]
    norm_num
  have h2' : (2.7182818283:ℝ) ^ 2 ≤ Real.exp 2 := by
    rw [h2]
    gcongr
  have h04 : (96:ℝ) / 100 ≤ Real.exp (-(4 / 100)) := by linarith [Real.add_one_le_exp (-(4 / 100))]
  have h196 : Real.exp 2 * Real.exp (-(4 / 100)) ≤ Real.exp (28 * kup) := by
    rw [← Real.exp_add]
    apply Real.exp_le_exp.2
    norm_num [kup]
  have hlin : 1 + kup * (R - 28) ≤ Real.exp (kup * (R - 28)) := by
    linarith [Real.add_one_le_exp (kup * (R - 28))]
  have hsplit : Real.exp (kup * R) = Real.exp (28 * kup) * Real.exp (kup * (R - 28)) := by
    rw [← Real.exp_add]
    congr 1
    ring
  have hR28 : R / 28 ≤ 1 + kup * (R - 28) := by
    unfold kup
    nlinarith
  have hA : (70934:ℝ) / 10000 ≤ Real.exp (28 * kup) := by
    have : (70934:ℝ) / 10000 ≤ (2.7182818283:ℝ) ^ 2 * (96 / 100) := by norm_num
    nlinarith [mul_le_mul h2' h04 (by norm_num) (Real.exp_pos 2).le]
  rw [hsplit]
  have hB : R / 28 ≤ Real.exp (kup * (R - 28)) := hR28.trans hlin
  have hBpos : 0 ≤ R / 28 := by linarith
  calc 70934 / 10000 / 28 * R = (70934 / 10000) * (R / 28) := by ring
    _ ≤ Real.exp (28 * kup) * Real.exp (kup * (R - 28)) :=
        mul_le_mul hA hB hBpos (Real.exp_pos _).le

lemma alg_tail (x y R g : ℝ) (hx : 0 < x) (hy : 0 < y) (hR : 0 < R) (hg : 0 < g)
    (h : x * (25 / 64 * R) + 2 * g ≤ x * g * y) :
    x / y ≤ 1 / 2 + x * x * (x * g / (4 * R * (25 / 256)))
      / (x * x + 2 * (x * g / (4 * R * (25 / 256)))) := by
  have hD : 0 < 4 * R * (25 / 256) := by positivity
  have hKpos : 0 < x * g / (4 * R * (25 / 256)) := by positivity
  have hKD : x * g / (4 * R * (25 / 256)) * (25 / 64 * R) = x * g := by
    field_simp
    ring
  have hKDx : x * g / (4 * R * (25 / 256)) * (25 / 64 * R) * x = x * g * x := by rw [hKD]
  have h3 : g * (x * x + 2 * (x * g / (4 * R * (25 / 256))))
      ≤ g * (x * y * (x * g / (4 * R * (25 / 256)))) := by
    nlinarith [mul_le_mul_of_nonneg_left h hKpos.le, hKDx]
  have h4 : x * x + 2 * (x * g / (4 * R * (25 / 256))) ≤ x * y * (x * g / (4 * R * (25 / 256))) :=
    le_of_mul_le_mul_left h3 hg
  have h5 : x / y ≤ x * x * (x * g / (4 * R * (25 / 256)))
      / (x * x + 2 * (x * g / (4 * R * (25 / 256)))) := by
    rw [div_le_div_iff₀ hy (by positivity)]
    nlinarith [mul_le_mul_of_nonneg_left h4 hx.le]
  linarith

/-- §4.3 (iii): TailLargeR holds. -/
theorem tailLargeR : TailLargeR := by
  intro R hR
  have hR0 : 0 < R := by linarith
  have hxpos := Real.exp_pos R
  have hx28 : 4096 ≤ Real.exp R := le_trans exp28_ge (Real.exp_le_exp.2 hR)
  have hgd : 15687 / 10000 ≤ gd R := le_trans gd28_ge (gd_mono28 hR)
  have hG : G R = Real.exp R * gd R / (4 * R) := G_eq hR0.ne'
  have hy := expk_ge hR
  have hypos := Real.exp_pos (kup * R)
  have hK : Real.exp R * (15687 / 10000) / (4 * R * xi0) ≤ 1 / xi0 * G R := by
    rw [hG]
    have h1 : Real.exp R * (15687 / 10000) ≤ Real.exp R * gd R :=
      mul_le_mul_of_nonneg_left hgd hxpos.le
    have h2 : 1 / xi0 * (Real.exp R * gd R / (4 * R)) = Real.exp R * gd R / (4 * R * xi0) := by
      unfold xi0
      field_simp
    rw [h2]
    exact div_le_div_of_nonneg_right h1 (by unfold xi0; positivity)
  have hK0 : 0 ≤ Real.exp R * (15687 / 10000) / (4 * R * xi0) := by
    unfold xi0
    positivity
  refine le_trans ?_ (crude_mono le_rfl hK0 hK)
  have hT : target R = Real.exp R / Real.exp (kup * R) := by
    unfold target
    rw [← Real.exp_sub]
    congr 1
    ring
  have hE : Real.exp (2 * R) = Real.exp R * Real.exp R := by
    rw [← Real.exp_add]
    congr 1
    ring
  rw [hT]
  unfold crude
  rw [hE]
  unfold xi0
  apply alg_tail _ _ _ _ hxpos hypos hR0 (by norm_num)
  have hRx : 4096 * 28 ≤ Real.exp R * R := mul_le_mul hx28 hR (by norm_num) hxpos.le
  nlinarith [mul_le_mul_of_nonneg_left hy (by positivity : (0:ℝ) ≤ Real.exp R * (15687 / 10000)), hRx]

/-- Main theorem with all analytic tails proved: from Cover, CertCells (1276 JSON cells) and the
crude strip inequality on [19/32, 28] (ceil.py, Arb). -/
theorem ceiling_bernard_main
    (hcov : Cover jsonCells) (hcert : CertCells jsonCells)
    (hstrip : ∀ R : ℝ, 19 / 32 ≤ R → R ≤ 28 → target R ≤ crude R (Cl 3 / xi0 * G R))
    (P : Polynomial ℝ) (hP0 : P.eval 0 = 0) (hP1 : P.eval 1 = 1) (M : ℕ) (lam : ℕ → ℝ)
    {R ξ : ℝ} (hR0 : 0 < R) (hξ : 0 < ξ) (hξ0 : ξ ≤ 25 / 256) :
    kappaC (fun x => P.eval x) (QB M lam) R ξ ≤ 0.07005633 :=
  ceiling_bernard_final hcov hcert hstrip tailLargeR P hP0 hP1 M lam hR0 hξ hξ0

/-- Same for real C¹ (in particular real polynomial) P, Q of the class. -/
theorem ceiling_main (cells : List Cell) (hcov : Cover cells) (hcert : CertCells cells)
    (hstrip : ∀ R : ℝ, 19 / 32 ≤ R → R ≤ 28 → target R ≤ crude R (Cl 3 / xi0 * G R))
    {P Q : ℝ → ℝ} (hP : PAdm P) (hQ : QAdm Q) {R ξ : ℝ} (hR0 : 0 < R) (hξ : 0 < ξ)
    (hξ0 : ξ ≤ 25 / 256) : kappa P Q R ξ ≤ 0.07005633 :=
  ceiling_final cells hcov hcert hstrip tailLargeR hP hQ hR0 hξ hξ0

end GL2Ceil


