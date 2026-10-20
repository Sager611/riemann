import Riemann.Deep.QRHCor.LOne
import Riemann.Deep.PNTError.InvZeta

/-!
# Lindelöf for Dirichlet L-functions right of 7/8, uniformly in q and t (rung 1b)

For η, ε > 0 there is C such that for every q, every χ ≠ 1 mod q and every s with
Re s ≥ 7/8 + η:
  ‖L(s, χ)‖ ≤ C (q(|t| + 8))^ε  and  ‖L(s, χ)⁻¹‖ ≤ C (q(|t| + 8))^ε   (`LFunction_lindelof`).
This is the Lindelöf hypothesis in the q- and t-aspects for all Dirichlet L-functions to the right
of 7/8, from OpenAI's zero-free half-plane (Riemann/Deep/QRH/).

Proof: the three-circles argument of `Fayed.PNTError.inv_zeta_bound_strip` (copied, then adapted)
for J = log L(c + ·, χ) - log L(c, χ), c = 3 + it. Inputs: ‖L(w, χ)‖ ≤ q‖w‖/Re w
(`norm_LFunction_le`, LOne.lean), so Re J ≤ log(3q(|t| + 8)) on |z| < 3 - 7/8 - η/4;
‖L'/L(w, χ)‖ ≤ Σ Λ(n)/n² on Re w ≥ 2 (uniform in χ); 2/3 ≤ ‖L(c, χ)‖ ≤ 2.
-/

open Complex Set Filter Topology
open scoped Real LSeries.notation ArithmeticFunction.vonMangoldt

namespace Fayed.QRHCor

open Fayed.PNTError

noncomputable section

variable {q : ℕ} [NeZero q]

/-- ‖L'/L(s, χ)‖ ≤ Σ Λ(n)/n² on Re s ≥ 2, uniformly in χ. -/
lemma norm_logDeriv_LFunction_le_K2 (χ : DirichletCharacter ℂ q) {s : ℂ} (hs : 2 ≤ s.re) :
    ‖deriv χ.LFunction s / χ.LFunction s‖ ≤ K2 := by
  have hs1 : 1 < s.re := by linarith
  rw [DirichletCharacter.deriv_LFunction_eq_deriv_LSeries χ hs1,
    DirichletCharacter.LFunction_eq_LSeries χ hs1, ← norm_neg, ← neg_div,
    ← DirichletCharacter.LSeries_twist_vonMangoldt_eq χ hs1]
  have hsum2 : Summable (fun n => ‖LSeries.term ↗Λ 2 n‖) :=
    summable_norm_iff.mpr (ArithmeticFunction.LSeriesSummable_vonMangoldt (by norm_num))
  have hle : ∀ n, ‖LSeries.term (↗χ * ↗Λ) s n‖ ≤ ‖LSeries.term ↗Λ 2 n‖ := by
    intro n
    refine (LSeries.norm_term_le s ?_).trans
      (LSeries.norm_term_le_of_re_le_re _ (by simpa using hs) n)
    rw [Pi.mul_apply, norm_mul]
    exact mul_le_of_le_one_left (norm_nonneg _) (DirichletCharacter.norm_le_one χ _)
  have hsums : Summable (fun n => ‖LSeries.term (↗χ * ↗Λ) s n‖) :=
    Summable.of_nonneg_of_le (fun _ => norm_nonneg _) hle hsum2
  calc ‖LSeries (↗χ * ↗Λ) s‖ = ‖∑' n, LSeries.term (↗χ * ↗Λ) s n‖ := rfl
    _ ≤ ∑' n, ‖LSeries.term (↗χ * ↗Λ) s n‖ := norm_tsum_le_tsum_norm hsums
    _ ≤ K2 := hsums.tsum_le_tsum hle hsum2

/-- ‖L(w, χ)‖ ≤ 2 on Re w ≥ 2. -/
lemma norm_LFunction_le_two (χ : DirichletCharacter ℂ q) {w : ℂ} (hw : 2 ≤ w.re) :
    ‖χ.LFunction w‖ ≤ 2 := by
  rw [DirichletCharacter.LFunction_eq_LSeries χ (by linarith)]
  have h := norm_LSeries_le_of_norm_le_one (f := fun n : ℕ => χ (n : ZMod q))
    (fun n => DirichletCharacter.norm_le_one χ _) hw
  have hπ : π ^ 2 / 6 < 2 := by nlinarith [Real.pi_lt_d2, Real.pi_pos]
  linarith

set_option maxHeartbeats 1000000 in
/-- Three circles for L(s, χ) on 7/8 + η ≤ Re s ≤ 2 (all t): both ‖L‖ and ‖1/L‖ are
≤ 2 exp(A (log Q)^θ₀), Q = q(|t| + 8), A = max K2 1 · (40/η + 1),
θ₀ = log(3 - 7/8 - η)/log(3 - 7/8 - η/2) < 1. -/
theorem LFunction_bound_strip {χ : DirichletCharacter ℂ q} (hχ : χ ≠ 1)
    {η : ℝ} (hη0 : 0 < η) (hη : η ≤ 1 / 2) {s : ℂ} (hs1 : 7 / 8 + η ≤ s.re) (hs2 : s.re ≤ 2) :
    ‖χ.LFunction s‖ ≤ 2 * Real.exp (max K2 1 * (40 / η + 1) *
        Real.log ((q : ℝ) * (|s.im| + 8)) ^
          (Real.log (3 - 7 / 8 - η) / Real.log (3 - 7 / 8 - η / 2))) ∧
    ‖(χ.LFunction s)⁻¹‖ ≤ 2 * Real.exp (max K2 1 * (40 / η + 1) *
        Real.log ((q : ℝ) * (|s.im| + 8)) ^
          (Real.log (3 - 7 / 8 - η) / Real.log (3 - 7 / 8 - η / 2))) := by
  obtain ⟨c, hc⟩ : ∃ c : ℂ, c = (3 : ℂ) + (s.im : ℂ) * I := ⟨_, rfl⟩
  have hcre : c.re = 3 := by rw [hc]; simp
  have hcim : c.im = s.im := by rw [hc]; simp
  obtain ⟨R, hR⟩ : ∃ R : ℝ, R = 3 - 7 / 8 - η / 4 := ⟨_, rfl⟩
  obtain ⟨r₃, hr₃⟩ : ∃ r : ℝ, r = 3 - 7 / 8 - η / 2 := ⟨_, rfl⟩
  obtain ⟨r₂, hr₂⟩ : ∃ r : ℝ, r = 3 - s.re := ⟨_, rfl⟩
  have hr₃1 : 1 < r₃ := by rw [hr₃]; linarith
  have hr₃R : r₃ < R := by rw [hr₃, hR]; linarith
  have hR52 : R ≤ 5 / 2 := by rw [hR]; linarith
  have hR0 : 0 < R := by linarith
  have hR1 : 1 < R := by linarith
  have hr₂1 : 1 ≤ r₂ := by rw [hr₂]; linarith
  have hr₂3 : r₂ ≤ r₃ - η / 2 := by rw [hr₂, hr₃]; linarith
  obtain ⟨Q, hQ⟩ : ∃ Q : ℝ, Q = (q : ℝ) * (|s.im| + 8) := ⟨_, rfl⟩
  have hq1 : (1 : ℝ) ≤ q := Nat.one_le_cast.mpr (NeZero.pos q)
  have hQ8 : 8 ≤ Q := by rw [hQ]; nlinarith [abs_nonneg s.im]
  have hcnorm : ‖c‖ ≤ 3 + |s.im| := by
    have h := Complex.norm_le_abs_re_add_abs_im c
    rw [hcre, hcim] at h
    norm_num at h
    linarith
  have hgeom : ∀ z ∈ Metric.closedBall (0 : ℂ) R,
      7 / 8 + η / 4 ≤ (c + z).re ∧ ‖c + z‖ ≤ |s.im| + 8 := by
    intro z hz
    rw [mem_closedBall_zero_iff] at hz
    have hre' := abs_le.mp ((Complex.abs_re_le_norm z).trans hz)
    refine ⟨?_, ?_⟩
    · rw [add_re, hcre]
      rw [hR] at hre'
      linarith [hre'.1]
    · calc ‖c + z‖ ≤ ‖c‖ + ‖z‖ := norm_add_le _ _
        _ ≤ (3 + |s.im|) + 5 / 2 := by linarith
        _ ≤ |s.im| + 8 := by linarith
  have hanal : AnalyticOnNhd ℂ (fun z => χ.LFunction (c + z)) (Metric.closedBall 0 R) := by
    intro z _
    exact ((DirichletCharacter.differentiable_LFunction hχ).analyticAt (c + z)).comp
      (analyticAt_const.add analyticAt_id)
  have hB0 : ∀ z ∈ Metric.closedBall (0 : ℂ) R, χ.LFunction (c + z) ≠ 0 := fun z hz =>
    OAI.DirichletCharacter.LFunction_ne_zero_of_seven_eighths_lt_re χ
      (by linarith [(hgeom z hz).1]) (fun h => hχ h.1)
  obtain ⟨J, hJan, hJ0, hJd, hJre⟩ :=
    LogOfAnalyticFunction (r := 1) one_pos hR1 hanal hB0
  have hlow : 2 / 3 ≤ ‖χ.LFunction c‖ := by
    have h := norm_LFunction_ge_of_one_lt_re χ (s := c) (by rw [hcre]; norm_num)
    rw [hcre] at h
    norm_num at h ⊢
    linarith
  have hcpos : 0 < ‖χ.LFunction c‖ := lt_of_lt_of_le (by norm_num) hlow
  have hc2 : ‖χ.LFunction c‖ ≤ 2 := norm_LFunction_le_two χ (by rw [hcre]; norm_num)
  obtain ⟨M, hM⟩ : ∃ M : ℝ, M = Real.log (3 * Q) := ⟨_, rfl⟩
  have hM0 : 0 < M := by rw [hM]; exact Real.log_pos (by linarith)
  have hJM : ∀ z ∈ Metric.ball (0 : ℂ) R, (J z).re ≤ M := by
    intro z hz
    have hz' := Metric.ball_subset_closedBall hz
    rw [← hJre z hz]
    obtain ⟨g1, g2⟩ := hgeom z hz'
    have hBz : 0 < ‖χ.LFunction (c + z)‖ := norm_pos_iff.mpr (hB0 z hz')
    have hup : ‖χ.LFunction (c + z)‖ ≤ 2 * Q := by
      have h := norm_LFunction_le hχ (s := c + z) (by linarith)
      refine h.trans ?_
      rw [div_le_iff₀ (by linarith), hQ]
      have hZ : 0 ≤ (q : ℝ) * (|s.im| + 8) := by positivity
      have e1 : (q : ℝ) * ‖c + z‖ ≤ q * (|s.im| + 8) :=
        mul_le_mul_of_nonneg_left g2 (by linarith)
      have e2 : (q : ℝ) * (|s.im| + 8) * (7 / 8) ≤ q * (|s.im| + 8) * (c + z).re :=
        mul_le_mul_of_nonneg_left (by linarith) hZ
      nlinarith
    show Real.log ‖χ.LFunction (c + z)‖ - Real.log ‖χ.LFunction (c + 0)‖ ≤ M
    rw [add_zero, ← Real.log_div hBz.ne' hcpos.ne', hM]
    apply Real.log_le_log (div_pos hBz hcpos)
    rw [div_le_iff₀ hcpos]
    have := mul_le_mul_of_nonneg_left hlow (by linarith : (0 : ℝ) ≤ 3 * Q)
    nlinarith
  -- Borel–Carathéodory on |z| ≤ r₃
  have hbig : ∀ z : ℂ, ‖z‖ ≤ r₃ → ‖J z‖ ≤ 2 * M * r₃ / (R - r₃) := fun z hz =>
    borelCaratheodory' hM0 hR0 hr₃R hJan.analyticOn hJ0 hJM (by rwa [mem_closedBall_zero_iff])
  -- the small disc |z| ≤ 1 lies in Re ≥ 2
  have hsmall : ∀ z : ℂ, ‖z‖ ≤ 1 → ‖J z‖ ≤ K2 := by
    intro z hz
    have hd : ∀ w ∈ Metric.closedBall (0 : ℂ) 1, DifferentiableAt ℂ J w := fun w hw =>
      (hJan w (Metric.closedBall_subset_ball (by linarith) hw)).differentiableAt
    have hbd : ∀ w ∈ Metric.closedBall (0 : ℂ) 1, ‖deriv J w‖ ≤ K2 := by
      intro w hw
      rw [hJd w hw, deriv_comp_const_add]
      apply norm_logDeriv_LFunction_le_K2
      rw [mem_closedBall_zero_iff] at hw
      rw [add_re, hcre]
      linarith [(abs_le.mp ((Complex.abs_re_le_norm w).trans hw)).1]
    have hmv := (convex_closedBall (0 : ℂ) 1).norm_image_sub_le_of_norm_deriv_le hd hbd
      (Metric.mem_closedBall_self zero_le_one) (by rwa [mem_closedBall_zero_iff])
    rw [hJ0, sub_zero, sub_zero] at hmv
    calc ‖J z‖ ≤ K2 * ‖z‖ := hmv
      _ ≤ K2 * 1 := mul_le_mul_of_nonneg_left hz K2_nonneg
      _ = K2 := mul_one _
  -- three circles at radius r₂
  have h3 := three_circles (J := J) (R := R) (r₂ := r₂) (r₃ := r₃) (a := max K2 1)
    (b := max (2 * M * r₃ / (R - r₃)) 1) hJan.differentiableOn hr₃1 hr₃R hr₂1 (by linarith)
    (fun z hz => (hsmall z hz.le).trans (le_max_left _ _))
    (fun z hz => (hbig z hz).trans (le_max_left _ _))
  have hsc : -(r₂ : ℂ) = s - c := by
    apply Complex.ext <;> simp [hr₂, hcre, hcim]
  rw [hsc] at h3
  -- exponents
  set θ := Real.log r₂ / Real.log r₃ with hθ
  set θ₀ := Real.log (3 - 7 / 8 - η) / Real.log (3 - 7 / 8 - η / 2) with hθ₀
  have hlr₃ : 0 < Real.log r₃ := Real.log_pos hr₃1
  have hθ0 : 0 ≤ θ := div_nonneg (Real.log_nonneg hr₂1) hlr₃.le
  have hθθ₀ : θ ≤ θ₀ := by
    rw [hθ, hθ₀, ← hr₃]
    apply div_le_div_of_nonneg_right _ hlr₃.le
    apply Real.log_le_log (by linarith)
    rw [hr₃] at hr₂3
    linarith
  have hθ₀0 : 0 ≤ θ₀ := le_trans hθ0 hθθ₀
  have hθ₀1 : θ₀ ≤ 1 := by
    rw [hθ₀, div_le_one (by rw [← hr₃]; exact hlr₃)]
    exact Real.log_le_log (by linarith) (by linarith)
  -- the log Q bound for b
  have hLQ : 1 ≤ Real.log Q := by
    rw [Real.le_log_iff_exp_le (by linarith)]
    have := Real.exp_one_lt_d9
    linarith
  have hM2 : M ≤ 2 * Real.log Q := by
    rw [hM, Real.log_mul (by norm_num) (by linarith)]
    have : Real.log 3 ≤ Real.log Q := Real.log_le_log (by norm_num) (by linarith)
    linarith
  have hRr₃ : R - r₃ = η / 4 := by rw [hR, hr₃]; ring
  have hbb : max (2 * M * r₃ / (R - r₃)) 1 ≤ (40 / η + 1) * Real.log Q := by
    have hr₃52 : r₃ ≤ 5 / 2 := by linarith
    have hq : 2 * M * r₃ / (R - r₃) ≤ 40 / η * Real.log Q := by
      rw [hRr₃, div_le_iff₀ (by positivity)]
      have : 2 * M * r₃ ≤ 2 * (2 * Real.log Q) * (5 / 2) := by
        apply mul_le_mul (by linarith) hr₃52 (by linarith) (by linarith)
      have e : 40 / η * Real.log Q * (η / 4) = 10 * Real.log Q := by
        field_simp
        ring
      rw [e]
      linarith
    have hpos : 0 ≤ 40 / η * Real.log Q := mul_nonneg (by positivity) (by linarith)
    apply max_le
    · nlinarith
    · nlinarith
  -- assemble the bound on ‖J (s - c)‖
  have hX : ‖J (s - c)‖ ≤ max K2 1 * (40 / η + 1) * Real.log Q ^ θ₀ := by
    have ha1 : 1 ≤ max K2 1 := le_max_right _ _
    have hb1 : 1 ≤ max (2 * M * r₃ / (R - r₃)) 1 := le_max_right _ _
    have hθ1 : θ ≤ 1 := le_trans hθθ₀ hθ₀1
    have e1 : max K2 1 ^ (1 - θ) ≤ max K2 1 := by
      calc max K2 1 ^ (1 - θ) ≤ max K2 1 ^ (1 : ℝ) :=
            Real.rpow_le_rpow_of_exponent_le ha1 (by linarith)
        _ = max K2 1 := Real.rpow_one _
    have e2 : max (2 * M * r₃ / (R - r₃)) 1 ^ θ ≤ max (2 * M * r₃ / (R - r₃)) 1 ^ θ₀ :=
      Real.rpow_le_rpow_of_exponent_le hb1 hθθ₀
    have e3 : max (2 * M * r₃ / (R - r₃)) 1 ^ θ₀ ≤ ((40 / η + 1) * Real.log Q) ^ θ₀ :=
      Real.rpow_le_rpow (by linarith) hbb hθ₀0
    have e4 : ((40 / η + 1) * Real.log Q) ^ θ₀ ≤ (40 / η + 1) * Real.log Q ^ θ₀ := by
      rw [Real.mul_rpow (by positivity) (by linarith)]
      apply mul_le_mul_of_nonneg_right _ (Real.rpow_nonneg (by linarith) _)
      calc (40 / η + 1) ^ θ₀ ≤ (40 / η + 1) ^ (1 : ℝ) :=
            Real.rpow_le_rpow_of_exponent_le
              (by have h40 : (0 : ℝ) ≤ 40 / η := (by positivity); linarith) hθ₀1
        _ = 40 / η + 1 := Real.rpow_one _
    have hbθ0 : 0 ≤ max (2 * M * r₃ / (R - r₃)) 1 ^ θ := Real.rpow_nonneg (by linarith) _
    calc ‖J (s - c)‖ ≤ max K2 1 ^ (1 - θ) * max (2 * M * r₃ / (R - r₃)) 1 ^ θ := h3
      _ ≤ max K2 1 * ((40 / η + 1) * Real.log Q ^ θ₀) :=
          mul_le_mul e1 (e2.trans (e3.trans e4)) hbθ0 (by linarith)
      _ = max K2 1 * (40 / η + 1) * Real.log Q ^ θ₀ := by ring
  -- from J to L and 1/L
  have hsball : s - c ∈ Metric.ball (0 : ℂ) R := by
    rw [mem_ball_zero_iff, ← hsc, norm_neg, Complex.norm_real, Real.norm_eq_abs,
      abs_of_pos (by linarith)]
    linarith
  have hre := hJre (s - c) hsball
  simp only [add_sub_cancel, add_zero] at hre
  have hzs : χ.LFunction s ≠ 0 :=
    OAI.DirichletCharacter.LFunction_ne_zero_of_seven_eighths_lt_re χ (by linarith)
      (fun h => hχ h.1)
  have hpos : 0 < ‖χ.LFunction s‖ := norm_pos_iff.mpr hzs
  set X := ‖J (s - c)‖ with hXdef
  have hreX : (J (s - c)).re ≤ X := Complex.re_le_norm _
  have hreX' : -X ≤ (J (s - c)).re := by
    have h1 := Complex.abs_re_le_norm (J (s - c))
    have h2 := neg_abs_le (J (s - c)).re
    linarith
  have hEX : Real.exp X ≤ Real.exp (max K2 1 * (40 / η + 1) * Real.log Q ^ θ₀) :=
    Real.exp_le_exp.mpr hX
  rw [← hQ]
  constructor
  · have hlogs : Real.log ‖χ.LFunction s‖ ≤ Real.log 2 + X := by
      have := Real.log_le_log hcpos hc2
      linarith
    calc ‖χ.LFunction s‖ = Real.exp (Real.log ‖χ.LFunction s‖) := (Real.exp_log hpos).symm
      _ ≤ Real.exp (Real.log 2 + X) := Real.exp_le_exp.mpr hlogs
      _ = 2 * Real.exp X := by rw [Real.exp_add, Real.exp_log (by norm_num)]
      _ ≤ 2 * Real.exp (max K2 1 * (40 / η + 1) * Real.log Q ^ θ₀) := by linarith
  · have hlogs : Real.log (2 / 3) - X ≤ Real.log ‖χ.LFunction s‖ := by
      have := Real.log_le_log (by norm_num) hlow
      linarith
    have hge : (2 / 3) * Real.exp (-X) ≤ ‖χ.LFunction s‖ := by
      calc (2 / 3) * Real.exp (-X) = Real.exp (Real.log (2 / 3) - X) := by
            rw [Real.exp_sub, Real.exp_log (by norm_num), Real.exp_neg]
            ring
        _ ≤ Real.exp (Real.log ‖χ.LFunction s‖) := Real.exp_le_exp.mpr hlogs
        _ = ‖χ.LFunction s‖ := Real.exp_log hpos
    rw [norm_inv]
    have h23 : 0 < (2 / 3) * Real.exp (-X) := by positivity
    have hEXp := Real.exp_pos X
    calc ‖χ.LFunction s‖⁻¹ ≤ ((2 / 3) * Real.exp (-X))⁻¹ := inv_anti₀ h23 hge
      _ = (3 / 2) * Real.exp X := by
          rw [mul_inv, Real.exp_neg, inv_inv]
          norm_num
      _ ≤ 2 * Real.exp (max K2 1 * (40 / η + 1) * Real.log Q ^ θ₀) := by linarith

/-- **Lindelöf for Dirichlet L-functions right of 7/8, uniformly in q and t:** for η, ε > 0
there is C with ‖L(s, χ)‖ ≤ C (q(|t| + 8))^ε and ‖L(s, χ)⁻¹‖ ≤ C (q(|t| + 8))^ε for all q,
all χ ≠ 1 mod q and all s with Re s ≥ 7/8 + η. -/
theorem LFunction_lindelof {η : ℝ} (hη0 : 0 < η) {ε : ℝ} (hε : 0 < ε) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ (q : ℕ) [NeZero q] (χ : DirichletCharacter ℂ q), χ ≠ 1 →
      ∀ s : ℂ, 7 / 8 + η ≤ s.re →
        ‖χ.LFunction s‖ ≤ C * ((q : ℝ) * (|s.im| + 8)) ^ ε ∧
          ‖(χ.LFunction s)⁻¹‖ ≤ C * ((q : ℝ) * (|s.im| + 8)) ^ ε := by
  obtain ⟨η', hη'⟩ : ∃ η' : ℝ, η' = min η (1 / 2) := ⟨_, rfl⟩
  have hη'0 : 0 < η' := by rw [hη']; exact lt_min hη0 (by norm_num)
  have hη'1 : η' ≤ 1 / 2 := by rw [hη']; exact min_le_right _ _
  have hη'η : η' ≤ η := by rw [hη']; exact min_le_left _ _
  obtain ⟨A, hA⟩ : ∃ A : ℝ, A = max K2 1 * (40 / η' + 1) := ⟨_, rfl⟩
  have hA0 : 0 ≤ A := by
    rw [hA]
    exact mul_nonneg (le_trans zero_le_one (le_max_right _ _)) (by positivity)
  obtain ⟨θ₀, hθ₀⟩ : ∃ θ : ℝ, θ = Real.log (3 - 7 / 8 - η') / Real.log (3 - 7 / 8 - η' / 2) :=
    ⟨_, rfl⟩
  have hl2 : 0 < Real.log (3 - 7 / 8 - η' / 2) := Real.log_pos (by linarith)
  have hθ₀0 : 0 ≤ θ₀ := by
    rw [hθ₀]
    exact div_nonneg (Real.log_nonneg (by linarith)) hl2.le
  have hθ₀1 : θ₀ < 1 := by
    rw [hθ₀, div_lt_one hl2]
    exact Real.log_lt_log (by linarith) (by linarith)
  obtain ⟨K, hK⟩ := rpow_le_mul_add hθ₀0 hθ₀1 hA0 hε
  have hK0 : 0 ≤ K := by
    have h := hK 0 le_rfl
    have : 0 ≤ A * (0 : ℝ) ^ θ₀ := mul_nonneg hA0 (Real.rpow_nonneg le_rfl _)
    linarith
  refine ⟨2 * Real.exp K, by positivity, ?_⟩
  intro q _ χ hχ s hs
  obtain ⟨Q, hQ⟩ : ∃ Q : ℝ, Q = (q : ℝ) * (|s.im| + 8) := ⟨_, rfl⟩
  rw [← hQ]
  have hq1 : (1 : ℝ) ≤ q := Nat.one_le_cast.mpr (NeZero.pos q)
  have hQ8 : 8 ≤ Q := by rw [hQ]; nlinarith [abs_nonneg s.im]
  have hQ0 : 0 < Q := by linarith
  have hlogQ : 0 ≤ Real.log Q := Real.log_nonneg (by linarith)
  have hQε : 1 ≤ Q ^ ε := Real.one_le_rpow (by linarith) hε.le
  have key : 2 * Real.exp (A * Real.log Q ^ θ₀) ≤ 2 * Real.exp K * Q ^ ε := by
    have h := hK (Real.log Q) hlogQ
    calc 2 * Real.exp (A * Real.log Q ^ θ₀) ≤ 2 * Real.exp (ε * Real.log Q + K) := by
          gcongr
      _ = 2 * Real.exp K * Q ^ ε := by
          rw [Real.exp_add, Real.rpow_def_of_pos hQ0, mul_comm ε]
          ring
  have hbase : 2 ≤ 2 * Real.exp K * Q ^ ε := by
    have h1 : 1 ≤ Real.exp K := Real.one_le_exp hK0
    nlinarith
  rcases le_or_gt s.re 2 with hs2 | hs2
  · have h := LFunction_bound_strip hχ hη'0 hη'1 (by linarith) hs2
    rw [← hA, ← hθ₀, ← hQ] at h
    exact ⟨h.1.trans key, h.2.trans key⟩
  · have h1 := norm_LFunction_le_two χ hs2.le
    have h2 := norm_LFunction_ge_of_one_lt_re χ (s := s) (by linarith)
    refine ⟨h1.trans hbase, ?_⟩
    have hhalf : 1 / 2 ≤ ‖χ.LFunction s‖ := by
      refine le_trans ?_ h2
      rw [le_div_iff₀ (by linarith)]
      linarith
    rw [norm_inv]
    have hinv : ‖χ.LFunction s‖⁻¹ ≤ 2 := by
      rw [inv_le_comm₀ (by linarith) (by norm_num)]
      linarith
    exact hinv.trans hbase

end

end Fayed.QRHCor
