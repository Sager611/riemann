import Riemann.Deep.PNTError.LogDeriv

/-!
# Power-saving PNT, part 4: 1/ζ(s) ≪ |t|^η to the right of a zero-free half-plane

If ζ has no zero in Re s > σ₀ (1/2 ≤ σ₀ < 1), then for every δ, η > 0 there is C with
‖1/ζ(s)‖ ≤ C |t|^η on Re s ≥ σ₀ + δ, |t| ≥ 4 (`inv_zeta_bound`).

Proof (Hadamard three circles, Titchmarsh §14.2). Around c = 3 + it let J = log ζ(c + ·) - log ζ(c)
(`LogOfAnalyticFunction`). On |z| ≤ 1 (inside Re ≥ 2) ‖J‖ ≤ Σ Λ(n)/n² (mean value theorem);
on |z| ≤ r₃ = 3 - σ₀ - δ/2, Borel–Carathéodory gives ‖J‖ ≪ log|t|. Three circles (`three_circles`,
from Mathlib's Hadamard three-lines theorem applied to J ∘ exp) at radius r₂ = 3 - Re s gives
‖J‖ ≪ (log|t|)^θ₀ with θ₀ = log(3 - σ₀ - δ)/log(3 - σ₀ - δ/2) < 1, hence
‖1/ζ(s)‖ ≤ 2 exp(A (log|t|)^θ₀) ≪_η |t|^η.
-/

open Complex Set Filter Topology
open scoped Real

namespace Fayed.PNTError

noncomputable section

/-- ‖ζ(w)‖ ≤ 15|Im w| + 2 for Re w ≥ 1/2, |Im w| ≥ 1. -/
lemma norm_zeta_le_crude' {w : ℂ} (h1 : 1 / 2 ≤ w.re) (him : 1 ≤ |w.im|) :
    ‖riemannZeta w‖ ≤ 15 * |w.im| + 2 := by
  rcases le_or_gt w.re 4 with h4 | h4
  · have := norm_zeta_le_crude h1 h4 him
    linarith
  · have hw : 2 ≤ w.re := by linarith
    have h1' : 1 < w.re := by linarith
    have h := norm_LSeries_le_of_norm_le_one (f := 1) (fun n => by simp) hw
    rw [LSeries_one_eq_riemannZeta h1'] at h
    have hπ : π ^ 2 / 6 < 2 := by nlinarith [Real.pi_lt_d2, Real.pi_pos]
    have : 0 ≤ 15 * |w.im| := by positivity
    linarith

/-- Hadamard's three circles theorem (radii 1 ≤ r₂ ≤ r₃), via Mathlib's three-lines theorem
for J ∘ exp on the strip 0 ≤ Re w ≤ log r₃. -/
lemma three_circles {J : ℂ → ℂ} {R r₂ r₃ a b : ℝ} (hJ : DifferentiableOn ℂ J (Metric.ball 0 R))
    (hr₃ : 1 < r₃) (hr₃R : r₃ < R) (hr₂1 : 1 ≤ r₂) (hr₂3 : r₂ ≤ r₃)
    (ha : ∀ z : ℂ, ‖z‖ = 1 → ‖J z‖ ≤ a) (hb : ∀ z : ℂ, ‖z‖ ≤ r₃ → ‖J z‖ ≤ b) :
    ‖J (-(r₂ : ℂ))‖ ≤ a ^ (1 - Real.log r₂ / Real.log r₃) * b ^ (Real.log r₂ / Real.log r₃) := by
  obtain ⟨u, hu⟩ : ∃ u : ℝ, u = Real.log r₃ := ⟨_, rfl⟩
  rw [← hu]
  have hu0 : 0 < u := by rw [hu]; exact Real.log_pos hr₃
  have hexpu : Real.exp u = r₃ := by rw [hu, Real.exp_log (by linarith)]
  have hdiff : ∀ w : ℂ, w.re ≤ u → DifferentiableAt ℂ (fun w => J (Complex.exp w)) w := by
    intro w hw
    have hnorm : ‖Complex.exp w‖ < R := by
      rw [Complex.norm_exp]
      calc Real.exp w.re ≤ Real.exp u := Real.exp_le_exp.mpr hw
        _ = r₃ := hexpu
        _ < R := hr₃R
    have hJd : DifferentiableAt ℂ J (Complex.exp w) :=
      hJ.differentiableAt (Metric.isOpen_ball.mem_nhds (by rwa [mem_ball_zero_iff]))
    exact hJd.comp w Complex.differentiableAt_exp
  have hcl : closure (Complex.HadamardThreeLines.verticalStrip 0 u) ⊆ {w : ℂ | w.re ≤ u} := by
    apply closure_minimal
    · intro w hw
      exact le_of_lt hw.2
    · exact isClosed_le Complex.continuous_re continuous_const
  have hd : DiffContOnCl ℂ (fun w => J (Complex.exp w))
      (Complex.HadamardThreeLines.verticalStrip 0 u) :=
    DifferentiableOn.diffContOnCl (fun w hw => (hdiff w (hcl hw)).differentiableWithinAt)
  have hbound : ∀ w : ℂ, w.re ≤ u → ‖J (Complex.exp w)‖ ≤ b := by
    intro w hw
    apply hb
    rw [Complex.norm_exp, ← hexpu]
    exact Real.exp_le_exp.mpr hw
  have hB : BddAbove ((norm ∘ (fun w => J (Complex.exp w))) ''
      Complex.HadamardThreeLines.verticalClosedStrip 0 u) := by
    refine ⟨b, ?_⟩
    rintro _ ⟨w, hw, rfl⟩
    exact hbound w hw.2
  have ha' : ∀ z ∈ re ⁻¹' {(0 : ℝ)}, ‖(fun w => J (Complex.exp w)) z‖ ≤ a := by
    intro z hz
    apply ha
    have hz' : z.re = 0 := hz
    rw [Complex.norm_exp, hz', Real.exp_zero]
  have hb' : ∀ z ∈ re ⁻¹' {u}, ‖(fun w => J (Complex.exp w)) z‖ ≤ b := fun z hz =>
    hbound z (le_of_eq hz)
  obtain ⟨w₀, hw₀⟩ : ∃ w : ℂ, w = (Real.log r₂ : ℂ) + π * I := ⟨_, rfl⟩
  have hw₀re : w₀.re = Real.log r₂ := by rw [hw₀]; simp
  have hmem : w₀ ∈ Complex.HadamardThreeLines.verticalClosedStrip 0 u := by
    show w₀.re ∈ Icc 0 u
    rw [hw₀re, hu]
    exact ⟨Real.log_nonneg hr₂1, Real.log_le_log (by linarith) hr₂3⟩
  have key := Complex.HadamardThreeLines.norm_le_interp_of_mem_verticalClosedStrip' hu0 hmem hd
    hB ha' hb'
  have hexp0 : Complex.exp w₀ = -(r₂ : ℂ) := by
    rw [hw₀, Complex.exp_add, ← Complex.ofReal_exp, Real.exp_log (by linarith),
      Complex.exp_pi_mul_I]
    ring
  simp only [hexp0, hw₀re, sub_zero] at key
  exact key

/-- If A ≥ 0, 0 ≤ θ < 1 and η > 0 then A u^θ ≤ η u + K for all u ≥ 0. -/
lemma rpow_le_mul_add {θ A η : ℝ} (hθ0 : 0 ≤ θ) (hθ1 : θ < 1) (hA : 0 ≤ A) (hη : 0 < η) :
    ∃ K : ℝ, ∀ u : ℝ, 0 ≤ u → A * u ^ θ ≤ η * u + K := by
  have h1θ : 0 < 1 - θ := by linarith
  obtain ⟨u₀, hu₀⟩ : ∃ u₀ : ℝ, u₀ = (A / η) ^ (1 / (1 - θ)) := ⟨_, rfl⟩
  have hAη : 0 ≤ A / η := div_nonneg hA hη.le
  have hu₀0 : 0 ≤ u₀ := by rw [hu₀]; exact Real.rpow_nonneg hAη _
  refine ⟨A * u₀ ^ θ, fun u hu => ?_⟩
  rcases le_or_gt u u₀ with h | h
  · have h1 : A * u ^ θ ≤ A * u₀ ^ θ :=
      mul_le_mul_of_nonneg_left (Real.rpow_le_rpow hu h hθ0) hA
    have h2 : 0 ≤ η * u := mul_nonneg hη.le hu
    linarith
  · have hu0 : 0 < u := lt_of_le_of_lt hu₀0 h
    have hsplit : u = u ^ θ * u ^ (1 - θ) := by
      rw [← Real.rpow_add hu0]
      simp
    have h3 : u₀ ^ (1 - θ) = A / η := by
      rw [hu₀, ← Real.rpow_mul hAη, show 1 / (1 - θ) * (1 - θ) = 1 by field_simp, Real.rpow_one]
    have h2 : u₀ ^ (1 - θ) ≤ u ^ (1 - θ) := Real.rpow_le_rpow hu₀0 h.le h1θ.le
    have hA' : A ≤ η * u ^ (1 - θ) := by
      rw [h3, div_le_iff₀ hη] at h2
      linarith
    have huθ : 0 ≤ u ^ θ := Real.rpow_nonneg hu _
    have h4 : A * u ^ θ ≤ η * u := by
      calc A * u ^ θ ≤ η * u ^ (1 - θ) * u ^ θ := mul_le_mul_of_nonneg_right hA' huθ
        _ = η * (u ^ θ * u ^ (1 - θ)) := by ring
        _ = η * u := by rw [← hsplit]
    have hK : 0 ≤ A * u₀ ^ θ := mul_nonneg hA (Real.rpow_nonneg hu₀0 _)
    linarith

lemma one_lt_log_four : 1 < Real.log 4 := by
  rw [Real.lt_log_iff_exp_lt (by norm_num)]
  have := Real.exp_one_lt_d9
  norm_num at this ⊢
  linarith

/-- The three-circles bound: ‖1/ζ(s)‖ ≤ 2 exp(A (log|t|)^θ₀) on σ₀ + δ ≤ Re s ≤ 2, |t| ≥ 4. -/
theorem inv_zeta_bound_strip {σ₀ : ℝ} (hZ : ZeroFree σ₀) (h12 : 1 / 2 ≤ σ₀) (h1 : σ₀ < 1)
    {δ : ℝ} (hδ0 : 0 < δ) (hδ : δ ≤ 1 / 2) {s : ℂ} (hs1 : σ₀ + δ ≤ s.re) (hs2 : s.re ≤ 2)
    (ht : 4 ≤ |s.im|) :
    ‖(riemannZeta s)⁻¹‖ ≤ 2 * Real.exp (max K2 1 * (100 / δ + 1) *
      Real.log |s.im| ^ (Real.log (3 - σ₀ - δ) / Real.log (3 - σ₀ - δ / 2))) := by
  obtain ⟨c, hc⟩ : ∃ c : ℂ, c = (3 : ℂ) + (s.im : ℂ) * I := ⟨_, rfl⟩
  have hcre : c.re = 3 := by rw [hc]; simp
  have hcim : c.im = s.im := by rw [hc]; simp
  obtain ⟨R, hR⟩ : ∃ R : ℝ, R = 3 - σ₀ - δ / 4 := ⟨_, rfl⟩
  obtain ⟨r₃, hr₃⟩ : ∃ r : ℝ, r = 3 - σ₀ - δ / 2 := ⟨_, rfl⟩
  obtain ⟨r₂, hr₂⟩ : ∃ r : ℝ, r = 3 - s.re := ⟨_, rfl⟩
  have hr₃1 : 1 < r₃ := by rw [hr₃]; linarith
  have hr₃R : r₃ < R := by rw [hr₃, hR]; linarith
  have hR52 : R ≤ 5 / 2 := by rw [hR]; linarith
  have hR0 : 0 < R := by linarith
  have hr₂1 : 1 ≤ r₂ := by rw [hr₂]; linarith
  have hr₂3 : r₂ ≤ r₃ - δ / 2 := by rw [hr₂, hr₃]; linarith
  have hgeom : ∀ z ∈ Metric.closedBall (0 : ℂ) R,
      σ₀ + δ / 4 ≤ (c + z).re ∧ |s.im| - 5 / 2 ≤ |(c + z).im| ∧
        |(c + z).im| ≤ |s.im| + 5 / 2 := by
    intro z hz
    rw [mem_closedBall_zero_iff] at hz
    have hre := Complex.abs_re_le_norm z
    have him := Complex.abs_im_le_norm z
    rw [add_re, add_im, hcre, hcim]
    have hre' := abs_le.mp (hre.trans hz)
    have him' : |z.im| ≤ 5 / 2 := (him.trans hz).trans hR52
    refine ⟨by rw [hR] at hre'; linarith [hre'.1], ?_, ?_⟩
    · have h := abs_add_le (s.im + z.im) (-z.im)
      rw [add_neg_cancel_right, abs_neg] at h
      linarith
    · have h := abs_add_le s.im z.im
      linarith
  have hne1 : ∀ z ∈ Metric.closedBall (0 : ℂ) R, c + z ≠ 1 := by
    intro z hz h
    have h' := (hgeom z hz).2.1
    rw [h, one_im, abs_zero] at h'
    linarith
  have hanal : AnalyticOnNhd ℂ (fun z => riemannZeta (c + z)) (Metric.closedBall 0 R) := by
    intro z hz
    exact (analyticAt_zeta (hne1 z hz)).comp (analyticAt_const.add analyticAt_id)
  have hB0 : ∀ z ∈ Metric.closedBall (0 : ℂ) R, riemannZeta (c + z) ≠ 0 := fun z hz =>
    hZ _ (by linarith [(hgeom z hz).1])
  obtain ⟨J, hJan, hJ0, hJd, hJre⟩ :=
    LogOfAnalyticFunction (r := 1) one_pos (by linarith) hanal hB0
  obtain ⟨M, hM⟩ : ∃ M : ℝ, M = Real.log (50 * |s.im|) := ⟨_, rfl⟩
  have hlow : 1 / 2 ≤ ‖riemannZeta c‖ := norm_zeta_ge_half (by rw [hcre]; norm_num)
  have hM0 : 0 < M := by rw [hM]; exact Real.log_pos (by linarith)
  have hJM : ∀ z ∈ Metric.ball (0 : ℂ) R, (J z).re ≤ M := by
    intro z hz
    have hz' := Metric.ball_subset_closedBall hz
    rw [← hJre z hz]
    obtain ⟨g1, g2, g3⟩ := hgeom z hz'
    have hBz : 0 < ‖riemannZeta (c + z)‖ := norm_pos_iff.mpr (hB0 z hz')
    have hup : ‖riemannZeta (c + z)‖ ≤ 25 * |s.im| := by
      have := norm_zeta_le_crude' (w := c + z) (by linarith) (by linarith)
      linarith
    have l1 := Real.log_le_log hBz hup
    have l2 : Real.log (1 / 2) ≤ Real.log ‖riemannZeta c‖ := Real.log_le_log (by norm_num) hlow
    rw [one_div, Real.log_inv] at l2
    have l3 : Real.log (50 * |s.im|) = Real.log 2 + Real.log (25 * |s.im|) := by
      rw [← Real.log_mul (by norm_num) (by positivity)]
      ring_nf
    show Real.log ‖riemannZeta (c + z)‖ - Real.log ‖riemannZeta (c + 0)‖ ≤ M
    rw [add_zero, hM, l3]
    linarith
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
      apply norm_logDeriv_le_K2
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
  set θ₀ := Real.log (3 - σ₀ - δ) / Real.log (3 - σ₀ - δ / 2) with hθ₀
  have hlr₃ : 0 < Real.log r₃ := Real.log_pos hr₃1
  have hθ0 : 0 ≤ θ := div_nonneg (Real.log_nonneg hr₂1) hlr₃.le
  have hθ1 : θ ≤ 1 := by
    rw [hθ, div_le_one hlr₃]
    exact Real.log_le_log (by linarith) (by linarith)
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
  -- the log|t| bound for b
  have hL4 : Real.log 4 ≤ Real.log |s.im| := Real.log_le_log (by norm_num) ht
  have hL1 : 1 ≤ Real.log |s.im| := le_trans one_lt_log_four.le hL4
  have hM5 : M ≤ 5 * Real.log |s.im| := by
    have ht0 : 0 < |s.im| := by linarith
    have h256 : (64 : ℝ) ≤ |s.im| ^ 3 := by
      have := pow_le_pow_left₀ (by norm_num : (0 : ℝ) ≤ 4) ht 3
      norm_num at this
      linarith
    have hA : 50 * |s.im| ≤ |s.im| ^ 5 := by
      have e5 : |s.im| ^ 5 = |s.im| * |s.im| * |s.im| ^ 3 := by ring
      rw [e5]
      nlinarith
    rw [hM]
    calc Real.log (50 * |s.im|) ≤ Real.log (|s.im| ^ 5) := Real.log_le_log (by positivity) hA
      _ = 5 * Real.log |s.im| := by rw [Real.log_pow]; norm_num
  have hRr₃ : R - r₃ = δ / 4 := by rw [hR, hr₃]; ring
  have hbb : max (2 * M * r₃ / (R - r₃)) 1 ≤ (100 / δ + 1) * Real.log |s.im| := by
    have hr₃52 : r₃ ≤ 5 / 2 := by linarith
    have hq : 2 * M * r₃ / (R - r₃) ≤ 100 / δ * Real.log |s.im| := by
      rw [hRr₃, div_le_iff₀ (by positivity)]
      have : 2 * M * r₃ ≤ 2 * (5 * Real.log |s.im|) * (5 / 2) := by
        apply mul_le_mul (by linarith) hr₃52 (by linarith) (by linarith)
      have e : 100 / δ * Real.log |s.im| * (δ / 4) = 25 * Real.log |s.im| := by
        field_simp
        ring
      rw [e]
      linarith
    have hpos : 0 ≤ 100 / δ * Real.log |s.im| := by positivity
    apply max_le
    · nlinarith
    · nlinarith
  -- assemble the bound on ‖J (s - c)‖
  have hX : ‖J (s - c)‖ ≤ max K2 1 * (100 / δ + 1) * Real.log |s.im| ^ θ₀ := by
    have ha1 : 1 ≤ max K2 1 := le_max_right _ _
    have hb1 : 1 ≤ max (2 * M * r₃ / (R - r₃)) 1 := le_max_right _ _
    have e1 : max K2 1 ^ (1 - θ) ≤ max K2 1 := by
      calc max K2 1 ^ (1 - θ) ≤ max K2 1 ^ (1 : ℝ) :=
            Real.rpow_le_rpow_of_exponent_le ha1 (by linarith)
        _ = max K2 1 := Real.rpow_one _
    have e2 : max (2 * M * r₃ / (R - r₃)) 1 ^ θ ≤ max (2 * M * r₃ / (R - r₃)) 1 ^ θ₀ :=
      Real.rpow_le_rpow_of_exponent_le hb1 hθθ₀
    have e3 : max (2 * M * r₃ / (R - r₃)) 1 ^ θ₀ ≤ ((100 / δ + 1) * Real.log |s.im|) ^ θ₀ :=
      Real.rpow_le_rpow (by linarith) hbb hθ₀0
    have e4 : ((100 / δ + 1) * Real.log |s.im|) ^ θ₀ ≤
        (100 / δ + 1) * Real.log |s.im| ^ θ₀ := by
      rw [Real.mul_rpow (by positivity) (by linarith)]
      apply mul_le_mul_of_nonneg_right _ (Real.rpow_nonneg (by linarith) _)
      calc (100 / δ + 1) ^ θ₀ ≤ (100 / δ + 1) ^ (1 : ℝ) :=
            Real.rpow_le_rpow_of_exponent_le (by have h100 : (0 : ℝ) ≤ 100 / δ := (by positivity); linarith)
              hθ₀1
        _ = 100 / δ + 1 := Real.rpow_one _
    have hbθ0 : 0 ≤ max (2 * M * r₃ / (R - r₃)) 1 ^ θ := Real.rpow_nonneg (by linarith) _
    calc ‖J (s - c)‖ ≤ max K2 1 ^ (1 - θ) * max (2 * M * r₃ / (R - r₃)) 1 ^ θ := h3
      _ ≤ max K2 1 * ((100 / δ + 1) * Real.log |s.im| ^ θ₀) :=
          mul_le_mul e1 (e2.trans (e3.trans e4)) hbθ0 (by linarith)
      _ = max K2 1 * (100 / δ + 1) * Real.log |s.im| ^ θ₀ := by ring
  -- from J to 1/ζ
  have hsball : s - c ∈ Metric.ball (0 : ℂ) R := by
    rw [mem_ball_zero_iff, ← hsc, norm_neg, Complex.norm_real, Real.norm_eq_abs,
      abs_of_pos (by linarith)]
    linarith
  have hre := hJre (s - c) hsball
  simp only [add_sub_cancel, add_zero] at hre
  have hreJ : -‖J (s - c)‖ ≤ (J (s - c)).re := by
    have := Complex.abs_re_le_norm (J (s - c))
    linarith [neg_abs_le (J (s - c)).re]
  have l2 : Real.log (1 / 2) ≤ Real.log ‖riemannZeta c‖ := Real.log_le_log (by norm_num) hlow
  rw [one_div, Real.log_inv] at l2
  have hzs : riemannZeta s ≠ 0 := hZ s (by linarith)
  have hpos : 0 < ‖riemannZeta s‖ := norm_pos_iff.mpr hzs
  set X := ‖J (s - c)‖ with hXdef
  have hlogs : -Real.log 2 - X ≤ Real.log ‖riemannZeta s‖ := by linarith
  have hexp : Real.exp (-Real.log 2 - X) ≤ ‖riemannZeta s‖ := by
    calc Real.exp (-Real.log 2 - X) ≤ Real.exp (Real.log ‖riemannZeta s‖) :=
          Real.exp_le_exp.mpr hlogs
      _ = ‖riemannZeta s‖ := Real.exp_log hpos
  rw [norm_inv]
  calc ‖riemannZeta s‖⁻¹ ≤ (Real.exp (-Real.log 2 - X))⁻¹ := inv_anti₀ (Real.exp_pos _) hexp
    _ = 2 * Real.exp X := by
        rw [← Real.exp_neg, neg_sub, sub_neg_eq_add, Real.exp_add, Real.exp_log (by norm_num)]
        ring
    _ ≤ 2 * Real.exp (max K2 1 * (100 / δ + 1) * Real.log |s.im| ^ θ₀) := by
        gcongr

/-- **1/ζ(s) ≪_η |t|^η** on Re s ≥ σ₀ + δ, |t| ≥ 4. -/
theorem inv_zeta_bound {σ₀ : ℝ} (hZ : ZeroFree σ₀) (h12 : 1 / 2 ≤ σ₀) (h1 : σ₀ < 1)
    {δ : ℝ} (hδ0 : 0 < δ) {η : ℝ} (hη : 0 < η) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ s : ℂ, σ₀ + δ ≤ s.re → 4 ≤ |s.im| →
      ‖(riemannZeta s)⁻¹‖ ≤ C * |s.im| ^ η := by
  obtain ⟨δ', hδ'⟩ : ∃ δ' : ℝ, δ' = min δ (1 / 2) := ⟨_, rfl⟩
  have hδ'0 : 0 < δ' := by rw [hδ']; exact lt_min hδ0 (by norm_num)
  have hδ'1 : δ' ≤ 1 / 2 := by rw [hδ']; exact min_le_right _ _
  have hδ'δ : δ' ≤ δ := by rw [hδ']; exact min_le_left _ _
  obtain ⟨A, hA⟩ : ∃ A : ℝ, A = max K2 1 * (100 / δ' + 1) := ⟨_, rfl⟩
  have hA0 : 0 ≤ A := by
    rw [hA]
    exact mul_nonneg (le_trans zero_le_one (le_max_right _ _)) (by positivity)
  obtain ⟨θ₀, hθ₀⟩ : ∃ θ : ℝ, θ = Real.log (3 - σ₀ - δ') / Real.log (3 - σ₀ - δ' / 2) :=
    ⟨_, rfl⟩
  have hl2 : 0 < Real.log (3 - σ₀ - δ' / 2) := Real.log_pos (by linarith)
  have hθ₀0 : 0 ≤ θ₀ := by
    rw [hθ₀]
    exact div_nonneg (Real.log_nonneg (by linarith)) hl2.le
  have hθ₀1 : θ₀ < 1 := by
    rw [hθ₀, div_lt_one hl2]
    exact Real.log_lt_log (by linarith) (by linarith)
  obtain ⟨K, hK⟩ := rpow_le_mul_add hθ₀0 hθ₀1 hA0 hη
  refine ⟨2 * Real.exp K, by positivity, ?_⟩
  intro s hs ht
  have ht0 : 0 < |s.im| := by linarith
  have hlogt : 0 ≤ Real.log |s.im| := Real.log_nonneg (by linarith)
  have key : 2 * Real.exp (A * Real.log |s.im| ^ θ₀) ≤ 2 * Real.exp K * |s.im| ^ η := by
    have h := hK (Real.log |s.im|) hlogt
    calc 2 * Real.exp (A * Real.log |s.im| ^ θ₀) ≤ 2 * Real.exp (η * Real.log |s.im| + K) := by
          gcongr
      _ = 2 * Real.exp K * |s.im| ^ η := by
          rw [Real.exp_add, Real.rpow_def_of_pos ht0, mul_comm η]
          ring
  rcases le_or_gt s.re 2 with hs2 | hs2
  · have h := inv_zeta_bound_strip hZ h12 h1 hδ'0 hδ'1 (by linarith) hs2 ht
    rw [← hA, ← hθ₀] at h
    exact h.trans key
  · have h := norm_zeta_ge_half hs2.le
    have hpos : 0 < ‖riemannZeta s‖ := by linarith
    rw [norm_inv]
    have h2 : ‖riemannZeta s‖⁻¹ ≤ 2 := by
      rw [inv_le_comm₀ hpos (by norm_num)]
      linarith
    have h3 : 1 ≤ Real.exp (A * Real.log |s.im| ^ θ₀) :=
      Real.one_le_exp (mul_nonneg hA0 (Real.rpow_nonneg hlogt _))
    linarith

end

end Fayed.PNTError
