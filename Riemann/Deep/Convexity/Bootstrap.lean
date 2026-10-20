import Riemann.Deep.Convexity.FE

/-!
# Bootstrap maximum principle on the half-strip 0 ≤ Re s ≤ 1+δ, Im s ≥ 1

g(s) = ζ(s)·(-is)^{(s-1-δ)/2}·e^{εs²}. On Re s = 0, |g(it)| ≤ (4/5)|g(1+it)| (left_edge);
on Re s = 1+δ, g is bounded; the Gaussian kills the top edge (a-priori bound apriori).
Maximum modulus on rectangles + iteration ⇒ |g| ≤ 5K on the whole half-strip, K independent of ε.
-/

open Complex Filter Topology Set MeasureTheory Asymptotics
open scoped Real ComplexConjugate

namespace Fayed.Convexity

noncomputable section

/-- Weight (-is)^{(s-1-δ)/2}. -/
def wt (δ : ℝ) (s : ℂ) : ℂ := (-I * s) ^ ((s - (1 + δ)) / 2)

/-- Gaussian e^{εs²}. -/
def gs (ε : ℝ) (s : ℂ) : ℂ := Complex.exp (ε * s ^ 2)

/-- The auxiliary function. -/
def gf (δ ε : ℝ) (s : ℂ) : ℂ := riemannZeta s * wt δ s * gs ε s

/-- Half-strip. -/
def H (δ : ℝ) : Set ℂ := {s | 0 ≤ s.re ∧ s.re ≤ 1 + δ ∧ 1 ≤ s.im}

lemma norm_gs (ε : ℝ) (s : ℂ) : ‖gs ε s‖ = Real.exp (ε * (s.re ^ 2 - s.im ^ 2)) := by
  rw [gs, Complex.norm_exp]; congr 1; simp [sq]

lemma negI_mul_re (s : ℂ) : (-I * s).re = s.im := by simp
lemma negI_mul_im (s : ℂ) : (-I * s).im = -s.re := by simp

lemma arg_negI_le {s : ℂ} (h0 : 0 ≤ s.re) (h1 : 0 < s.im) : arg (-I * s) ≤ 0 := by
  rcases h0.lt_or_eq with h | h
  · exact (Complex.arg_neg_iff.mpr (by rw [negI_mul_im]; linarith)).le
  · exact (Complex.arg_eq_zero_iff.mpr ⟨by rw [negI_mul_re]; linarith,
      by rw [negI_mul_im]; linarith⟩).le

lemma arg_negI_ge {s : ℂ} (h0 : 0 ≤ s.re) (h1 : 0 < s.im) :
    -(s.re / s.im) ≤ arg (-I * s) := by
  have hle := arg_negI_le h0 h1
  have habs : |arg (-I * s)| < π / 2 :=
    Complex.abs_arg_lt_pi_div_two_iff.mpr (Or.inl (by rw [negI_mul_re]; exact h1))
  have htan : Real.tan (arg (-I * s)) = -(s.re / s.im) := by
    rw [Complex.tan_arg, negI_mul_re, negI_mul_im, neg_div]
  rcases hle.lt_or_eq with h | h
  · have hlt : -arg (-I * s) < Real.tan (-arg (-I * s)) :=
      Real.lt_tan (by linarith) (by linarith [(abs_lt.mp habs).1])
    rw [Real.tan_neg, htan] at hlt; linarith
  · rw [h]; have : 0 ≤ s.re / s.im := div_nonneg h0 h1.le; linarith

lemma norm_wt (δ : ℝ) {s : ℂ} (h1 : 0 < s.im) :
    ‖wt δ s‖ = ‖s‖ ^ ((s.re - (1 + δ)) / 2) / Real.exp (arg (-I * s) * (s.im / 2)) := by
  have hB : -I * s ≠ 0 := by
    intro h; have := congrArg Complex.re h; rw [negI_mul_re] at this; simp at this; linarith
  rw [wt, Complex.norm_cpow_of_ne_zero hB]
  have hn : ‖-I * s‖ = ‖s‖ := by simp
  have hre : ((s - (1 + δ)) / 2).re = (s.re - (1 + δ)) / 2 := by simp
  have him : ((s - (1 + δ)) / 2).im = s.im / 2 := by simp
  rw [hn, hre, him]

lemma wt_ge (δ : ℝ) {s : ℂ} (h0 : 0 ≤ s.re) (h1 : 0 < s.im) :
    ‖s‖ ^ ((s.re - (1 + δ)) / 2) ≤ ‖wt δ s‖ := by
  rw [norm_wt δ h1]
  refine le_div_self (Real.rpow_nonneg (norm_nonneg _) _) (Real.exp_pos _) ?_
  rw [Real.exp_le_one_iff]
  exact mul_nonpos_of_nonpos_of_nonneg (arg_negI_le h0 h1) (by linarith)

lemma wt_le (δ : ℝ) {s : ℂ} (h0 : 0 ≤ s.re) (h1 : 0 < s.im) :
    ‖wt δ s‖ ≤ ‖s‖ ^ ((s.re - (1 + δ)) / 2) * Real.exp (s.re / 2) := by
  rw [norm_wt δ h1]
  have hge : Real.exp (-(s.re / 2)) ≤ Real.exp (arg (-I * s) * (s.im / 2)) := by
    apply Real.exp_le_exp.mpr
    have := arg_negI_ge h0 h1
    have e : -(s.re / s.im) * (s.im / 2) = -(s.re / 2) := by field_simp
    nlinarith
  calc ‖s‖ ^ ((s.re - (1 + δ)) / 2) / Real.exp (arg (-I * s) * (s.im / 2))
      ≤ ‖s‖ ^ ((s.re - (1 + δ)) / 2) / Real.exp (-(s.re / 2)) :=
        div_le_div_of_nonneg_left (Real.rpow_nonneg (norm_nonneg _) _) (Real.exp_pos _) hge
    _ = _ := by rw [Real.exp_neg, div_inv_eq_mul]

lemma one_le_norm_of_im {s : ℂ} (h : 1 ≤ s.im) : 1 ≤ ‖s‖ :=
  le_trans (by rw [abs_of_nonneg (by linarith)]; exact h) (abs_im_le_norm s)

lemma wt_le_H {δ : ℝ} (hδ : δ ≤ 1) {s : ℂ} (hs : s ∈ H δ) : ‖wt δ s‖ ≤ Real.exp 1 := by
  obtain ⟨h0, h1, h2⟩ := hs
  refine (wt_le δ h0 (by linarith)).trans ?_
  have hr : ‖s‖ ^ ((s.re - (1 + δ)) / 2) ≤ 1 :=
    Real.rpow_le_one_of_one_le_of_nonpos (one_le_norm_of_im h2) (by linarith)
  have he : Real.exp (s.re / 2) ≤ Real.exp 1 := Real.exp_le_exp.mpr (by linarith)
  calc ‖s‖ ^ ((s.re - (1 + δ)) / 2) * Real.exp (s.re / 2) ≤ 1 * Real.exp 1 := by
        gcongr
    _ = Real.exp 1 := one_mul _

lemma gs_le_H {δ ε : ℝ} (hδ : δ ≤ 1) (_hε0 : 0 < ε) (hε1 : ε ≤ 1) {s : ℂ} (hs : s ∈ H δ) :
    ‖gs ε s‖ ≤ Real.exp (4 - ε * s.im ^ 2) := by
  obtain ⟨h0, h1, h2⟩ := hs
  rw [norm_gs]; apply Real.exp_le_exp.mpr
  have : ε * s.re ^ 2 ≤ 4 := by nlinarith
  nlinarith

lemma gf_le_H {δ ε : ℝ} (hδ : δ ≤ 1) (hε0 : 0 < ε) (hε1 : ε ≤ 1) {s : ℂ} (hs : s ∈ H δ) :
    ‖gf δ ε s‖ ≤ ‖riemannZeta s‖ * Real.exp 1 * Real.exp (4 - ε * s.im ^ 2) := by
  rw [gf, norm_mul, norm_mul]
  gcongr
  · exact wt_le_H hδ hs
  · exact gs_le_H hδ hε0 hε1 hs

lemma gf_gen {δ ε : ℝ} (hδ : δ ≤ 1) (hε0 : 0 < ε) (hε1 : ε ≤ 1) {s : ℂ} (hs : s ∈ H δ) :
    ‖gf δ ε s‖ ≤ 24 * Real.exp (5 + 3 * s.im - ε * s.im ^ 2) := by
  refine (gf_le_H hδ hε0 hε1 hs).trans ?_
  have hz := apriori hs.1 (by linarith [hs.2.1]) hs.2.2
  calc ‖riemannZeta s‖ * Real.exp 1 * Real.exp (4 - ε * s.im ^ 2)
      ≤ 24 * Real.exp (3 * s.im) * Real.exp 1 * Real.exp (4 - ε * s.im ^ 2) := by gcongr
    _ = 24 * Real.exp (5 + 3 * s.im - ε * s.im ^ 2) := by
        rw [mul_assoc, mul_assoc, ← Real.exp_add, ← Real.exp_add]; ring_nf

lemma differentiableAt_gf (δ ε : ℝ) {s : ℂ} (h : 0 < s.im) : DifferentiableAt ℂ (gf δ ε) s := by
  have h1 : s ≠ 1 := by intro e; rw [e] at h; simp at h
  have hz := differentiableAt_riemannZeta h1
  have hw : DifferentiableAt ℂ (wt δ) s := by
    unfold wt
    refine DifferentiableAt.cpow (by fun_prop) (by fun_prop) ?_
    exact Complex.mem_slitPlane_iff.mpr (Or.inl (by rw [negI_mul_re]; exact h))
  have hg : DifferentiableAt ℂ (gs ε) s := by unfold gs; fun_prop
  show DifferentiableAt ℂ (fun s => riemannZeta s * wt δ s * gs ε s) s
  exact (hz.mul hw).mul hg

lemma zeta_right_bound {δ : ℝ} (hδ : 0 < δ) {s : ℂ} (hs : s.re = 1 + δ) :
    ‖riemannZeta s‖ ≤ ∑' n : ℕ, 1 / ((n : ℝ) + 1) ^ (1 + δ) := by
  have hre : 1 < s.re := by linarith
  rw [zeta_eq_tsum_one_div_nat_add_one_cpow hre]
  have hterm : ∀ n : ℕ, ‖1 / ((n : ℂ) + 1) ^ s‖ = 1 / ((n : ℝ) + 1) ^ (1 + δ) := by
    intro n
    rw [norm_div, norm_one, show ((n : ℂ) + 1) = (((n : ℝ) + 1 : ℝ) : ℂ) by push_cast; ring,
      Complex.norm_cpow_eq_rpow_re_of_pos (by positivity), hs]
  have hsum : Summable (fun n : ℕ => 1 / ((n : ℝ) + 1) ^ (1 + δ)) := by
    have := (Real.summable_one_div_nat_add_rpow 1 (1 + δ)).mpr (by linarith)
    refine this.congr (fun n => ?_)
    rw [abs_of_nonneg (by positivity)]
  refine (norm_tsum_le_tsum_norm ?_).trans (le_of_eq ?_)
  · simpa only [hterm] using hsum
  · exact tsum_congr hterm

/-- Left edge transfer: |g(it)| ≤ (4/5)|g(1+it)|. -/
lemma gf_left {δ ε : ℝ} (hδ0 : 0 < δ) (hδ1 : δ ≤ 1) (hε : 0 ≤ ε) {t : ℝ} (ht : 1 ≤ t) :
    ‖gf δ ε (t * I)‖ ≤ (4 / 5) * ‖gf δ ε (1 + t * I)‖ := by
  have ht0 : 0 < t := by linarith
  have hre0 : (↑t * I : ℂ).re = 0 := by simp
  have him0 : (↑t * I : ℂ).im = t := by simp
  have hre1 : (1 + ↑t * I : ℂ).re = 1 := by simp
  have him1 : (1 + ↑t * I : ℂ).im = t := by simp
  have hw0 : ‖wt δ (t * I)‖ = t ^ (-(1 + δ) / 2) := by
    rw [norm_wt δ (by rw [him0]; exact ht0)]
    have harg : arg (-I * (t * I)) = 0 := Complex.arg_eq_zero_iff.mpr
      ⟨by rw [negI_mul_re, him0]; exact ht0.le, by rw [negI_mul_im, hre0]; simp⟩
    rw [harg, zero_mul, Real.exp_zero, div_one, hre0]
    have hn : ‖(↑t * I : ℂ)‖ = t := by simp [abs_of_pos ht0]
    rw [hn]; congr 1; ring
  have hw1 := wt_ge δ (s := 1 + t * I) (by rw [hre1]; norm_num) (by rw [him1]; exact ht0)
  rw [hre1] at hw1
  have hg0 : ‖gs ε (t * I)‖ = Real.exp (ε * (0 ^ 2 - t ^ 2)) := by rw [norm_gs, hre0, him0]
  have hg1 : Real.exp (ε * (0 ^ 2 - t ^ 2)) ≤ ‖gs ε (1 + t * I)‖ := by
    rw [norm_gs, hre1, him1]
    apply Real.exp_le_exp.mpr
    nlinarith
  have hz := left_edge ht0
  have hsq : √((4 / 25) * t) = (2 / 5) * t ^ (1 / 2 : ℝ) := by
    rw [Real.sqrt_mul (by norm_num), Real.sqrt_eq_rpow t]
    congr 1
    rw [show (4 / 25 : ℝ) = (2 / 5) ^ 2 by norm_num]; exact Real.sqrt_sq (by norm_num)
  rw [hsq] at hz
  have hn1 : ‖(1 + ↑t * I : ℂ)‖ ≤ 2 * t := by
    refine (norm_le_abs_re_add_abs_im _).trans ?_
    rw [hre1, him1, abs_one, abs_of_pos ht0]; linarith
  have hn1pos : 0 < ‖(1 + ↑t * I : ℂ)‖ := by
    have := one_le_norm_of_im (s := 1 + t * I) (by rw [him1]; exact ht); linarith
  have hA : (2 * t) ^ (-δ / 2) ≤ ‖(1 + ↑t * I : ℂ)‖ ^ ((1 - (1 + δ)) / 2) := by
    rw [show (1 - (1 + δ)) / 2 = -δ / 2 by ring]
    exact Real.rpow_le_rpow_of_nonpos hn1pos hn1 (by linarith)
  have hB : (2 * t) ^ (-δ / 2) = (2 : ℝ) ^ (-δ / 2) * t ^ (-δ / 2) :=
    Real.mul_rpow (by norm_num) ht0.le
  have hC : (1 / 2 : ℝ) ≤ (2 : ℝ) ^ (-δ / 2) := by
    have := Real.rpow_le_rpow_of_exponent_le (by norm_num : (1 : ℝ) ≤ 2)
      (by linarith : (-1 : ℝ) ≤ -δ / 2)
    rwa [Real.rpow_neg_one, ← one_div] at this
  have hD : t ^ (1 / 2 : ℝ) * t ^ (-(1 + δ) / 2) = t ^ (-δ / 2) := by
    rw [← Real.rpow_add ht0]; congr 1; ring
  set Z := ‖riemannZeta (1 + ↑t * I)‖
  set E := Real.exp (ε * (0 ^ 2 - t ^ 2))
  have hZ : 0 ≤ Z := norm_nonneg _
  have hE : 0 < E := Real.exp_pos _
  have htd : 0 < t ^ (-δ / 2) := Real.rpow_pos_of_pos ht0 _
  simp only [gf, norm_mul]
  rw [hw0, hg0]
  calc ‖riemannZeta (↑t * I)‖ * t ^ (-(1 + δ) / 2) * E
      ≤ ((2 / 5) * t ^ (1 / 2 : ℝ) * Z) * t ^ (-(1 + δ) / 2) * E := by gcongr
    _ = (2 / 5) * (Z * t ^ (-δ / 2) * E) := by rw [← hD]; ring
    _ ≤ ((4 / 5) * (2 : ℝ) ^ (-δ / 2)) * (Z * t ^ (-δ / 2) * E) :=
        mul_le_mul_of_nonneg_right (by linarith) (by positivity)
    _ = (4 / 5) * Z * (2 * t) ^ (-δ / 2) * E := by rw [hB]; ring
    _ ≤ (4 / 5) * Z * ‖(1 + ↑t * I : ℂ)‖ ^ ((1 - (1 + δ)) / 2) * E := by gcongr
    _ ≤ (4 / 5) * Z * ‖wt δ (1 + ↑t * I)‖ * ‖gs ε (1 + ↑t * I)‖ := by gcongr
    _ = 4 / 5 * (Z * ‖wt δ (1 + ↑t * I)‖ * ‖gs ε (1 + ↑t * I)‖) := by ring

end

end Fayed.Convexity
