import Riemann.Deep.Convexity.Bootstrap

/-!
# Max-modulus step, iteration, and the convexity bound

* step: if |g| ≤ M on H then |g| ≤ K + (4/5)M on H (rectangles [0,1+δ]×[1,T], T large).
* gf_bound: |g| ≤ 5K on H, K = K δ independent of ε.
* convexity_holds: the exact statement of Fayed.Inputs.convexity.
-/

open Complex Filter Topology Set MeasureTheory Asymptotics
open scoped Real ComplexConjugate

namespace Fayed.Convexity

noncomputable section

/-- Bound for the right edge and bottom edge (independent of ε). -/
def Kc (δ : ℝ) : ℝ :=
  max ((∑' n : ℕ, 1 / ((n : ℝ) + 1) ^ (1 + δ)) * Real.exp 1 * Real.exp 4)
    (max (24 * Real.exp 3 * Real.exp 1 * Real.exp 4) 1)

lemma one_le_Kc (δ : ℝ) : 1 ≤ Kc δ := le_trans (le_max_right _ _) (le_max_right _ _)

lemma iterate_bound {X : Type*} (P : X → Prop) (f : X → ℝ) {K θ M0 : ℝ} (hK : 0 ≤ K)
    (hθ0 : 0 ≤ θ) (hθ1 : θ < 1) (h0 : ∀ x, P x → f x ≤ M0)
    (hstep : ∀ M, (∀ x, P x → f x ≤ M) → ∀ x, P x → f x ≤ K + θ * M) :
    ∀ x, P x → f x ≤ K / (1 - θ) := by
  set L := K / (1 - θ) with hL_def
  have h1θ : (1 - θ) ≠ 0 := by linarith
  have hL : K + θ * L = L := by rw [hL_def]; field_simp; ring
  have hL0 : 0 ≤ L := div_nonneg hK (by linarith)
  have hn : ∀ n : ℕ, ∀ x, P x → f x ≤ L + θ ^ n * |M0| := by
    intro n
    induction n with
    | zero =>
      intro x hx
      have := h0 x hx
      simp only [pow_zero, one_mul]
      linarith [le_abs_self M0]
    | succ n ih =>
      intro x hx
      have := hstep _ ih x hx
      calc f x ≤ K + θ * (L + θ ^ n * |M0|) := this
        _ = L + θ ^ (n + 1) * |M0| := by rw [pow_succ]; linear_combination hL
  intro x hx
  have ht : Tendsto (fun n : ℕ => L + θ ^ n * |M0|) atTop (𝓝 (L + 0 * |M0|)) :=
    tendsto_const_nhds.add ((tendsto_pow_atTop_nhds_zero_of_lt_one hθ0 hθ1).mul_const _)
  rw [zero_mul, add_zero] at ht
  exact ge_of_tendsto' ht (fun n => hn n x hx)

lemma step {δ ε : ℝ} (hδ0 : 0 < δ) (hδ1 : δ ≤ 1) (hε0 : 0 < ε) (hε1 : ε ≤ 1) (M : ℝ)
    (hM : ∀ s ∈ H δ, ‖gf δ ε s‖ ≤ M) : ∀ s ∈ H δ, ‖gf δ ε s‖ ≤ Kc δ + (4 / 5) * M := by
  intro z hz
  have hM0 : 0 ≤ M := (norm_nonneg _).trans (hM z hz)
  have hK1 := one_le_Kc δ
  set L := 5 + Real.log 24 with hL_def
  have hL : 0 ≤ L := by have := Real.log_nonneg (by norm_num : (1 : ℝ) ≤ 24); linarith
  set T := max (z.im + 1) ((3 + L) / ε) with hT_def
  have hT1 : z.im + 1 ≤ T := le_max_left _ _
  have hT2 : (3 + L) / ε ≤ T := le_max_right _ _
  have hz1 := hz.2.2
  have hT : 1 < T := by linarith
  have hεT : 3 + L ≤ ε * T := by rwa [div_le_iff₀ hε0, mul_comm] at hT2
  -- top edge
  have htop : ∀ u : ℂ, u ∈ H δ → u.im = T → ‖gf δ ε u‖ ≤ 1 := by
    intro u hu hut
    refine (gf_gen hδ1 hε0 hε1 hu).trans ?_
    rw [hut]
    have hq : 5 + 3 * T - ε * T ^ 2 ≤ -Real.log 24 := by nlinarith
    calc 24 * Real.exp (5 + 3 * T - ε * T ^ 2) ≤ 24 * Real.exp (-Real.log 24) := by gcongr
      _ = 1 := by rw [Real.exp_neg, Real.exp_log (by norm_num)]; norm_num
  set R := Ioo (0 : ℝ) (1 + δ) ×ℂ Ioo (1 : ℝ) T with hR_def
  have hRb : Bornology.IsBounded R :=
    (Metric.isBounded_Ioo _ _).reProdIm (Metric.isBounded_Ioo _ _)
  have hcl : closure R = Icc 0 (1 + δ) ×ℂ Icc 1 T := by
    rw [hR_def, closure_reProdIm, closure_Ioo (by linarith), closure_Ioo (by linarith)]
  have hd : DiffContOnCl ℂ (gf δ ε) R := by
    apply DifferentiableOn.diffContOnCl
    rw [hcl]; intro u hu
    have hu2 : 1 ≤ u.im := hu.2.1
    exact (differentiableAt_gf δ ε (by linarith)).differentiableWithinAt
  have hfr : ∀ u ∈ frontier R, ‖gf δ ε u‖ ≤ Kc δ + 4 / 5 * M := by
    intro u hu
    rw [hR_def, frontier_reProdIm, closure_Ioo (by linarith), closure_Ioo (by linarith),
      frontier_Ioo (by linarith), frontier_Ioo (by linarith)] at hu
    rcases hu with ⟨hre, him⟩ | ⟨hre, him⟩
    · have hre' : u.re ∈ Icc 0 (1 + δ) := hre
      have him' : u.im = 1 ∨ u.im = T := him
      rcases him' with h | h
      · -- bottom
        have huH : u ∈ H δ := ⟨hre'.1, hre'.2, by rw [h]⟩
        refine (gf_le_H hδ1 hε0 hε1 huH).trans ?_
        have hzb := apriori hre'.1 (by linarith [hre'.2]) (by rw [h])
        rw [h] at hzb
        have : ‖riemannZeta u‖ * Real.exp 1 * Real.exp (4 - ε * u.im ^ 2) ≤
            24 * Real.exp (3 * 1) * Real.exp 1 * Real.exp 4 := by
          gcongr
          nlinarith
        have hk : 24 * Real.exp (3 * 1) * Real.exp 1 * Real.exp 4 ≤ Kc δ := by
          rw [mul_one]; exact le_trans (le_max_left _ _) (le_max_right _ _)
        nlinarith
      · -- top
        have huH : u ∈ H δ := ⟨hre'.1, hre'.2, by rw [h]; linarith⟩
        have := htop u huH h
        nlinarith
    · have him' : u.im ∈ Icc 1 T := him
      have hre' : u.re = 0 ∨ u.re = 1 + δ := hre
      rcases hre' with h | h
      · -- left edge
        have hu : u = (u.im : ℂ) * I := by apply Complex.ext <;> simp [h]
        rw [hu]
        refine (gf_left hδ0 hδ1 hε0.le him'.1).trans ?_
        have hH : (1 + (u.im : ℂ) * I) ∈ H δ := ⟨by simp, by simp; linarith, by simpa using him'.1⟩
        have := hM _ hH
        nlinarith
      · -- right edge
        have huH : u ∈ H δ := ⟨by rw [h]; linarith, by rw [h], him'.1⟩
        refine (gf_le_H hδ1 hε0 hε1 huH).trans ?_
        have hzr := zeta_right_bound hδ0 h
        have : ‖riemannZeta u‖ * Real.exp 1 * Real.exp (4 - ε * u.im ^ 2) ≤
            (∑' n : ℕ, 1 / ((n : ℝ) + 1) ^ (1 + δ)) * Real.exp 1 * Real.exp 4 := by
          gcongr
          nlinarith
        have hk : (∑' n : ℕ, 1 / ((n : ℝ) + 1) ^ (1 + δ)) * Real.exp 1 * Real.exp 4 ≤ Kc δ :=
          le_max_left _ _
        nlinarith
  exact Complex.norm_le_of_forall_mem_frontier_norm_le hRb hd hfr
    (by rw [hcl]; exact ⟨⟨hz.1, hz.2.1⟩, ⟨hz1, by linarith⟩⟩)

/-- Uniform bound on the half-strip, independent of ε ∈ (0,1]. -/
lemma gf_bound {δ ε : ℝ} (hδ0 : 0 < δ) (hδ1 : δ ≤ 1) (hε0 : 0 < ε) (hε1 : ε ≤ 1) :
    ∀ s ∈ H δ, ‖gf δ ε s‖ ≤ 5 * Kc δ := by
  have h0 : ∀ s ∈ H δ, ‖gf δ ε s‖ ≤ 24 * Real.exp (5 + 9 / (4 * ε)) := by
    intro s hs
    refine (gf_gen hδ1 hε0 hε1 hs).trans ?_
    gcongr
    have : 3 * s.im - ε * s.im ^ 2 ≤ 9 / (4 * ε) := by
      rw [le_div_iff₀ (by positivity)]; nlinarith [sq_nonneg (2 * ε * s.im - 3)]
    linarith
  have := iterate_bound (fun s => s ∈ H δ) (fun s => ‖gf δ ε s‖) (K := Kc δ) (θ := 4 / 5)
    (by linarith [one_le_Kc δ]) (by norm_num) (by norm_num) h0
    (fun M hM => step hδ0 hδ1 hε0 hε1 M hM)
  intro s hs
  have h := this s hs
  have e : Kc δ / (1 - 4 / 5) = 5 * Kc δ := by ring
  linarith

/-- The convexity bound on Re s = 1-η, for Im s ≥ 1. -/
lemma convexity_pos {η ε : ℝ} (hη0 : 0 < η) (hη1 : η < 1 / 2) (hε : 0 < ε) {t : ℝ}
    (ht : 1 ≤ t) :
    ‖riemannZeta (((1 - η : ℝ) : ℂ) + t * I)‖ ≤
      10 * Real.exp 1 * Kc (min 1 (2 * ε)) * t ^ (η / 2 + ε) := by
  set δ := min 1 (2 * ε) with hδ_def
  have hδ0 : 0 < δ := lt_min one_pos (by linarith)
  have hδ1 : δ ≤ 1 := min_le_left _ _
  have hδ2 : δ ≤ 2 * ε := min_le_right _ _
  have ht0 : 0 < t := by linarith
  set s : ℂ := ((1 - η : ℝ) : ℂ) + t * I with hs_def
  have hre : s.re = 1 - η := by simp [s]
  have him : s.im = t := by simp [s]
  have hsH : s ∈ H δ := ⟨by rw [hre]; linarith, by rw [hre]; linarith, by rw [him]; exact ht⟩
  set ε' : ℝ := 1 / t ^ 2 with hε'_def
  have hε'0 : 0 < ε' := by positivity
  have hε'1 : ε' ≤ 1 := by
    rw [hε'_def, div_le_one (by positivity)]; nlinarith
  have hb := gf_bound hδ0 hδ1 hε'0 hε'1 s hsH
  have hw := wt_ge δ (s := s) (by rw [hre]; linarith) (by rw [him]; exact ht0)
  rw [hre, show (1 - η - (1 + δ)) / 2 = -((η + δ) / 2) by ring] at hw
  have hg : Real.exp (-1) ≤ ‖gs ε' s‖ := by
    rw [norm_gs, hre, him]
    apply Real.exp_le_exp.mpr
    have : ε' * t ^ 2 = 1 := by rw [hε'_def]; field_simp
    nlinarith [sq_nonneg (1 - η), hε'0]
  have hn1 : 1 ≤ ‖s‖ := one_le_norm_of_im (by rw [him]; exact ht)
  have hn2 : ‖s‖ ≤ 2 * t := by
    refine (norm_le_abs_re_add_abs_im _).trans ?_
    rw [hre, him, abs_of_pos (by linarith), abs_of_pos ht0]; linarith
  set P := ‖s‖ ^ ((η + δ) / 2) with hP_def
  have hP0 : 0 < P := Real.rpow_pos_of_pos (by linarith) _
  have hPinv : ‖s‖ ^ (-((η + δ) / 2)) = P⁻¹ := Real.rpow_neg (by linarith) _
  rw [hPinv] at hw
  -- ‖ζ s‖ * P⁻¹ * e^{-1} ≤ 5K
  have hmain : ‖riemannZeta s‖ * P⁻¹ * Real.exp (-1) ≤ 5 * Kc δ := by
    refine le_trans ?_ hb
    simp only [gf, norm_mul]
    gcongr
  have hζ : ‖riemannZeta s‖ ≤ 5 * Kc δ * P * Real.exp 1 := by
    have he : Real.exp (-1) * Real.exp 1 = 1 := by rw [← Real.exp_add]; simp
    have := mul_le_mul_of_nonneg_right hmain (by positivity : 0 ≤ P * Real.exp 1)
    calc ‖riemannZeta s‖ = ‖riemannZeta s‖ * P⁻¹ * Real.exp (-1) * (P * Real.exp 1) := by
          field_simp; rw [mul_assoc, he, mul_one]
      _ ≤ _ := this
      _ = 5 * Kc δ * P * Real.exp 1 := by ring
  -- P ≤ 2 t^{η/2+ε}
  have hP : P ≤ 2 * t ^ (η / 2 + ε) := by
    have h1 : P ≤ (2 * t) ^ ((η + δ) / 2) :=
      Real.rpow_le_rpow (by linarith) hn2 (by linarith)
    rw [Real.mul_rpow (by norm_num) ht0.le] at h1
    have h2 : (2 : ℝ) ^ ((η + δ) / 2) ≤ 2 := by
      have := Real.rpow_le_rpow_of_exponent_le (by norm_num : (1 : ℝ) ≤ 2)
        (by linarith : (η + δ) / 2 ≤ 1)
      rwa [Real.rpow_one] at this
    have h3 : t ^ ((η + δ) / 2) ≤ t ^ (η / 2 + ε) :=
      Real.rpow_le_rpow_of_exponent_le ht (by linarith)
    have h4 : 0 ≤ t ^ ((η + δ) / 2) := Real.rpow_nonneg ht0.le _
    nlinarith
  have hK0 : 0 ≤ Kc δ := by linarith [one_le_Kc δ]
  calc ‖riemannZeta s‖ ≤ 5 * Kc δ * P * Real.exp 1 := hζ
    _ ≤ 5 * Kc δ * (2 * t ^ (η / 2 + ε)) * Real.exp 1 := by gcongr
    _ = 10 * Real.exp 1 * Kc δ * t ^ (η / 2 + ε) := by ring

end

end Fayed.Convexity
