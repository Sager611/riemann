import Riemann.Deep.Hardy.Limit

/-!
# Hardy's theorem: ζ has infinitely many zeros on the critical line

Proof (Hardy 1914, Titchmarsh §10.2, with the kernel normalisation of Kernel.lean).
Write L(t) = Λ(1/2 + i t) (real, even, continuous; L(t) = 0 ⇔ ζ(1/2+it) = 0).
For |α| < π/4 and every n (Moments.lean, Limit.lean)

  J (2n) α = ∫ L(t) t^{2n} e^{αt} dt = 2π (-4)^n (Rk^{(2n)}(-2αi) - 16^{-n} · 2cos(α/2)),

and Rk^{(2n)}(-2αi) → 0 as α → π/4⁻. If ζ had only finitely many zeros on the line, then
L would have a constant sign s for |t| > T; for n large the region [T+1, T+2] dominates
[-T, T], so s · J (2n) α ≥ c_n > 0 uniformly in α ∈ [0, π/4). Choosing n with
s (-1)^n > 0 and α close to π/4 makes s · J (2n) α < 0. Contradiction.
-/

open Complex Real MeasureTheory Set Filter Topology
open scoped ComplexConjugate Nat

namespace Fayed.Hardy

/-- Λ on the critical line, as a real function. -/
noncomputable def Lr (t : ℝ) : ℝ := (completedRiemannZeta (1 / 2 + t * I)).re

lemma Lr_eq (t : ℝ) : ((Lr t : ℝ) : ℂ) = completedRiemannZeta (1 / 2 + t * I) :=
  Complex.conj_eq_iff_re.mp (Lam_conj t)

lemma continuous_Lr : Continuous Lr := continuous_re.comp continuous_Lam

lemma Lr_neg (t : ℝ) : Lr (-t) = Lr t := by
  simp only [Lr]; rw [Lam_neg]

lemma moment_re (m : ℕ) (α t : ℝ) :
    (completedRiemannZeta (1 / 2 + t * I) * (t : ℂ) ^ m * (rexp (α * t) : ℂ)).re =
      Lr t * t ^ m * rexp (α * t) := by
  rw [← Lr_eq, ← ofReal_pow, ← ofReal_mul, ← ofReal_mul, ofReal_re]

lemma integrable_Lr_moment {α : ℝ} (hα : |α| < π / 4) (m : ℕ) :
    Integrable (fun t : ℝ => Lr t * t ^ m * rexp (α * t)) :=
  (integrable_moment hα m).re.congr (Eventually.of_forall fun t => moment_re m α t)

lemma re_J {α : ℝ} (hα : |α| < π / 4) (m : ℕ) :
    (J m α).re = ∫ t : ℝ, Lr t * t ^ m * rexp (α * t) := by
  have := integral_re (integrable_moment hα m)
  simp only [RCLike.re_to_complex] at this
  rw [J, ← this]
  congr 1; ext t; exact moment_re m α t

lemma re_J_even (n : ℕ) {α : ℝ} (hα0 : 0 ≤ α) (hα : α < π / 4) :
    (J (2 * n) α).re = 2 * π * (-4) ^ n *
      ((iteratedDeriv (2 * n) Rk (-2 * α * I)).re - (1 / 16) ^ n * (2 * Real.cos (α / 2))) := by
  have hα' : |α| < π / 4 := by rwa [abs_of_nonneg hα0]
  have hmem : (-2 * (α : ℂ) * I) ∈ strip := by
    show |(-2 * (α : ℂ) * I).im| < π / 2
    simp [abs_mul]; linarith
  have h : J (2 * n) α = ((2 * π * (-4) ^ n : ℝ) : ℂ) * (iteratedDeriv (2 * n) Rk (-2 * α * I) -
      (((1 / 16) ^ n * (2 * Real.cos (α / 2)) : ℝ) : ℂ)) := by
    rw [J_eq (2 * n) hα', iteratedDeriv_Hk_eq _ hmem, iteratedDeriv_Ck_even]
    beta_reduce
    rw [Ck_neg_two_mul_I, pow_mul, show (-2 * I) ^ 2 = (-4 : ℂ) by ring_nf; rw [I_sq]; ring]
    push_cast; ring
  rw [h, re_ofReal_mul, sub_re, ofReal_re]

/-- **Hardy's theorem**: ζ has infinitely many zeros on the critical line. -/
theorem hardy : Set.Infinite {t : ℝ | riemannZeta ((1 / 2 : ℂ) + t * Complex.I) = 0} := by
  intro hfin
  have hπ3 := Real.pi_gt_three
  have hπ4 := Real.pi_le_four
  obtain ⟨T₀, hT₀⟩ := (hfin.image (fun t => |t|)).bddAbove
  set T := max T₀ 1 with hT
  have hT1 : 1 ≤ T := le_max_right _ _
  -- no zeros of L beyond T
  have hnz : ∀ t : ℝ, T < |t| → Lr t ≠ 0 := by
    intro t ht h0
    have hΛ : completedRiemannZeta (1 / 2 + t * I) = 0 := by rw [← Lr_eq, h0, ofReal_zero]
    have hs0 : (1 / 2 + t * I : ℂ) ≠ 0 := by
      intro h; have := congrArg Complex.re h; simp at this
    have hζ : riemannZeta (1 / 2 + t * I) = 0 := by
      rw [riemannZeta_def_of_ne_zero hs0, hΛ, zero_div]
    have := hT₀ ⟨t, hζ, rfl⟩
    linarith [le_max_left T₀ 1]
  -- constant sign beyond T
  have hprod : ∀ a b : ℝ, T < a → T < b → 0 < Lr a * Lr b := by
    intro a b ha hb
    by_contra hneg
    push_neg at hneg
    have ha0 := hnz a (by rw [abs_of_pos (by linarith)]; exact ha)
    have hb0 := hnz b (by rw [abs_of_pos (by linarith)]; exact hb)
    have h0 : (0 : ℝ) ∈ uIcc (Lr a) (Lr b) := by
      rw [mem_uIcc]
      rcases lt_or_gt_of_ne ha0 with h | h
      · left; refine ⟨h.le, ?_⟩
        by_contra hb'; push_neg at hb'; nlinarith [mul_pos_of_neg_of_neg h hb']
      · right; refine ⟨?_, h.le⟩
        by_contra hb'; push_neg at hb'; nlinarith [mul_pos h hb']
    obtain ⟨c, hc, hc0⟩ := intermediate_value_uIcc continuous_Lr.continuousOn h0
    have hcT : T < c := lt_of_lt_of_le (lt_min ha hb) hc.1
    exact hnz c (by rw [abs_of_pos (by linarith)]; exact hcT) hc0
  set s := Lr (T + 1) with hs_def
  have hs : ∀ t : ℝ, T < |t| → 0 < s * Lr t := by
    intro t ht
    rcases le_or_gt 0 t with h | h
    · rw [abs_of_nonneg h] at ht; exact hprod _ _ (by linarith) ht
    · rw [abs_of_neg h] at ht; rw [← Lr_neg]; exact hprod _ _ (by linarith) ht
  have hs0 : s ≠ 0 := hnz (T + 1) (by rw [abs_of_pos (by linarith)]; linarith)
  -- compactness constants
  obtain ⟨t₀, ht₀, hmin⟩ := isCompact_Icc.exists_isMinOn
    (nonempty_Icc.mpr (by linarith : T + 1 ≤ T + 2))
    ((continuous_Lr.const_mul s).continuousOn (s := Icc (T + 1) (T + 2)))
  set m₀ := s * Lr t₀ with hm₀_def
  have hm₀ : 0 < m₀ := hs t₀ (by rw [abs_of_pos (by linarith [ht₀.1])]; linarith [ht₀.1])
  obtain ⟨Mx, hMx⟩ := isCompact_Icc.exists_bound_of_continuousOn
    (continuous_Lr.continuousOn (s := Icc (-T) T))
  have hMx0 : 0 ≤ Mx := (norm_nonneg _).trans (hMx 0 ⟨by linarith, by linarith⟩)
  set A := |s| * Mx * rexp T * (2 * T) with hA
  -- uniform lower bound
  have hlow : ∀ n : ℕ, ∀ α : ℝ, 0 ≤ α → α < π / 4 →
      m₀ * (T + 1) ^ (2 * n) - A * T ^ (2 * n) ≤ s * (J (2 * n) α).re := by
    intro n α hα0 hα
    have hα' : |α| < π / 4 := by rwa [abs_of_nonneg hα0]
    rw [re_J hα', ← integral_const_mul]
    set c₁ := m₀ * (T + 1) ^ (2 * n)
    set c₂ := |s| * Mx * T ^ (2 * n) * rexp T
    set g : ℝ → ℝ := fun t => (Icc (T + 1) (T + 2)).indicator (fun _ => c₁) t -
      (Icc (-T) T).indicator (fun _ => c₂) t with hg
    have hg1 : IntegrableOn (fun _ : ℝ => c₁) (Icc (T + 1) (T + 2)) :=
      integrableOn_const (by simp [Real.volume_Icc])
    have hg2 : IntegrableOn (fun _ : ℝ => c₂) (Icc (-T) T) :=
      integrableOn_const (by simp [Real.volume_Icc])
    have hgint : Integrable g :=
      (hg1.integrable_indicator measurableSet_Icc).sub (hg2.integrable_indicator measurableSet_Icc)
    have hgval : ∫ t, g t = c₁ - c₂ * (2 * T) := by
      rw [hg, integral_sub (hg1.integrable_indicator measurableSet_Icc)
        (hg2.integrable_indicator measurableSet_Icc), integral_indicator_const _ measurableSet_Icc,
        integral_indicator_const _ measurableSet_Icc]
      simp only [smul_eq_mul, Measure.real, Real.volume_Icc]
      rw [ENNReal.toReal_ofReal (by linarith), ENNReal.toReal_ofReal (by linarith)]
      ring
    have hle : ∀ t, g t ≤ s * (Lr t * t ^ (2 * n) * rexp (α * t)) := by
      intro t
      have hev : 0 ≤ t ^ (2 * n) := (even_two_mul n).pow_nonneg t
      by_cases h1 : t ∈ Icc (T + 1) (T + 2)
      · have h2 : t ∉ Icc (-T) T := fun h => by linarith [h.2, h1.1]
        simp only [hg, indicator_of_mem h1, indicator_of_notMem h2, sub_zero]
        have hsL : m₀ ≤ s * Lr t := hmin h1
        have hp : (T + 1) ^ (2 * n) ≤ t ^ (2 * n) :=
          pow_le_pow_left₀ (by linarith) h1.1 _
        have he : 1 ≤ rexp (α * t) := Real.one_le_exp (by nlinarith [h1.1])
        calc c₁ = m₀ * (T + 1) ^ (2 * n) * 1 := by ring
          _ ≤ (s * Lr t) * t ^ (2 * n) * rexp (α * t) :=
            mul_le_mul (mul_le_mul hsL hp (by positivity) (by linarith)) he zero_le_one
              (mul_nonneg (by linarith) hev)
          _ = _ := by ring
      · by_cases h2 : t ∈ Icc (-T) T
        · simp only [hg, indicator_of_notMem h1, indicator_of_mem h2, zero_sub]
          have habs : |s * (Lr t * t ^ (2 * n) * rexp (α * t))| ≤ c₂ := by
            rw [abs_mul, abs_mul, abs_mul, abs_pow, abs_of_pos (Real.exp_pos _)]
            have ht : |t| ≤ T := abs_le.mpr ⟨h2.1, h2.2⟩
            have hL : |Lr t| ≤ Mx := by simpa using hMx t h2
            have he : rexp (α * t) ≤ rexp T := Real.exp_le_exp.mpr (by
              have : α * t ≤ |α * t| := le_abs_self _
              rw [abs_mul, abs_of_nonneg hα0] at this
              nlinarith [abs_nonneg t])
            have hp : |t| ^ (2 * n) ≤ T ^ (2 * n) := pow_le_pow_left₀ (abs_nonneg _) ht _
            calc |s| * (|Lr t| * |t| ^ (2 * n) * rexp (α * t))
                ≤ |s| * (Mx * T ^ (2 * n) * rexp T) := by gcongr
              _ = c₂ := by ring
          linarith [neg_abs_le (s * (Lr t * t ^ (2 * n) * rexp (α * t)))]
        · simp only [hg, indicator_of_notMem h1, indicator_of_notMem h2, sub_zero]
          have ht : T < |t| := by
            by_contra hc; push_neg at hc
            exact h2 ⟨(abs_le.mp hc).1, (abs_le.mp hc).2⟩
          have := hs t ht
          have : 0 ≤ (s * Lr t) * t ^ (2 * n) * rexp (α * t) :=
            mul_nonneg (mul_nonneg this.le hev) (Real.exp_pos _).le
          linarith
    have hfint := (integrable_Lr_moment hα' (2 * n)).const_mul s
    calc m₀ * (T + 1) ^ (2 * n) - A * T ^ (2 * n) = c₁ - c₂ * (2 * T) := by
          simp only [c₁, c₂, hA]; ring
      _ = ∫ t, g t := hgval.symm
      _ ≤ _ := integral_mono hgint hfint hle
  -- for large n the lower bound is positive
  have hr : 1 < ((T + 1) / T) ^ 2 := by
    have : 1 < (T + 1) / T := by rw [lt_div_iff₀ (by linarith)]; linarith
    nlinarith
  have hev := (tendsto_pow_atTop_atTop_of_one_lt hr).eventually (eventually_gt_atTop (A / m₀))
  obtain ⟨N, hN⟩ := eventually_atTop.mp hev
  have hpos : ∀ n ≥ N, 0 < m₀ * (T + 1) ^ (2 * n) - A * T ^ (2 * n) := by
    intro n hn
    have h1 := hN n hn
    have hT0 : 0 < T ^ (2 * n) := by positivity
    have h2 : (T + 1) ^ (2 * n) = (((T + 1) / T) ^ 2) ^ n * T ^ (2 * n) := by
      rw [← pow_mul, ← mul_pow, div_mul_cancel₀ _ (by linarith)]
    rw [h2]
    rw [div_lt_iff₀ hm₀] at h1
    nlinarith
  -- choose the parity of n against the sign s
  obtain ⟨n, hn, hsn⟩ : ∃ n ≥ N, 0 < s * (-4 : ℝ) ^ n := by
    rcases lt_or_gt_of_ne hs0 with h | h
    · refine ⟨2 * N + 1, by omega, ?_⟩
      rw [pow_succ, pow_mul]; norm_num; nlinarith [pow_pos (by norm_num : (0 : ℝ) < 16) N]
    · refine ⟨2 * N, by omega, ?_⟩
      rw [pow_mul]; norm_num; positivity
  -- choose α close to π/4
  have hlim := Metric.tendsto_nhds.mp (tendsto_iteratedDeriv_Rk (2 * n)) ((1 / 16) ^ n)
    (by positivity)
  obtain ⟨α, hαD, hαI⟩ := (hlim.and (Ioo_mem_nhdsLT (by positivity : (0 : ℝ) < π / 4))).exists
  rw [dist_zero_right] at hαD
  have hcos : 1 ≤ 2 * Real.cos (α / 2) := by
    have := Real.one_sub_sq_div_two_le_cos (x := α / 2)
    nlinarith [hαI.1, hαI.2]
  have hneg : s * (J (2 * n) α).re < 0 := by
    rw [re_J_even n hαI.1.le hαI.2]
    have hD : (iteratedDeriv (2 * n) Rk (-2 * α * I)).re - (1 / 16) ^ n * (2 * Real.cos (α / 2))
        < 0 := by
      have h1 := (re_le_norm _).trans_lt hαD
      have h2 : (1 / 16 : ℝ) ^ n ≤ (1 / 16) ^ n * (2 * Real.cos (α / 2)) :=
        le_mul_of_one_le_right (by positivity) hcos
      linarith
    have : s * (2 * π * (-4) ^ n * ((iteratedDeriv (2 * n) Rk (-2 * α * I)).re -
        (1 / 16) ^ n * (2 * Real.cos (α / 2)))) = (2 * π) * (s * (-4) ^ n) *
        ((iteratedDeriv (2 * n) Rk (-2 * α * I)).re - (1 / 16) ^ n * (2 * Real.cos (α / 2))) := by
      ring
    rw [this]
    exact mul_neg_of_pos_of_neg (mul_pos (by positivity) hsn) hD
  have := (hpos n hn).trans_le (hlow n α hαI.1.le hαI.2)
  linarith

end Fayed.Hardy

#print axioms Fayed.Hardy.hardy
