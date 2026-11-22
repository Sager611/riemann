module

public import Riemann.Deep.LeanLi.Lit
public import Mathlib.Analysis.Calculus.Deriv.Star
public import Mathlib.Analysis.Complex.CauchyIntegral

/-!
# LeanLi: basic facts on the nontrivial zeros (strip, symmetries, multiplicity, finiteness)
-/

@[expose] public section

open Complex Filter Topology Set ComplexConjugate

namespace LeanLi

lemma mem_nontrivialZeros {s : ℂ} :
    s ∈ NontrivialZeros ↔ riemannZeta s = 0 ∧ (¬∃ n : ℕ, s = -2 * (n + 1)) ∧ s ≠ 1 := Iff.rfl

/-- Nontrivial zeros have real part `< 1` (Mathlib: no zeros on `Re s ≥ 1`). -/
lemma re_lt_one_of_mem {s : ℂ} (hs : s ∈ NontrivialZeros) : s.re < 1 := by
  by_contra h
  exact riemannZeta_ne_zero_of_one_le_re (not_lt.1 h) hs.1

/-- Nontrivial zeros have real part `> 0` (functional equation). -/
lemma re_pos_of_mem {s : ℂ} (hs : s ∈ NontrivialZeros) : 0 < s.re := by
  obtain ⟨hz, htriv, _⟩ := hs
  by_contra h
  push Not at h
  have hs0 : s ≠ 0 := by
    rintro rfl; rw [riemannZeta_zero] at hz; norm_num at hz
  set t := 1 - s with ht
  have htre : 1 ≤ t.re := by simp only [ht, sub_re, one_re]; linarith
  have h1 : ∀ n : ℕ, t ≠ -n := by
    intro n hn
    have := congrArg Complex.re hn
    simp at this
    have : (0 : ℝ) ≤ n := Nat.cast_nonneg n
    linarith
  have h2 : t ≠ 1 := by intro h; apply hs0; linear_combination -h
  have hfe := riemannZeta_one_sub h1 h2
  rw [show 1 - t = s by ring, hz] at hfe
  have hζ : riemannZeta t ≠ 0 := riemannZeta_ne_zero_of_one_le_re htre
  have hΓ : Gamma t ≠ 0 := Complex.Gamma_ne_zero h1
  have hp : (2 * (Real.pi : ℂ)) ^ (-t) ≠ 0 := by
    rw [Ne, cpow_eq_zero_iff]; push Not; intro h; exfalso
    have : (2 : ℂ) * Real.pi ≠ 0 := by
      exact mul_ne_zero two_ne_zero (ofReal_ne_zero.2 Real.pi_ne_zero)
    exact this h
  have hcos : Complex.cos (Real.pi * t / 2) = 0 := by
    have := hfe.symm
    simp only [mul_eq_zero, OfNat.ofNat_ne_zero, false_or, hp, hΓ, hζ, or_false] at this
    exact this
  obtain ⟨k, hk⟩ := Complex.cos_eq_zero_iff.1 hcos
  have hpi : (Real.pi : ℂ) ≠ 0 := ofReal_ne_zero.2 Real.pi_ne_zero
  have htk : t = 2 * k + 1 := by
    have := hk
    field_simp at this
    linear_combination this
  have hsk : s = -2 * (k : ℂ) := by linear_combination -htk + ht
  -- s.re ≤ 0, s ≠ 0 ⇒ k ≥ 1
  have hkre : s.re = -2 * k := by rw [hsk]; simp
  have hk1 : (1 : ℤ) ≤ k := by
    have hk0 : (0 : ℝ) ≤ k := by linarith
    have hkne : k ≠ 0 := by rintro rfl; apply hs0; rw [hsk]; simp
    have : (0 : ℤ) ≤ k := by exact_mod_cast hk0
    omega
  apply htriv
  refine ⟨(k - 1).toNat, ?_⟩
  rw [hsk]
  have : ((k - 1).toNat : ℤ) = k - 1 := Int.toNat_of_nonneg (by omega)
  have hc : (((k - 1).toNat : ℕ) : ℂ) = (k : ℂ) - 1 := by
    have := congrArg (fun z : ℤ => (z : ℂ)) this
    push_cast at this
    exact_mod_cast this
  rw [hc]; ring

/-- `ζ(s̄) = conj ζ(s)` (identity theorem on `ℂ ∖ {1}`, from real values on `(1, ∞)`). -/
lemma riemannZeta_conj (s : ℂ) : riemannZeta (conj s) = conj (riemannZeta s) := by
  by_cases hs1 : s = 1
  · subst hs1
    have h := riemannZeta_one
    have h4 : (4 * (Real.pi : ℂ)) = ((4 * Real.pi : ℝ) : ℂ) := by push_cast; ring
    have hl : (Complex.log (4 * (Real.pi : ℂ))).im = 0 := by
      rw [h4, Complex.log_im, Complex.arg_ofReal_of_nonneg (by positivity)]
    have him : (riemannZeta 1).im = 0 := by rw [h]; simp [hl]
    simp only [map_one]
    exact (Complex.conj_eq_iff_im.2 him).symm
  set U : Set ℂ := {1}ᶜ
  have hU : IsPreconnected U := (isConnected_compl_singleton_of_one_lt_rank (by simp) 1).isPreconnected
  have hf : AnalyticOnNhd ℂ riemannZeta U := analyticOn_riemannZeta
  have hg : AnalyticOnNhd ℂ (fun z => conj (riemannZeta (conj z))) U := by
    apply DifferentiableOn.analyticOnNhd _ isOpen_compl_singleton
    intro z hz
    have hz' : conj z ≠ 1 := by
      intro h; apply hz; rw [mem_singleton_iff]; rw [← conj_conj z, h, map_one]
    have := (differentiableAt_riemannZeta hz').conj_conj
    rw [conj_conj] at this
    exact this.differentiableWithinAt
  have h2 : (2 : ℂ) ∈ U := by simp only [U, mem_compl_iff, mem_singleton_iff]; norm_num
  have hfreq : ∃ᶠ z in 𝓝[≠] (2 : ℂ), riemannZeta z = conj (riemannZeta (conj z)) := by
    have ht : Tendsto (fun k : ℕ => ((2 + 1 / ((k : ℝ) + 1) : ℝ) : ℂ)) atTop (𝓝[≠] 2) := by
      rw [tendsto_nhdsWithin_iff]
      constructor
      · have : Tendsto (fun k : ℕ => (2 + 1 / ((k : ℝ) + 1) : ℝ)) atTop (𝓝 (2 + 0)) :=
          tendsto_const_nhds.add tendsto_one_div_add_atTop_nhds_zero_nat
        rw [add_zero] at this
        have h := (Complex.continuous_ofReal.tendsto 2).comp this
        rw [Complex.ofReal_ofNat] at h
        exact h
      · refine Eventually.of_forall fun k => ?_
        simp only [mem_compl_iff, mem_singleton_iff]
        intro h
        have h' : (2 + 1 / ((k : ℝ) + 1) : ℝ) = 2 := by exact_mod_cast h
        have : (0 : ℝ) < 1 / ((k : ℝ) + 1) := by positivity
        linarith
    refine ht.frequently (Frequently.of_forall fun k => ?_)
    set x : ℝ := 2 + 1 / ((k : ℝ) + 1)
    have hx : 1 < x := by have : (0 : ℝ) < 1 / ((k : ℝ) + 1) := by positivity
                          linarith
    rw [Complex.conj_ofReal]
    exact (Complex.conj_eq_iff_im.2 (riemannZeta_im_eq_zero_of_one_lt hx)).symm
  have := hf.eqOn_of_preconnected_of_frequently_eq hg hU h2 hfreq
  have hs' : conj s ∈ U := by
    simp only [U, mem_compl_iff, mem_singleton_iff]
    intro h; apply hs1; rw [← conj_conj s, h, map_one]
  have h := this hs'
  simp only [conj_conj] at h
  exact h

lemma conj_mem {s : ℂ} (hs : s ∈ NontrivialZeros) : conj s ∈ NontrivialZeros := by
  obtain ⟨hz, ht, h1⟩ := hs
  refine ⟨by rw [riemannZeta_conj, hz, map_zero], ?_, ?_⟩
  · rintro ⟨n, hn⟩
    apply ht
    refine ⟨n, ?_⟩
    have := congrArg conj hn
    simp only [conj_conj, map_neg, map_mul, map_add, map_one, map_natCast, map_ofNat] at this
    exact this
  · intro h; apply h1; rw [← conj_conj s, h, map_one]

/-- The analytic order is invariant under `f ↦ conj ∘ f ∘ conj`. -/
lemma analyticOrderAt_conj_conj {f : ℂ → ℂ} {z₀ : ℂ} (hf : AnalyticAt ℂ f z₀) :
    analyticOrderAt (fun z => conj (f (conj z))) (conj z₀) = analyticOrderAt f z₀ := by
  have hcont : Tendsto (fun z : ℂ => conj z) (𝓝 (conj z₀)) (𝓝 z₀) := by
    have := (Complex.continuous_conj.tendsto (conj z₀)); rwa [conj_conj] at this
  have hF : AnalyticAt ℂ (fun z => conj (f (conj z))) (conj z₀) := by
    rw [Complex.analyticAt_iff_eventually_differentiableAt] at hf ⊢
    filter_upwards [hcont.eventually hf] with z hz
    have := hz.conj_conj
    rw [conj_conj] at this
    exact this
  cases h : analyticOrderAt f z₀ with
  | top =>
    rw [analyticOrderAt_eq_top] at h ⊢
    filter_upwards [hcont.eventually h] with z hz
    simp [hz]
  | coe n =>
    obtain ⟨g, hg, hg0, hfg⟩ := (hf.analyticOrderAt_eq_natCast).1 h
    refine (hF.analyticOrderAt_eq_natCast).2 ⟨fun z => conj (g (conj z)), ?_, ?_, ?_⟩
    · rw [Complex.analyticAt_iff_eventually_differentiableAt] at hg ⊢
      filter_upwards [hcont.eventually hg] with z hz
      have := hz.conj_conj
      rw [conj_conj] at this
      exact this
    · simpa using hg0
    · filter_upwards [hcont.eventually hfg] with z hz
      rw [hz]
      simp [smul_eq_mul, map_mul, map_pow, map_sub, conj_conj]

lemma analyticAt_riemannZeta {s : ℂ} (hs : s ≠ 1) : AnalyticAt ℂ riemannZeta s :=
  analyticOn_riemannZeta s hs

lemma mult_conj {s : ℂ} (hs : s ≠ 1) : mult (conj s) = mult s := by
  unfold mult analyticOrderNatAt
  have hs' : conj s ≠ 1 := by intro h; apply hs; rw [← conj_conj s, h, map_one]
  have hev : (fun z => conj (riemannZeta (conj z))) =ᶠ[𝓝 (conj s)] riemannZeta := by
    filter_upwards [isOpen_compl_singleton.mem_nhds hs'] with z _
    rw [riemannZeta_conj, conj_conj]
  rw [← analyticOrderAt_congr hev, analyticOrderAt_conj_conj (analyticAt_riemannZeta hs)]

/-- ζ is not locally zero anywhere on `ℂ ∖ {1}`. -/
lemma analyticOrderAt_ne_top {s : ℂ} (hs : s ≠ 1) : analyticOrderAt riemannZeta s ≠ ⊤ := by
  have hU : IsPreconnected ({1}ᶜ : Set ℂ) :=
    (isConnected_compl_singleton_of_one_lt_rank (by simp) 1).isPreconnected
  have h2 : (2 : ℂ) ∈ ({1}ᶜ : Set ℂ) := by rw [mem_compl_iff, mem_singleton_iff]; norm_num
  refine analyticOn_riemannZeta.analyticOrderAt_ne_top_of_isPreconnected hU h2 hs ?_
  rw [Ne, analyticOrderAt_eq_top]
  intro h
  have := h.self_of_nhds
  exact riemannZeta_ne_zero_of_one_lt_re (by norm_num) this

lemma mult_pos {s : ℂ} (hs : s ∈ NontrivialZeros) : 0 < mult s := by
  unfold mult analyticOrderNatAt
  have h0 : analyticOrderAt riemannZeta s ≠ 0 :=
    analyticOrderAt_ne_zero.2 ⟨analyticAt_riemannZeta hs.2.2, hs.1⟩
  have ht := analyticOrderAt_ne_top hs.2.2
  revert h0 ht
  cases analyticOrderAt riemannZeta s with
  | top => simp
  | coe n => intro h _; simp at h ⊢; omega

lemma one_sub_mem {s : ℂ} (hs : s ∈ NontrivialZeros) : 1 - s ∈ NontrivialZeros := by
  have hre := re_pos_of_mem hs
  have hre1 := re_lt_one_of_mem hs
  have h1 : ∀ n : ℕ, s ≠ -n := by
    intro n hn
    have := congrArg Complex.re hn; simp at this
    have : (0 : ℝ) ≤ n := Nat.cast_nonneg n
    linarith
  refine ⟨?_, ?_, ?_⟩
  · rw [riemannZeta_one_sub h1 hs.2.2, hs.1, mul_zero]
  · rintro ⟨n, hn⟩
    have := congrArg Complex.re hn; simp at this
    have : (0 : ℝ) ≤ n := Nat.cast_nonneg n
    linarith
  · intro h
    have := congrArg Complex.re h; simp at this
    linarith

/-- No nontrivial zero is real (input: `ζ ≠ 0` on `(0, 1)`). -/
lemma im_ne_zero (h : ZetaNoZerosUnitInterval) {s : ℂ} (hs : s ∈ NontrivialZeros) : s.im ≠ 0 := by
  intro him
  have hs' : s = (s.re : ℂ) := by
    apply Complex.ext <;> simp [him]
  apply h s.re (re_pos_of_mem hs) (re_lt_one_of_mem hs)
  rw [← hs']
  exact hs.1

/-- The zeros with `|ρ| ≤ T` form a finite set. -/
lemma finite_norm_le (T : ℝ) : {ρ : ℂ | ρ ∈ NontrivialZeros ∧ ‖ρ‖ ≤ T}.Finite := by
  refine ((isCompact_closedBall (0 : ℂ) T).inter_riemannZetaZeros_finite).subset ?_
  intro ρ hρ
  exact ⟨by simpa using hρ.2, hρ.1.1⟩

lemma norm_le_of_mem_zerosUpTo {T : ℝ} {ρ : ℂ} (hρ : ρ ∈ zerosUpTo T) : ‖ρ‖ ≤ T + 1 := by
  obtain ⟨hz, h0, hT⟩ := hρ
  have h1 := re_pos_of_mem hz
  have h2 := re_lt_one_of_mem hz
  calc ‖ρ‖ ≤ |ρ.re| + |ρ.im| := Complex.norm_le_abs_re_add_abs_im ρ
    _ ≤ 1 + T := by rw [abs_of_pos h1, abs_of_pos h0]; linarith
    _ = T + 1 := by ring

lemma finite_zerosUpTo (T : ℝ) : (zerosUpTo T).Finite :=
  (finite_norm_le (T + 1)).subset fun ρ hρ => ⟨hρ.1, norm_le_of_mem_zerosUpTo hρ⟩

lemma finite_zerosBelow (T : ℝ) : (zerosBelow T).Finite :=
  (finite_zerosUpTo T).subset fun ρ hρ => ⟨hρ.1, hρ.2.1, hρ.2.2.le⟩

end LeanLi
