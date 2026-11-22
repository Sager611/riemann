import Riemann.Deep.LeanNBRate.ExpSumOps

/-!
# Nikolskii inequality for exponential sums, and the abstract form of [N] Lemma B

* `nikolskii`: `|r(u)|² ≤ 2(2c)^{2(M−1)} · avg_{[a,b]} |r|²` on `[a,b]` (Chebyshev set + Turán–Nazarov).
* `lemmaB_abstract`: from Turán, Nikolskii and a Bernstein bound `|r′| ≤ W sup|r|`, with
  `q = icr + r′`, `|r(u)|² ≤ 4T₁²T₂/(ℓc²) ∫_{ℓ/4}^{3ℓ/4} |q|²` on `[0, ℓ]` ([N] Lemma B, its absorption step).
-/

open MeasureTheory Filter Topology Set Complex

noncomputable section

namespace Riemann.Deep.LeanNBRate

/-- Chebyshev: on `[a,b]` the set where `g ≤ 2A'` (`A' > 0`, `A' ≥` the average of `g`) has
measure at least `(b − a)/2`. -/
theorem volume_le_two_avg {g : ℝ → ℝ} {a b A' : ℝ} (hab : a < b) (hg : Continuous g)
    (hg0 : ∀ x, 0 ≤ g x) (hA' : 0 < A') (havg : ∫ x in Icc a b, g x ≤ A' * (b - a)) :
    ENNReal.ofReal ((b - a) / 2) ≤ volume {x ∈ Icc a b | g x ≤ 2 * A'} := by
  set μ := volume.restrict (Icc a b)
  have hint : Integrable g μ :=
    (hg.continuousOn.integrableOn_compact isCompact_Icc)
  have hmarkov := mul_meas_ge_le_integral_of_nonneg (μ := μ) (Eventually.of_forall hg0) hint (2 * A')
  have hμ : μ {x | 2 * A' ≤ g x} ≤ ENNReal.ofReal ((b - a) / 2) := by
    have h1 : (2 * A') * (μ {x | 2 * A' ≤ g x}).toReal ≤ A' * (b - a) := hmarkov.trans havg
    have hfin : μ {x | 2 * A' ≤ g x} ≠ ⊤ := by
      refine ne_top_of_le_ne_top ?_ (measure_mono (subset_univ _))
      rw [Measure.restrict_apply_univ, Real.volume_Icc]; exact ENNReal.ofReal_ne_top
    rw [← ENNReal.ofReal_toReal hfin]
    apply ENNReal.ofReal_le_ofReal
    nlinarith
  have hmeas : MeasurableSet {x | 2 * A' ≤ g x} := measurableSet_le measurable_const hg.measurable
  have htot : volume (Icc a b) = ENNReal.ofReal (b - a) := Real.volume_Icc
  -- volume(Icc ∖ E) ≤ (b−a)/2
  have hsub : Icc a b \ {x ∈ Icc a b | g x ≤ 2 * A'} ⊆ Icc a b ∩ {x | 2 * A' ≤ g x} := by
    intro x hx
    simp only [mem_diff, mem_setOf_eq, not_and, not_le] at hx
    exact ⟨hx.1, (hx.2 hx.1).le⟩
  have heq : volume (Icc a b ∩ {x | 2 * A' ≤ g x}) = μ {x | 2 * A' ≤ g x} := by
    rw [Measure.restrict_apply hmeas, inter_comm]
  have hcomp : volume (Icc a b \ {x ∈ Icc a b | g x ≤ 2 * A'}) ≤ ENNReal.ofReal ((b - a) / 2) :=
    (measure_mono hsub).trans (heq.le.trans hμ)
  have h1 := (le_measure_diff (μ := volume) (s₁ := Icc a b)
    (s₂ := {x ∈ Icc a b | g x ≤ 2 * A'})).trans hcomp
  rw [htot, tsub_le_iff_right] at h1
  have h3 : ENNReal.ofReal (b - a) = ENNReal.ofReal ((b - a) / 2) + ENNReal.ofReal ((b - a) / 2) := by
    rw [← ENNReal.ofReal_add (by linarith) (by linarith)]; ring_nf
  rw [h3] at h1
  exact (ENNReal.add_le_add_iff_left ENNReal.ofReal_ne_top).mp h1


/-- **Nikolskii** for exponential sums with jets: `|r(u)|² ≤ 2(2c)^{2(M−1)} avg_{[a,b]} |r|²`. -/
theorem nikolskii (hT : Nazarov1993_Turan) : ∃ c : ℝ, 0 < c ∧
    ∀ (r : ExpSum) (a b : ℝ), 1 ≤ r.M → a < b → ∀ u ∈ Icc a b,
      ‖r.eval u‖ ^ 2 ≤ 2 * (2 * c) ^ (2 * (r.M - 1)) *
        ((b - a)⁻¹ * ∫ x in Icc a b, ‖r.eval x‖ ^ 2) := by
  obtain ⟨c, hc, hT⟩ := turan_expSum hT
  refine ⟨c, hc, fun r a b hM hab u hu => ?_⟩
  set A := (b - a)⁻¹ * ∫ x in Icc a b, ‖r.eval x‖ ^ 2
  have hba : 0 < b - a := by linarith
  have hgc : Continuous fun x => ‖r.eval x‖ ^ 2 := (r.continuous_eval.norm.pow 2)
  have hA0 : 0 ≤ A := by
    apply mul_nonneg (inv_nonneg.mpr hba.le)
    exact setIntegral_nonneg measurableSet_Icc fun x _ => by positivity
  have key : ∀ ε : ℝ, 0 < ε → ‖r.eval u‖ ^ 2 ≤ 2 * (2 * c) ^ (2 * (r.M - 1)) * (A + ε) := by
    intro ε hε
    set A' := A + ε
    have hA' : 0 < A' := by positivity
    have havg : ∫ x in Icc a b, ‖r.eval x‖ ^ 2 ≤ A' * (b - a) := by
      have : ∫ x in Icc a b, ‖r.eval x‖ ^ 2 = A * (b - a) := by
        simp only [A]; field_simp
      rw [this]; nlinarith
    have hvol := volume_le_two_avg hab hgc (fun x => by positivity) hA' havg
    set E := {x ∈ Icc a b | ‖r.eval x‖ ^ 2 ≤ 2 * A'}
    have hEm : MeasurableSet E :=
      measurableSet_Icc.inter (measurableSet_le hgc.measurable measurable_const)
    have hEpos : 0 < volume E :=
      lt_of_lt_of_le (ENNReal.ofReal_pos.mpr (by linarith)) hvol
    have hEfin : volume E ≠ ⊤ := ne_top_of_le_ne_top (by rw [Real.volume_Icc]; exact ENNReal.ofReal_ne_top)
      (measure_mono fun x hx => hx.1)
    have hEreal : (b - a) / 2 ≤ (volume E).toReal := by
      rw [← ENNReal.ofReal_le_iff_le_toReal hEfin]; exact hvol
    have hK : ∀ x ∈ E, ‖r.eval x‖ ≤ Real.sqrt (2 * A') := by
      intro x hx
      exact Real.le_sqrt_of_sq_le hx.2
    have h := hT r a b E (Real.sqrt (2 * A')) hM hab hEm (fun x hx => hx.1) hEpos hK u hu
    have hratio : c * (b - a) / (volume E).toReal ≤ 2 * c := by
      rw [div_le_iff₀ (by linarith)]; nlinarith
    have h2 : ‖r.eval u‖ ≤ (2 * c) ^ (r.M - 1) * Real.sqrt (2 * A') :=
      h.trans (mul_le_mul_of_nonneg_right (pow_le_pow_left₀ (by positivity) hratio _)
        (Real.sqrt_nonneg _))
    calc ‖r.eval u‖ ^ 2 ≤ ((2 * c) ^ (r.M - 1) * Real.sqrt (2 * A')) ^ 2 := by
          gcongr
      _ = 2 * (2 * c) ^ (2 * (r.M - 1)) * A' := by
          rw [mul_pow, Real.sq_sqrt (by positivity), ← pow_mul, mul_comm (r.M - 1) 2]; ring
  refine le_of_forall_pos_lt_add fun δ hδ => ?_
  set C := 2 * (2 * c) ^ (2 * (r.M - 1))
  have hC : 0 < C := by positivity
  have := key (δ / (2 * C)) (by positivity)
  calc ‖r.eval u‖ ^ 2 ≤ C * (A + δ / (2 * C)) := this
    _ = C * A + δ / 2 := by field_simp
    _ < C * A + δ := by linarith

/-- **[N] Lemma B, abstract form** (the absorption step). On `I = [ℓ/4, 3ℓ/4]`: Turán `T₁`,
Nikolskii `T₂`, Bernstein `W`; `q = icr + r′`. -/
theorem lemmaB_abstract {r rd q : ℝ → ℂ} {ℓ c T₁ T₂ W : ℝ} (hℓ : 0 < ℓ) (hT₂ : 0 < T₂)
    (hT₁ : 0 ≤ T₁) (hW : 0 ≤ W)
    (hr : Continuous r) (hrd : Continuous rd) (hqc : Continuous q)
    (hTur : ∀ X, (∀ u ∈ Icc (ℓ / 4) (3 * ℓ / 4), ‖r u‖ ≤ X) → ∀ u ∈ Icc 0 ℓ, ‖r u‖ ≤ T₁ * X)
    (hNik : ∀ u ∈ Icc (ℓ / 4) (3 * ℓ / 4), ‖r u‖ ^ 2 ≤ T₂ / ℓ * ∫ x in Icc (ℓ / 4) (3 * ℓ / 4), ‖r x‖ ^ 2)
    (hBern : ∀ X, (∀ u ∈ Icc 0 ℓ, ‖r u‖ ≤ X) → ∀ u ∈ Icc 0 ℓ, ‖rd u‖ ≤ W * X)
    (hq : ∀ u ∈ Icc 0 ℓ, q u = I * c * r u + rd u)
    (hsmall : T₂ * W ^ 2 * T₁ ^ 2 ≤ c ^ 2 / 2) (hc : c ≠ 0) :
    ∀ u ∈ Icc 0 ℓ, ‖r u‖ ^ 2 ≤ 4 * T₁ ^ 2 * T₂ / (ℓ * c ^ 2) *
      ∫ x in Icc (ℓ / 4) (3 * ℓ / 4), ‖q x‖ ^ 2 := by
  set Iv := Icc (ℓ / 4) (3 * ℓ / 4)
  set Jr := ∫ x in Iv, ‖r x‖ ^ 2
  set Jq := ∫ x in Iv, ‖q x‖ ^ 2
  set X₀ := Real.sqrt (T₂ / ℓ * Jr)
  have hIsub : Iv ⊆ Icc 0 ℓ := Icc_subset_Icc (by linarith) (by linarith)
  have hJr0 : 0 ≤ Jr := setIntegral_nonneg measurableSet_Icc fun x _ => by positivity
  have hX₀ : ∀ u ∈ Iv, ‖r u‖ ≤ X₀ := fun u hu => Real.le_sqrt_of_sq_le (hNik u hu)
  have hrb : ∀ u ∈ Icc 0 ℓ, ‖r u‖ ≤ T₁ * X₀ := hTur X₀ hX₀
  have hrdb : ∀ u ∈ Icc 0 ℓ, ‖rd u‖ ≤ W * (T₁ * X₀) := hBern _ hrb
  have hint : ∀ f : ℝ → ℂ, Continuous f → IntegrableOn (fun x => ‖f x‖ ^ 2) Iv :=
    fun f hf => (hf.norm.pow 2).continuousOn.integrableOn_compact isCompact_Icc
  -- c² ∫|r|² ≤ 2∫|q|² + 2∫|rd|²
  have hpt : ∀ x ∈ Iv, c ^ 2 * ‖r x‖ ^ 2 ≤ 2 * ‖q x‖ ^ 2 + 2 * ‖rd x‖ ^ 2 := by
    intro x hx
    have hqx := hq x (hIsub hx)
    have : (I * c * r x : ℂ) = q x - rd x := by rw [hqx]; ring
    have h1 : ‖(I * c * r x : ℂ)‖ = |c| * ‖r x‖ := by
      rw [norm_mul, norm_mul, Complex.norm_I, one_mul, Complex.norm_real, Real.norm_eq_abs]
    have h2 : ‖q x - rd x‖ ≤ ‖q x‖ + ‖rd x‖ := norm_sub_le _ _
    rw [this] at h1
    have h3 : |c| * ‖r x‖ ≤ ‖q x‖ + ‖rd x‖ := h1 ▸ h2
    have h4 : (|c| * ‖r x‖) ^ 2 ≤ (‖q x‖ + ‖rd x‖) ^ 2 := by gcongr
    rw [mul_pow, sq_abs] at h4
    nlinarith [sq_nonneg (‖q x‖ - ‖rd x‖)]
  have hrd2 : ∫ x in Iv, ‖rd x‖ ^ 2 ≤ ℓ / 2 * (W * (T₁ * X₀)) ^ 2 := by
    have : ∫ x in Iv, ‖rd x‖ ^ 2 ≤ ∫ x in Iv, (W * (T₁ * X₀)) ^ 2 := by
      refine setIntegral_mono_on (hint rd hrd) (integrableOn_const (by
        rw [Real.volume_Icc]; exact ENNReal.ofReal_ne_top) ) measurableSet_Icc fun x hx => ?_
      gcongr
      exact hrdb x (hIsub hx)
    refine this.trans (le_of_eq ?_)
    rw [setIntegral_const, Real.volume_real_Icc_of_le (by linarith), smul_eq_mul]; ring
  have hmain : c ^ 2 * Jr ≤ 2 * Jq + 2 * (ℓ / 2 * (W * (T₁ * X₀)) ^ 2) := by
    have h1 : ∫ x in Iv, c ^ 2 * ‖r x‖ ^ 2 ≤ ∫ x in Iv, (2 * ‖q x‖ ^ 2 + 2 * ‖rd x‖ ^ 2) :=
      setIntegral_mono_on ((hint r hr).const_mul _)
        (((hint q hqc).const_mul _).add ((hint rd hrd).const_mul _)) measurableSet_Icc hpt
    rw [integral_const_mul, integral_add ((hint q hqc).const_mul _) ((hint rd hrd).const_mul _),
      integral_const_mul, integral_const_mul] at h1
    linarith
  -- X₀² = (T₂/ℓ) Jr
  have hX₀sq : X₀ ^ 2 = T₂ / ℓ * Jr := Real.sq_sqrt (by positivity)
  have hJr : Jr = ℓ / T₂ * X₀ ^ 2 := by rw [hX₀sq]; field_simp
  have hc2 : 0 < c ^ 2 := by positivity
  have hX₀b : X₀ ^ 2 ≤ 4 * T₂ / (ℓ * c ^ 2) * Jq := by
    rw [hJr] at hmain
    have h5 : ℓ * X₀ ^ 2 * (c ^ 2 / T₂ - W ^ 2 * T₁ ^ 2) ≤ 2 * Jq := by
      have e1 : c ^ 2 * (ℓ / T₂ * X₀ ^ 2) = ℓ * X₀ ^ 2 * (c ^ 2 / T₂) := by ring
      have e2 : 2 * (ℓ / 2 * (W * (T₁ * X₀)) ^ 2) = ℓ * X₀ ^ 2 * (W ^ 2 * T₁ ^ 2) := by ring
      rw [e1, e2] at hmain
      rw [mul_sub]; linarith
    have h6 : c ^ 2 / (2 * T₂) ≤ c ^ 2 / T₂ - W ^ 2 * T₁ ^ 2 := by
      rw [le_sub_iff_add_le, div_add' _ _ _ (by positivity), div_le_div_iff₀ (by positivity) hT₂]
      nlinarith
    have h7 : ℓ * X₀ ^ 2 * (c ^ 2 / (2 * T₂)) ≤ 2 * Jq :=
      le_trans (by gcongr) h5
    rw [div_mul_eq_mul_div, le_div_iff₀ (by positivity)]
    have : ℓ * X₀ ^ 2 * c ^ 2 ≤ 4 * T₂ * Jq := by
      have := h7; field_simp at this; nlinarith
    linarith
  intro u hu
  calc ‖r u‖ ^ 2 ≤ (T₁ * X₀) ^ 2 := by gcongr; exact hrb u hu
    _ = T₁ ^ 2 * X₀ ^ 2 := by ring
    _ ≤ T₁ ^ 2 * (4 * T₂ / (ℓ * c ^ 2) * Jq) := by gcongr
    _ = 4 * T₁ ^ 2 * T₂ / (ℓ * c ^ 2) * Jq := by ring

end Riemann.Deep.LeanNBRate
