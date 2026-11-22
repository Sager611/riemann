import Riemann.Deep.LeanNBRate.LitInputs
import Mathlib.NumberTheory.LSeries.ZetaZeros
import Mathlib.NumberTheory.Harmonic.ZetaAsymp
import Mathlib.Analysis.Calculus.Deriv.Star
import Mathlib.Analysis.Complex.CauchyIntegral

/-!
# Zeros: finiteness in boxes, conjugation symmetry (with multiplicity)
-/

open MeasureTheory Filter Topology Set Complex ComplexConjugate

noncomputable section

namespace Riemann.Deep.LeanNBRate

theorem NTZ_subset_zeros : NTZ ⊆ riemannZetaZeros := fun ρ h => h.1

theorem zerosIn_finite (a b : ℝ) : (zerosIn a b).Finite := by
  have hK : IsCompact (Metric.closedBall (0 : ℂ) (1 + |a| + |b|)) := isCompact_closedBall _ _
  refine (hK.inter_riemannZetaZeros_finite).subset fun ρ hρ => ⟨?_, hρ.1.1⟩
  obtain ⟨⟨_, h0, h1⟩, ha, hb⟩ := hρ
  rw [Metric.mem_closedBall, dist_zero_right]
  have h2 : |ρ.im| ≤ |a| + |b| := by
    have := abs_nonneg a
    have := abs_nonneg b
    rw [abs_le]; constructor
    · have := neg_abs_le a; linarith
    · have := le_abs_self b; linarith
  calc ‖ρ‖ ≤ |ρ.re| + |ρ.im| := Complex.norm_le_abs_re_add_abs_im ρ
    _ ≤ 1 + (|a| + |b|) := by
        gcongr
        rw [abs_of_pos h0]; linarith
    _ = 1 + |a| + |b| := by ring

theorem NTZ_ne_one {ρ : ℂ} (h : ρ ∈ NTZ) : ρ ≠ 1 := by
  intro h1; rw [h1] at h; exact absurd h.2.2 (by simp)

theorem conj_mem_NTZ {ρ : ℂ} (h : ρ ∈ NTZ) : conj ρ ∈ NTZ := by
  obtain ⟨hz, h0, h1⟩ := h
  refine ⟨?_, by simpa using h0, by simpa using h1⟩
  rw [riemannZeta_conj, hz, map_zero]

theorem AnalyticAt.conj_conj' {f : ℂ → ℂ} {z : ℂ} (hf : AnalyticAt ℂ f z) :
    AnalyticAt ℂ (fun w => conj (f (conj w))) (conj z) := by
  set V := {x | AnalyticAt ℂ f x}
  have hV : V ∈ 𝓝 z := (isOpen_analyticAt ℂ f).mem_nhds hf
  have hs : (fun w => conj w) ⁻¹' V ∈ 𝓝 (conj z) := by
    have := Complex.continuous_conj.continuousAt (x := conj z)
    rw [ContinuousAt, Complex.conj_conj] at this
    exact this hV
  refine DifferentiableOn.analyticAt (s := (fun w => conj w) ⁻¹' V) (fun w hw => ?_) hs
  have hd : DifferentiableAt ℂ f (conj w) := (show AnalyticAt ℂ f (conj w) from hw).differentiableAt
  have := hd.conj_conj
  rw [Complex.conj_conj] at this
  exact this.differentiableWithinAt

theorem analyticOrderAt_conj_conj {f : ℂ → ℂ} {z : ℂ} (hf : AnalyticAt ℂ f z) :
    analyticOrderAt (fun w => conj (f (conj w))) (conj z) = analyticOrderAt f z := by
  have hcont : Tendsto (fun w : ℂ => conj w) (𝓝 (conj z)) (𝓝 z) := by
    have := Complex.continuous_conj.tendsto (conj z)
    rwa [Complex.conj_conj] at this
  rcases eq_or_ne (analyticOrderAt f z) ⊤ with htop | hfin
  · rw [htop, analyticOrderAt_eq_top]
    rw [analyticOrderAt_eq_top] at htop
    filter_upwards [hcont.eventually htop] with w hw
    rw [hw, map_zero]
  · obtain ⟨g, hg, hg0, hfg⟩ := (hf.analyticOrderAt_ne_top).mp hfin
    rw [← Nat.cast_analyticOrderNatAt hfin]
    refine ((AnalyticAt.conj_conj' hf).analyticOrderAt_eq_natCast).mpr
      ⟨fun w => conj (g (conj w)), AnalyticAt.conj_conj' hg,
      by simpa using hg0, ?_⟩
    filter_upwards [hcont.eventually hfg] with w hw
    rw [hw, smul_eq_mul, smul_eq_mul, map_mul, map_pow, map_sub, Complex.conj_conj]

theorem mult_conj (ρ : ℂ) (hρ : ρ ≠ 1) : mult (conj ρ) = mult ρ := by
  have hζ : (fun w => conj (riemannZeta (conj w))) = riemannZeta := by
    funext w; rw [riemannZeta_conj, Complex.conj_conj]
  have han : AnalyticAt ℂ riemannZeta ρ := analyticOn_riemannZeta ρ hρ
  unfold mult analyticOrderNatAt
  rw [← analyticOrderAt_conj_conj han, hζ]

/-- Under RH a nontrivial zero is ½ + iγ. -/
theorem NTZ_re_half (hRH : RiemannHypothesis) {ρ : ℂ} (h : ρ ∈ NTZ) : ρ.re = 1 / 2 := by
  obtain ⟨hz, h0, h1⟩ := h
  refine hRH ρ hz ?_ (NTZ_ne_one ⟨hz, h0, h1⟩)
  rintro ⟨n, rfl⟩
  simp at h0
  have : (0 : ℝ) ≤ n := n.cast_nonneg
  linarith

theorem NTZ_eq (hRH : RiemannHypothesis) {ρ : ℂ} (h : ρ ∈ NTZ) : ρ = (1 / 2 : ℂ) + ρ.im * I := by
  apply Complex.ext <;> simp [NTZ_re_half hRH h]

theorem Ncount_eq (a b : ℝ) : (Ncount a b : ℕ) = ∑ ρ ∈ (zerosIn_finite a b).toFinset, mult ρ := by
  unfold Ncount
  exact finsum_mem_eq_finite_toFinset_sum _ (zerosIn_finite a b)

/-- A finite family of distinct zeros with ordinates in (a, b] counts at most Ncount a b. -/
theorem sum_le_Ncount {a b : ℝ} (u : Finset ℂ) (hu : ∀ ρ ∈ u, ρ ∈ zerosIn a b) :
    ∑ ρ ∈ u, mult ρ ≤ Ncount a b := by
  rw [Ncount_eq]
  exact Finset.sum_le_sum_of_subset (fun ρ hρ => (Set.Finite.mem_toFinset _).mpr (hu ρ hρ))

/-- The shell counts c_k = 2 (N(k) - N(k-1) + N(k+1) - N(k)). -/
def shellC (k : ℕ) : ℕ := 2 * (Ncount ((k : ℝ) - 1) k + Ncount k ((k : ℝ) + 1))

/-- Shell bound for a finite set of nontrivial zeros. -/
theorem shell_sum_le (u : Finset ℂ) (hu : ∀ ρ ∈ u, ρ ∈ NTZ) (k : ℕ) :
    ∑ ρ ∈ u.filter (fun ρ => ⌊|ρ.im|⌋₊ = k), mult ρ ≤ shellC k := by
  classical
  set v := u.filter (fun ρ => ⌊|ρ.im|⌋₊ = k)
  have hsplit : ∑ ρ ∈ v, mult ρ = ∑ ρ ∈ v.filter (fun ρ => 0 ≤ ρ.im), mult ρ +
      ∑ ρ ∈ v.filter (fun ρ => ¬ 0 ≤ ρ.im), mult ρ := (Finset.sum_filter_add_sum_filter_not _ _ _).symm
  have hfl : ∀ ρ ∈ v, (k : ℝ) ≤ |ρ.im| ∧ |ρ.im| < k + 1 := by
    intro ρ hρ
    have hk := (Finset.mem_filter.mp hρ).2
    refine ⟨?_, ?_⟩
    · rw [← hk]; exact Nat.floor_le (abs_nonneg _)
    · rw [← hk]; exact Nat.lt_floor_add_one _
  -- positive part
  have hpos : ∑ ρ ∈ v.filter (fun ρ => 0 ≤ ρ.im), mult ρ ≤
      Ncount ((k : ℝ) - 1) k + Ncount k ((k : ℝ) + 1) := by
    set w := v.filter (fun ρ => 0 ≤ ρ.im)
    rw [← Finset.sum_filter_add_sum_filter_not w (fun ρ => ρ.im ≤ k)]
    gcongr
    · refine sum_le_Ncount _ fun ρ hρ => ?_
      obtain ⟨hρw, hle⟩ := Finset.mem_filter.mp hρ
      obtain ⟨hρv, h0⟩ := Finset.mem_filter.mp hρw
      have := hfl ρ hρv
      rw [abs_of_nonneg h0] at this
      exact ⟨hu ρ (Finset.mem_filter.mp hρv).1, by linarith [this.1], hle⟩
    · refine sum_le_Ncount _ fun ρ hρ => ?_
      obtain ⟨hρw, hle⟩ := Finset.mem_filter.mp hρ
      obtain ⟨hρv, h0⟩ := Finset.mem_filter.mp hρw
      have := hfl ρ hρv
      rw [abs_of_nonneg h0] at this
      exact ⟨hu ρ (Finset.mem_filter.mp hρv).1, lt_of_not_ge hle, this.2.le⟩
  -- negative part, via conjugation
  have hneg : ∑ ρ ∈ v.filter (fun ρ => ¬ 0 ≤ ρ.im), mult ρ ≤
      Ncount ((k : ℝ) - 1) k + Ncount k ((k : ℝ) + 1) := by
    set w := v.filter (fun ρ => ¬ 0 ≤ ρ.im)
    have hw : ∑ ρ ∈ w, mult ρ = ∑ ρ ∈ w.image (fun ρ => conj ρ), mult ρ := by
      rw [Finset.sum_image (fun a _ b _ h => (starRingEnd ℂ).injective h)]
      refine Finset.sum_congr rfl fun ρ hρ => ?_
      have hρu : ρ ∈ u := (Finset.mem_filter.mp (Finset.mem_filter.mp hρ).1).1
      exact (mult_conj ρ (NTZ_ne_one (hu ρ hρu))).symm
    rw [hw]
    set w' := w.image (fun ρ => conj ρ)
    have hw' : ∀ ρ ∈ w', ρ ∈ NTZ ∧ (k : ℝ) ≤ ρ.im ∧ ρ.im < k + 1 := by
      intro ρ hρ
      obtain ⟨σ, hσ, rfl⟩ := Finset.mem_image.mp hρ
      obtain ⟨hσv, h0⟩ := Finset.mem_filter.mp hσ
      have := hfl σ hσv
      rw [abs_of_neg (lt_of_not_ge h0)] at this
      refine ⟨conj_mem_NTZ (hu σ (Finset.mem_filter.mp hσv).1), ?_, ?_⟩ <;> simp <;> linarith [this.1, this.2]
    rw [← Finset.sum_filter_add_sum_filter_not w' (fun ρ => ρ.im ≤ k)]
    gcongr
    · refine sum_le_Ncount _ fun ρ hρ => ?_
      obtain ⟨hρw, hle⟩ := Finset.mem_filter.mp hρ
      obtain ⟨h1, h2, h3⟩ := hw' ρ hρw
      exact ⟨h1, by linarith, hle⟩
    · refine sum_le_Ncount _ fun ρ hρ => ?_
      obtain ⟨hρw, hle⟩ := Finset.mem_filter.mp hρ
      obtain ⟨h1, h2, h3⟩ := hw' ρ hρw
      exact ⟨h1, lt_of_not_ge hle, h3.le⟩
  rw [hsplit]
  unfold shellC
  omega

/-- Summability over the zeros from the shell counts. -/
theorem summable_NTZ_of_shells (F : ℝ → ℝ) (hF0 : ∀ t, 0 ≤ F t) (hFa : AntitoneOn F (Ici 0))
    (hS : Summable fun k : ℕ => F k * shellC k) :
    Summable fun ρ : NTZ => (mult (ρ : ℂ) : ℝ) * F |(ρ : ℂ).im| := by
  classical
  refine summable_of_sum_le (c := ∑' k : ℕ, F k * shellC k)
    (fun ρ => mul_nonneg (Nat.cast_nonneg _) (hF0 _)) ?_
  intro u
  set U : Finset ℂ := u.map (Function.Embedding.subtype _)
  have hU : ∀ ρ ∈ U, ρ ∈ NTZ := by
    intro ρ hρ
    obtain ⟨σ, _, rfl⟩ := Finset.mem_map.mp hρ
    exact σ.2
  have e1 : ∑ ρ ∈ u, (mult (ρ : ℂ) : ℝ) * F |(ρ : ℂ).im| = ∑ ρ ∈ U, (mult ρ : ℝ) * F |ρ.im| := by
    simp [U]
  rw [e1]
  set K := (U.image fun ρ => ⌊|ρ.im|⌋₊).sup id + 1
  have e2 : ∑ ρ ∈ U, (mult ρ : ℝ) * F |ρ.im| =
      ∑ k ∈ Finset.range K, ∑ ρ ∈ U.filter (fun ρ => ⌊|ρ.im|⌋₊ = k), (mult ρ : ℝ) * F |ρ.im| := by
    rw [Finset.sum_fiberwise_of_maps_to]
    intro ρ hρ
    rw [Finset.mem_range]
    have : ⌊|ρ.im|⌋₊ ≤ (U.image fun ρ => ⌊|ρ.im|⌋₊).sup id :=
      Finset.le_sup (f := id) (Finset.mem_image_of_mem _ hρ)
    omega
  rw [e2]
  have hterm : ∀ k ∈ Finset.range K, ∑ ρ ∈ U.filter (fun ρ => ⌊|ρ.im|⌋₊ = k), (mult ρ : ℝ) * F |ρ.im| ≤
      F k * shellC k := by
    intro k _
    have h1 : ∀ ρ ∈ U.filter (fun ρ => ⌊|ρ.im|⌋₊ = k), (mult ρ : ℝ) * F |ρ.im| ≤ (mult ρ : ℝ) * F k := by
      intro ρ hρ
      have hk := (Finset.mem_filter.mp hρ).2
      have hle : (k : ℝ) ≤ |ρ.im| := by rw [← hk]; exact Nat.floor_le (abs_nonneg _)
      have := hFa (Set.mem_Ici.mpr (Nat.cast_nonneg k)) (Set.mem_Ici.mpr (abs_nonneg ρ.im)) hle
      exact mul_le_mul_of_nonneg_left this (Nat.cast_nonneg _)
    refine (Finset.sum_le_sum h1).trans ?_
    rw [← Finset.sum_mul, mul_comm]
    gcongr
    · exact hF0 k
    · exact_mod_cast shell_sum_le U hU k
  refine (Finset.sum_le_sum hterm).trans ?_
  exact hS.sum_le_tsum _ fun k _ => mul_nonneg (hF0 k) (Nat.cast_nonneg _)

end Riemann.Deep.LeanNBRate


