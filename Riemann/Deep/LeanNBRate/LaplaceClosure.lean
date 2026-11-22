import Riemann.Deep.LeanNBRate.Closure

/-!
# Closure: Laplace transforms on [0, L] (Riemann sums, dominated convergence)
-/

open MeasureTheory Filter Topology Set Complex

noncomputable section

namespace Riemann.Deep.LeanNBRate

theorem riemann_sum_tendsto {φ : ℝ → ℂ} (hφ : Continuous φ) {L : ℝ} (hL : 0 < L) :
    Tendsto (fun n : ℕ => ∑ k ∈ Finset.range n, ((L / n : ℝ) : ℂ) * φ (k * L / n)) atTop
      (𝓝 (∫ u in (0 : ℝ)..L, φ u)) := by
  rw [Metric.tendsto_atTop]
  intro ε hε
  have huc := (isCompact_Icc (a := (0 : ℝ)) (b := L)).uniformContinuousOn_of_continuous hφ.continuousOn
  rw [Metric.uniformContinuousOn_iff] at huc
  obtain ⟨δ, hδ, hδφ⟩ := huc (ε / (2 * L)) (by positivity)
  obtain ⟨N, hN⟩ := exists_nat_gt (L / δ)
  refine ⟨N + 1, fun n hn => ?_⟩
  have hn1 : 1 ≤ n := by omega
  have hn0 : (0 : ℝ) < n := by exact_mod_cast (show 0 < n by omega)
  have hNn : (N : ℝ) < n := by exact_mod_cast (show N < n by omega)
  have hLn : L / n < δ := by
    rw [div_lt_iff₀ hn0]
    have : L / δ < n := lt_trans hN hNn
    rw [div_lt_iff₀ hδ] at this
    linarith
  set a : ℕ → ℝ := fun k => k * L / n with ha
  have ha0 : a 0 = 0 := by simp [a]
  have han : a n = L := by simp only [a]; field_simp
  have hint : ∀ k < n, IntervalIntegrable φ volume (a k) (a (k + 1)) :=
    fun k _ => hφ.intervalIntegrable _ _
  have hsplit := intervalIntegral.sum_integral_adjacent_intervals hint
  rw [ha0, han] at hsplit
  rw [dist_eq_norm, ← hsplit, ← Finset.sum_sub_distrib]
  have hLn0 : 0 < L / n := div_pos hL hn0
  have hterm : ∀ k ∈ Finset.range n, ‖((L / n : ℝ) : ℂ) * φ (k * L / n) - ∫ u in a k..a (k + 1), φ u‖ ≤
      L / n * (ε / (2 * L)) := by
    intro k hk
    have hk' : k < n := Finset.mem_range.mp hk
    have hlen : a (k + 1) - a k = L / n := by simp only [a]; push_cast; field_simp; ring
    have hle : a k ≤ a (k + 1) := by linarith
    have hak0 : 0 ≤ a k := by simp only [a]; positivity
    have hakL : a (k + 1) ≤ L := by
      simp only [a]
      rw [div_le_iff₀ hn0]
      have : ((k + 1 : ℕ) : ℝ) ≤ n := by exact_mod_cast hk'
      nlinarith
    have e : ((L / n : ℝ) : ℂ) * φ (k * L / n) = ∫ _u in a k..a (k + 1), φ (a k) := by
      rw [intervalIntegral.integral_const, hlen, Complex.real_smul]
    rw [e, ← intervalIntegral.integral_sub intervalIntegrable_const (hint k hk')]
    refine (intervalIntegral.norm_integral_le_of_norm_le_const (C := ε / (2 * L)) ?_).trans ?_
    · intro u hu
      rw [Set.uIoc_of_le hle] at hu
      have hu' : u ∈ Icc 0 L := ⟨by linarith [hu.1], by linarith [hu.2]⟩
      have hak : a k ∈ Icc 0 L := ⟨hak0, by linarith⟩
      have hd : dist (a k) u < δ := by
        rw [Real.dist_eq, abs_of_nonpos (by linarith [hu.1])]
        linarith [hu.2]
      have := hδφ (a k) hak u hu' hd
      rw [dist_eq_norm] at this
      exact this.le
    · rw [abs_of_nonneg (by linarith : 0 ≤ a (k + 1) - a k), hlen, mul_comm]
  calc ‖∑ k ∈ Finset.range n, (((L / n : ℝ) : ℂ) * φ (k * L / n) - ∫ u in a k..a (k + 1), φ u)‖
      ≤ ∑ k ∈ Finset.range n, ‖((L / n : ℝ) : ℂ) * φ (k * L / n) - ∫ u in a k..a (k + 1), φ u‖ :=
        norm_sum_le _ _
    _ ≤ ∑ k ∈ Finset.range n, L / n * (ε / (2 * L)) := Finset.sum_le_sum hterm
    _ = ε / 2 := by
        rw [Finset.sum_const, Finset.card_range, nsmul_eq_mul]
        field_simp
    _ < ε := by linarith

theorem memLp_inv_line : MemLp (fun t : ℝ => 1 / ((1 / 2 : ℂ) + t * I)) 2 := by
  have hcont : Continuous (fun t : ℝ => 1 / ((1 / 2 : ℂ) + t * I)) :=
    continuous_const.div (continuous_const.add (Complex.continuous_ofReal.mul continuous_const))
      fun t h => by
        have := congrArg Complex.re h
        simp at this
  have hm : AEStronglyMeasurable (fun t : ℝ => 1 / ((1 / 2 : ℂ) + t * I)) volume :=
    hcont.aestronglyMeasurable
  rw [memLp_two_iff_integrable_sq_norm hm]
  have hi : Integrable (fun t : ℝ => 4 * (1 + (2 * t) ^ 2)⁻¹) :=
    ((integrable_inv_one_add_sq).comp_mul_left' (by norm_num : (2 : ℝ) ≠ 0)).const_mul 4
  refine hi.congr (Filter.Eventually.of_forall fun t => ?_)
  dsimp only
  have hne : (1 / 2 : ℂ) + t * I ≠ 0 := fun h => by
    have := congrArg Complex.re h; simp at this
  rw [norm_div, norm_one, div_pow, one_pow]
  have hn : ‖(1 / 2 : ℂ) + t * I‖ ^ 2 = 1 / 4 + t ^ 2 := by
    rw [← Complex.normSq_eq_norm_sq, Complex.normSq_apply]; simp; ring
  rw [hn]
  field_simp
  ring

theorem memLp_zeta_line (hP : Rudin1987_Thm9_13) (h215 : Titchmarsh_2_1_5) :
    MemLp (fun t : ℝ => riemannZeta (1 / 2 + t * I) / ((1 / 2 : ℂ) + t * I)) 2 := by
  have h := memLp_line_nbF hP h215 (n := 1) (fun _ => -1) (fun _ => 1) (fun _ => ⟨one_pos, le_rfl⟩)
  have e : (fun t : ℝ => (1 + riemannZeta (1 / 2 + t * I) * ∑ k : Fin 1, (fun _ => (-1 : ℂ)) k *
      ((fun _ => (1 : ℝ)) k : ℂ) ^ (1 / 2 + t * I)) / (1 / 2 + t * I)) =
      fun t : ℝ => 1 / ((1 / 2 : ℂ) + t * I) - riemannZeta (1 / 2 + t * I) / ((1 / 2 : ℂ) + t * I) := by
    funext t
    simp only [Finset.univ_unique, Finset.sum_singleton, Complex.ofReal_one, Complex.one_cpow]
    ring
  rw [e] at h
  have h2 := memLp_inv_line.sub h
  convert h2 using 1
  funext t
  simp only [Pi.sub_apply]
  ring


theorem continuous_zeta_line : Continuous fun t : ℝ => riemannZeta (1 / 2 + t * I) := by
  rw [continuous_iff_continuousAt]
  intro t
  have hne : (1 / 2 : ℂ) + t * I ≠ 1 := fun h => by
    have := congrArg Complex.re h; norm_num at this
  have hc : Continuous (fun t : ℝ => (1 / 2 : ℂ) + t * I) :=
    continuous_const.add (Complex.continuous_ofReal.mul continuous_const)
  exact ContinuousAt.comp (f := fun t : ℝ => (1 / 2 : ℂ) + t * I)
    (differentiableAt_riemannZeta hne).continuousAt hc.continuousAt

theorem norm_cexp_line {u t : ℝ} (hu : 0 ≤ u) : ‖Complex.exp (-(u * ((1 / 2 : ℂ) + t * I)))‖ ≤ 1 := by
  rw [Complex.norm_exp]
  apply Real.exp_le_one_iff.mpr
  simp
  nlinarith

/-- **Closure for Laplace transforms on [0, L]**: a constant plus the Laplace transform of a continuous
function, vanishing at s = 1, bounds Nyman's distance at λ = e^{-L}. -/
theorem nymanDist_sq_le_laplace (hP : Rudin1987_Thm9_13) (h215 : Titchmarsh_2_1_5) {L : ℝ} (hL : 0 < L)
    {g : ℝ → ℂ} (hg : Continuous g) (c₀ : ℂ)
    (h1 : c₀ + ∫ u in (0 : ℝ)..L, g u * Complex.exp (-(u : ℂ)) = 0) :
    nymanDist (Real.exp (-L)) ^ 2 ≤ (2 * Real.pi)⁻¹ * ∫ t : ℝ,
      ‖(1 + riemannZeta (1 / 2 + t * I) *
        (c₀ + ∫ u in (0 : ℝ)..L, g u * Complex.exp (-(u * (1 / 2 + t * I))))) / (1 / 2 + t * I)‖ ^ 2 := by
  classical
  set A : ℂ → ℂ := fun s => c₀ + ∫ u in (0 : ℝ)..L, g u * Complex.exp (-(u * s)) with hA
  set S : ℕ → ℂ → ℂ := fun n s => ∑ k ∈ Finset.range n, ((L / n : ℝ) : ℂ) * (g (k * L / n) *
    Complex.exp (-(((k * L / n : ℝ) : ℂ) * s))) with hS
  -- the finite combinations
  have hfin : ∀ n : ℕ, 1 ≤ n → nymanDist (Real.exp (-L)) ^ 2 ≤ (2 * Real.pi)⁻¹ * ∫ t : ℝ,
      ‖(1 + riemannZeta (1 / 2 + t * I) * (S n (1 / 2 + t * I) - S n 1)) / (1 / 2 + t * I)‖ ^ 2 := by
    intro n hn
    set θ : Fin n → ℝ := fun k => Real.exp (-((k : ℕ) * L / n))
    set c : Fin n → ℂ := fun k => ((L / n : ℝ) : ℂ) * g ((k : ℕ) * L / n) - (if (k : ℕ) = 0 then S n 1 else 0)
    have hn0 : (0 : ℝ) < n := by exact_mod_cast (show 0 < n by omega)
    have hθ : ∀ k, θ k ∈ Icc (Real.exp (-L)) 1 := by
      intro k
      have hk : ((k : ℕ) : ℝ) < n := by exact_mod_cast k.2
      have hk0 : (0 : ℝ) ≤ (k : ℕ) := Nat.cast_nonneg _
      constructor
      · apply Real.exp_le_exp.mpr
        have : ((k : ℕ) : ℝ) * L / n ≤ L := by
          rw [div_le_iff₀ hn0]; nlinarith
        linarith
      · apply Real.exp_le_one_iff.mpr
        have : 0 ≤ ((k : ℕ) : ℝ) * L / n := by positivity
        linarith
    have hsum : ∀ s : ℂ, ∑ k, c k * (θ k : ℂ) ^ s = S n s - S n 1 := by
      intro s
      simp only [c, θ, sub_mul, Finset.sum_sub_distrib]
      congr 1
      · rw [Fin.sum_univ_eq_sum_range (fun k : ℕ => ((L / n : ℝ) : ℂ) * g (k * L / n) *
          ((Real.exp (-(k * L / n)) : ℝ) : ℂ) ^ s) n]
        refine Finset.sum_congr rfl fun k _ => ?_
        rw [ofReal_exp_cpow]
        push_cast
        ring_nf
      · rw [Finset.sum_eq_single ⟨0, by omega⟩]
        · simp
        · intro b _ hb
          have : (b : ℕ) ≠ 0 := fun h => hb (Fin.ext h)
          simp [this]
        · intro h; exact absurd (Finset.mem_univ _) h
    have hc : (∑ k, c k * (θ k : ℂ)).re = 0 := by
      have := hsum 1
      simp only [Complex.cpow_one] at this
      rw [this, sub_self, Complex.zero_re]
    have h := nymanDist_sq_le_line hP h215 (Real.exp_pos _) c θ hθ hc
    simp_rw [hsum] at h
    exact h
  -- convergence of the finite combinations
  have hconv : ∀ s : ℂ, Tendsto (fun n => S n s) atTop (𝓝 (∫ u in (0 : ℝ)..L, g u * Complex.exp (-(u * s)))) := by
    intro s
    have hφ : Continuous fun u : ℝ => g u * Complex.exp (-((u : ℂ) * s)) :=
      hg.mul (Complex.continuous_exp.comp ((Complex.continuous_ofReal.mul continuous_const).neg))
    have := riemann_sum_tendsto hφ hL
    refine this.congr (fun n => ?_)
    simp only [S]
  have hA1 : ∫ u in (0 : ℝ)..L, g u * Complex.exp (-(u * (1 : ℂ))) = -c₀ := by
    simp only [mul_one]; linear_combination h1
  -- bounds
  obtain ⟨M, hM⟩ := (isCompact_Icc (a := (0 : ℝ)) (b := L)).exists_bound_of_continuousOn hg.continuousOn
  set M' := max M 0
  have hM' : ∀ u ∈ Icc (0 : ℝ) L, ‖g u‖ ≤ M' := fun u hu => (hM u hu).trans (le_max_left _ _)
  have hM'0 : 0 ≤ M' := le_max_right _ _
  have hSb : ∀ n : ℕ, ∀ s : ℂ, (∀ u : ℝ, 0 ≤ u → ‖Complex.exp (-((u : ℂ) * s))‖ ≤ 1) →
      ‖S n s‖ ≤ L * M' := by
    intro n s hs
    rcases Nat.eq_zero_or_pos n with rfl | hn
    · simp [S]; positivity
    have hn0 : (0 : ℝ) < n := by exact_mod_cast hn
    calc ‖S n s‖ ≤ ∑ k ∈ Finset.range n, ‖((L / n : ℝ) : ℂ) * (g (k * L / n) *
          Complex.exp (-(((k * L / n : ℝ) : ℂ) * s)))‖ := norm_sum_le _ _
      _ ≤ ∑ k ∈ Finset.range n, L / n * M' := by
          refine Finset.sum_le_sum fun k hk => ?_
          have hk : (k : ℝ) < n := by exact_mod_cast Finset.mem_range.mp hk
          have hu : (k : ℝ) * L / n ∈ Icc (0 : ℝ) L := ⟨by positivity, by
            rw [div_le_iff₀ hn0]; nlinarith⟩
          rw [norm_mul, norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_pos (div_pos hL hn0)]
          have h1 := hM' _ hu
          have h2 := hs ((k : ℝ) * L / n) hu.1
          calc L / n * (‖g (k * L / n)‖ * ‖Complex.exp (-(((k * L / n : ℝ) : ℂ) * s))‖)
              ≤ L / n * (M' * 1) := by gcongr
            _ = L / n * M' := by ring
      _ = L * M' := by rw [Finset.sum_const, Finset.card_range, nsmul_eq_mul]; field_simp
  have hline : ∀ t : ℝ, ∀ u : ℝ, 0 ≤ u → ‖Complex.exp (-((u : ℂ) * (1 / 2 + t * I)))‖ ≤ 1 := by
    intro t u hu
    exact norm_cexp_line (t := t) hu
  have hone : ∀ u : ℝ, 0 ≤ u → ‖Complex.exp (-((u : ℂ) * 1))‖ ≤ 1 := by
    intro u hu
    rw [mul_one, Complex.norm_exp]; simp; linarith
  -- dominated convergence
  set F : ℕ → ℝ → ℝ := fun n t =>
    ‖(1 + riemannZeta (1 / 2 + t * I) * (S n (1 / 2 + t * I) - S n 1)) / (1 / 2 + t * I)‖ ^ 2
  set f : ℝ → ℝ := fun t =>
    ‖(1 + riemannZeta (1 / 2 + t * I) * A (1 / 2 + t * I)) / (1 / 2 + t * I)‖ ^ 2
  set bound : ℝ → ℝ := fun t => 2 * ‖1 / ((1 / 2 : ℂ) + t * I)‖ ^ 2 +
    2 * (2 * L * M') ^ 2 * ‖riemannZeta (1 / 2 + t * I) / ((1 / 2 : ℂ) + t * I)‖ ^ 2
  have hbound_int : Integrable bound := by
    have h1 := (memLp_two_iff_integrable_sq_norm memLp_inv_line.aestronglyMeasurable).mp memLp_inv_line
    have h2 := (memLp_two_iff_integrable_sq_norm (memLp_zeta_line hP h215).aestronglyMeasurable).mp
      (memLp_zeta_line hP h215)
    exact (h1.const_mul 2).add (h2.const_mul _)
  have hmeas : ∀ n, AEStronglyMeasurable (F n) volume := by
    intro n
    refine Continuous.aestronglyMeasurable ?_
    refine (Continuous.div ?_ ?_ ?_).norm.pow 2
    · refine continuous_const.add (continuous_zeta_line.mul ?_)
      refine Continuous.sub ?_ continuous_const
      simp only [S]
      refine continuous_finsetSum _ fun k _ => continuous_const.mul (continuous_const.mul ?_)
      exact Complex.continuous_exp.comp ((continuous_const.mul
        (continuous_const.add (Complex.continuous_ofReal.mul continuous_const))).neg)
    · exact continuous_const.add (Complex.continuous_ofReal.mul continuous_const)
    · intro t h; have := congrArg Complex.re h; simp at this
  have hbd : ∀ n, ∀ᵐ t ∂volume, ‖F n t‖ ≤ bound t := by
    intro n
    refine Filter.Eventually.of_forall fun t => ?_
    have hne : (1 / 2 : ℂ) + t * I ≠ 0 := fun h => by have := congrArg Complex.re h; simp at this
    have hb1 := hSb n (1 / 2 + t * I) (hline t)
    have hb2 := hSb n 1 hone
    have hdiff : ‖S n (1 / 2 + t * I) - S n 1‖ ≤ 2 * L * M' := by
      have := norm_sub_le (S n (1 / 2 + t * I)) (S n 1); linarith
    simp only [F, bound]
    rw [Real.norm_eq_abs, abs_of_nonneg (sq_nonneg _)]
    have e : (1 + riemannZeta (1 / 2 + t * I) * (S n (1 / 2 + t * I) - S n 1)) / (1 / 2 + t * I) =
        1 / ((1 / 2 : ℂ) + t * I) + riemannZeta (1 / 2 + t * I) / ((1 / 2 : ℂ) + t * I) *
          (S n (1 / 2 + t * I) - S n 1) := by
      field_simp
    rw [e]
    have hx := norm_add_le (1 / ((1 / 2 : ℂ) + t * I)) (riemannZeta (1 / 2 + t * I) / ((1 / 2 : ℂ) + t * I) *
      (S n (1 / 2 + t * I) - S n 1))
    rw [norm_mul] at hx
    set a := ‖1 / ((1 / 2 : ℂ) + t * I)‖
    set b := ‖riemannZeta (1 / 2 + t * I) / ((1 / 2 : ℂ) + t * I)‖
    have ha : 0 ≤ a := norm_nonneg _
    have hb : 0 ≤ b := norm_nonneg _
    have h3 : b * ‖S n (1 / 2 + t * I) - S n 1‖ ≤ b * (2 * L * M') := mul_le_mul_of_nonneg_left hdiff hb
    have h0 : 0 ≤ ‖1 / ((1 / 2 : ℂ) + t * I) + riemannZeta (1 / 2 + t * I) / ((1 / 2 : ℂ) + t * I) *
      (S n (1 / 2 + t * I) - S n 1)‖ := norm_nonneg _
    have h4 : ‖1 / ((1 / 2 : ℂ) + t * I) + riemannZeta (1 / 2 + t * I) / ((1 / 2 : ℂ) + t * I) *
      (S n (1 / 2 + t * I) - S n 1)‖ ≤ a + b * (2 * L * M') := by linarith
    calc _ ≤ (a + b * (2 * L * M')) ^ 2 := pow_le_pow_left₀ h0 h4 2
      _ ≤ 2 * a ^ 2 + 2 * (2 * L * M') ^ 2 * b ^ 2 := by nlinarith [sq_nonneg (a - b * (2 * L * M'))]
  have hlim : ∀ᵐ t ∂volume, Tendsto (fun n => F n t) atTop (𝓝 (f t)) := by
    refine Filter.Eventually.of_forall fun t => ?_
    have h1 := (hconv (1 / 2 + t * I)).sub (hconv 1)
    rw [hA1, sub_neg_eq_add] at h1
    have e : (∫ u in (0 : ℝ)..L, g u * Complex.exp (-(u * (1 / 2 + t * I)))) + c₀ = A (1 / 2 + t * I) := by
      simp only [A]; ring
    rw [e] at h1
    have h2 := (tendsto_const_nhds (x := (1 : ℂ))).add
      ((tendsto_const_nhds (x := riemannZeta (1 / 2 + t * I))).mul h1)
    have h3 := h2.div_const ((1 / 2 : ℂ) + t * I)
    exact (h3.norm.pow 2)
  have hT := tendsto_integral_of_dominated_convergence bound hmeas hbound_int hbd hlim
  have hT' := hT.const_mul (2 * Real.pi)⁻¹
  refine ge_of_tendsto hT' ?_
  filter_upwards [Filter.eventually_ge_atTop 1] with n hn
  exact hfin n hn

end Riemann.Deep.LeanNBRate

