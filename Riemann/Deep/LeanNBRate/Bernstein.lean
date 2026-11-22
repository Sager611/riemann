import Riemann.Deep.LeanNBRate.MLBound
import Riemann.Deep.LeanNBRate.ExpSumOps

/-!
# Bernstein inequality for exponential sums on `[0, ℓ]`

For `r` with nodes in `[−w, w]`, `η > 0` with `(4A₀/η)(M+1)² ≤ ℓ`:
`sup_{[0,ℓ]} |r′| ≤ (2c)^{M−1} C_B (M+1)³ ((w+η)²/η) sup_{[0,ℓ]} |r|`
(the multiplier lemma with `G = iθ·cutoff(−w, w, η)`).
-/

open MeasureTheory Filter Topology Set Complex

noncomputable section

namespace Riemann.Deep.LeanNBRate

/-- (Ported from zeta23 `Taper/GevreyRamps.lean`, Apache-2.0.) `iteratedDeriv` commutes with `ℝ → ℂ`. -/
theorem iteratedDeriv_ofReal_comp {f : ℝ → ℝ} (hf : ContDiff ℝ ((⊤ : ℕ∞) : WithTop ℕ∞) f) (k : ℕ) :
    iteratedDeriv k (fun u => (f u : ℂ)) = fun u => ((iteratedDeriv k f u : ℝ) : ℂ) := by
  induction k with
  | zero => simp
  | succ k ih =>
    rw [iteratedDeriv_succ, ih, iteratedDeriv_succ]
    funext u
    exact (((hf.differentiable_iteratedDeriv k
      (WithTop.coe_lt_coe.mpr (ENat.coe_lt_top k))) u).hasDerivAt.ofReal_comp).deriv

/-- Jets of `θ ↦ cθ`. -/
theorem iteratedDeriv_linear (c : ℂ) (n : ℕ) (θ : ℝ) :
    iteratedDeriv n (fun θ : ℝ => c * (θ : ℂ)) θ =
      if n = 0 then c * θ else if n = 1 then c else 0 := by
  rw [iteratedDeriv_const_mul_field]
  have h := iteratedDeriv_ofReal_comp (f := id) contDiff_id n
  show c * iteratedDeriv n (fun u : ℝ => ((id u : ℝ) : ℂ)) θ = _
  rw [h]
  simp only [iteratedDeriv_id]
  split_ifs <;> simp

/-- The Bernstein multiplier `G(θ) = iθ·cutoff(−w, w, η)(θ)`. -/
def bernG (w η : ℝ) (θ : ℝ) : ℂ := I * θ * (cutoff (-w) w η θ : ℂ)

theorem cutoff_hasCompactSupport {a b η : ℝ} (hη : 0 < η) :
    HasCompactSupport (cutoff a b η) := by
  refine HasCompactSupport.intro (K := Icc (a - η) (b + η)) isCompact_Icc fun θ hθ => ?_
  apply cutoff_eq_zero hη
  rw [mem_Icc, not_and_or, not_le, not_le] at hθ
  rcases hθ with h | h
  · left; exact h.le
  · right; exact h.le

theorem bernG_isTest {w η : ℝ} (hη : 0 < η) : IsTest (bernG w η) := by
  have hχ := (cutoff_gevB (-w) w η hη).smooth
  refine ⟨?_, ?_⟩
  · unfold bernG
    exact (contDiff_const.mul ofRealCLM.contDiff).mul (ofRealCLM.contDiff.comp hχ)
  · unfold bernG
    have : HasCompactSupport (fun θ : ℝ => (cutoff (-w) w η θ : ℂ)) :=
      (cutoff_hasCompactSupport hη).comp_left ofReal_zero
    exact this.mul_left

/-- Exact jets: `G^{(n)} = iθ χ^{(n)} + n i χ^{(n−1)}`. -/
theorem iteratedDeriv_bernG {w η : ℝ} (hη : 0 < η) (n : ℕ) (θ : ℝ) :
    iteratedDeriv n (bernG w η) θ =
      I * θ * ((iteratedDeriv n (cutoff (-w) w η) θ : ℝ) : ℂ) +
        (n : ℂ) * I * ((iteratedDeriv (n - 1) (cutoff (-w) w η) θ : ℝ) : ℂ) := by
  have hχ := (cutoff_gevB (-w) w η hη).smooth
  have hf : ContDiff ℝ (⊤ : ℕ∞) (fun θ : ℝ => I * (θ : ℂ)) := contDiff_const.mul ofRealCLM.contDiff
  have hg : ContDiff ℝ (⊤ : ℕ∞) (fun θ : ℝ => (cutoff (-w) w η θ : ℂ)) := ofRealCLM.contDiff.comp hχ
  have e : bernG w η = (fun θ : ℝ => I * (θ : ℂ)) * (fun θ : ℝ => (cutoff (-w) w η θ : ℂ)) := rfl
  rw [e, iteratedDeriv_mul (hf.of_le (WithTop.coe_le_coe.mpr le_top)).contDiffAt
    (hg.of_le (WithTop.coe_le_coe.mpr le_top)).contDiffAt]
  simp_rw [iteratedDeriv_ofReal_comp hχ, iteratedDeriv_linear]
  rcases n with _ | m
  · simp
  · rw [Finset.sum_range_succ', Finset.sum_range_succ']
    have hz : ∀ i ∈ Finset.range m, ((m + 1).choose (i + 1 + 1) : ℂ) *
        (if i + 1 + 1 = 0 then I * θ else if i + 1 + 1 = 1 then I else 0) *
        ((iteratedDeriv (m + 1 - (i + 1 + 1)) (cutoff (-w) w η) θ : ℝ) : ℂ) = 0 := by
      intro i _; simp
    rw [Finset.sum_eq_zero hz]
    simp only [zero_add, Nat.choose_zero_right, Nat.cast_one, one_mul, ite_true, Nat.choose_one_right,
      Nat.add_sub_cancel, Nat.sub_zero, one_ne_zero, and_false, ite_false]
    push_cast
    ring

/-- At points of `[−w, w]` the jets of `bernG` are those of `iθ`. -/
theorem iteratedDeriv_bernG_of_mem {w η δ : ℝ} (hη : 0 < η) (hδ : δ ∈ Icc (-w) w) (l : ℕ) :
    iteratedDeriv l (bernG w η) δ = if l = 0 then I * δ else if l = 1 then I else 0 := by
  have hev : bernG w η =ᶠ[𝓝 δ] fun θ : ℝ => I * (θ : ℂ) := by
    filter_upwards [cutoff_eventuallyEq_one hη hδ] with θ hθ
    simp only [bernG, hθ, ofReal_one, mul_one]
  rw [hev.iteratedDeriv_eq, iteratedDeriv_linear]

/-- `Λ_r(e^{iθx} bernG) = r′(x)` when all nodes lie in `[−w, w]`. -/
theorem ExpSum.Lam_bernG (r : ExpSum) {w η : ℝ} (hη : 0 < η) (hr : ∀ j, r.node j ∈ Icc (-w) w)
    (x : ℝ) : r.Lam x (bernG w η) = (r.dmod 0).eval x := by
  rw [ExpSum.dmod_eval_eq]
  unfold ExpSum.Lam
  refine Finset.sum_congr rfl fun j _ => Finset.sum_congr rfl fun k _ => ?_
  have hinner : ∑ l ∈ Finset.range (k + 1), (Nat.choose k l : ℂ) * (I * x) ^ (k - l) *
      iteratedDeriv l (bernG w η) (r.node j) =
      (I * x) ^ k * (I * r.node j) + (k : ℂ) * (I * x) ^ (k - 1) * I := by
    simp_rw [iteratedDeriv_bernG_of_mem hη (hr j)]
    rcases k with _ | m
    · simp
    · rw [Finset.sum_range_succ', Finset.sum_range_succ']
      have hz : ∀ i ∈ Finset.range m, ((m + 1).choose (i + 1 + 1) : ℂ) * (I * x) ^ (m + 1 - (i + 1 + 1)) *
          (if i + 1 + 1 = 0 then I * r.node j else if i + 1 + 1 = 1 then I else 0) = 0 := by
        intro i _; simp
      rw [Finset.sum_eq_zero hz]
      simp only [zero_add, Nat.choose_zero_right, Nat.cast_one, one_mul, ite_true,
        Nat.choose_one_right, Nat.add_sub_cancel, Nat.sub_zero, one_ne_zero, and_false,
        ite_false]
      push_cast
      ring
  rw [hinner]
  push_cast
  ring


theorem iteratedDeriv_cutoff_eq_zero {a b η θ : ℝ} (hη : 0 < η) (hθ : θ ∉ Icc (a - η) (b + η))
    (n : ℕ) : iteratedDeriv n (cutoff a b η) θ = 0 := by
  have hev : cutoff a b η =ᶠ[𝓝 θ] fun _ => 0 := by
    rw [mem_Icc, not_and_or, not_le, not_le] at hθ
    rcases hθ with h | h
    · filter_upwards [Iio_mem_nhds h] with x hx
      exact cutoff_eq_zero hη (Or.inl (le_of_lt hx))
    · filter_upwards [Ioi_mem_nhds h] with x hx
      exact cutoff_eq_zero hη (Or.inr (le_of_lt hx))
  rw [hev.iteratedDeriv_eq]
  rcases n with _ | n
  · simp
  · rw [iteratedDeriv_succ']; simp

theorem twelve_le_A₀ : 12 ≤ A₀ := by
  unfold A₀
  rw [le_div_iff₀ (Real.exp_pos 1)]
  have := Real.exp_one_lt_d9
  linarith

theorem pow_two_mul_le {n : ℕ} : ((n - 1 : ℕ) : ℝ) ^ (2 * (n - 1)) ≤ (n : ℝ) ^ (2 * n) := by
  rcases n with _ | n
  · simp
  · simp only [Nat.add_sub_cancel]
    calc (n : ℝ) ^ (2 * n) ≤ ((n + 1 : ℕ) : ℝ) ^ (2 * n) := by gcongr; push_cast; linarith
      _ ≤ ((n + 1 : ℕ) : ℝ) ^ (2 * (n + 1)) := by
          apply pow_le_pow_right₀ (by push_cast; linarith) (by omega)

/-- `∫ |G^{(n)}| ≤ 2 s² C_χ (1+n) A_χⁿ n^{2n}`, `s = w + η`, `C_χ = B₀²`, `A_χ = 4A₀/η`. -/
theorem integral_norm_iteratedDeriv_bernG_le {w η : ℝ} (hw : 0 ≤ w) (hη : 0 < η) (n : ℕ) :
    ∫ θ : ℝ, ‖iteratedDeriv n (bernG w η) θ‖ ≤
      2 * (w + η) ^ 2 * (B₀ * B₀) * (1 + n) * (2 * (2 / η * A₀)) ^ n * (n : ℝ) ^ (2 * n) := by
  set s := w + η
  set Cχ := B₀ * B₀
  set Aχ := 2 * (2 / η * A₀)
  have hgev := cutoff_gevB (-w) w η hη
  have hs : 0 < s := by positivity
  have hCχ : 0 < Cχ := by have := B₀_pos; positivity
  have hAχ : 0 < Aχ := by have := A₀_pos; positivity
  have hsA : 1 ≤ s * Aχ := by
    have h12 := twelve_le_A₀
    have : s * Aχ = 4 * A₀ * (s / η) := by simp only [Aχ]; ring
    rw [this]
    have : 1 ≤ s / η := by rw [le_div_iff₀ hη]; linarith
    nlinarith
  -- pointwise bound on the support
  have hpt : ∀ θ ∈ Icc (-s) s, ‖iteratedDeriv n (bernG w η) θ‖ ≤
      (1 + n) * s * Cχ * Aχ ^ n * (n : ℝ) ^ (2 * n) := by
    intro θ hθ
    rw [iteratedDeriv_bernG hη]
    have hθs : |θ| ≤ s := abs_le.mpr hθ
    have h1 := hgev.bound n θ
    have h2 := hgev.bound (n - 1) θ
    have hA1 : (n : ℝ) * (Aχ ^ (n - 1) * ((n - 1 : ℕ) : ℝ) ^ (2 * (n - 1))) ≤
        (n : ℝ) * (s * Aχ ^ n * (n : ℝ) ^ (2 * n)) := by
      rcases n with _ | m
      · simp
      · simp only [Nat.add_sub_cancel]
        have h0 : Aχ ^ m ≤ s * Aχ ^ (m + 1) := by rw [pow_succ]; nlinarith [pow_pos hAχ m]
        have h3 : ((m : ℕ) : ℝ) ^ (2 * m) ≤ ((m + 1 : ℕ) : ℝ) ^ (2 * (m + 1)) := by
          have := pow_two_mul_le (n := m + 1); simpa using this
        have h4 : Aχ ^ m * (m : ℝ) ^ (2 * m) ≤ (s * Aχ ^ (m + 1)) * ((m + 1 : ℕ) : ℝ) ^ (2 * (m + 1)) := by
          gcongr
        exact mul_le_mul_of_nonneg_left (by simpa [mul_assoc] using h4) (by positivity)
    have hfirst : |θ| * |iteratedDeriv n (cutoff (-w) w η) θ| ≤ s * (Cχ * Aχ ^ n * (n : ℝ) ^ (2 * n)) :=
      mul_le_mul hθs h1 (abs_nonneg _) hs.le
    have hsecond : (n : ℝ) * |iteratedDeriv (n - 1) (cutoff (-w) w η) θ| ≤
        Cχ * ((n : ℝ) * (s * Aχ ^ n * (n : ℝ) ^ (2 * n))) := by
      calc (n : ℝ) * |iteratedDeriv (n - 1) (cutoff (-w) w η) θ|
          ≤ (n : ℝ) * (Cχ * Aχ ^ (n - 1) * ((n - 1 : ℕ) : ℝ) ^ (2 * (n - 1))) :=
            mul_le_mul_of_nonneg_left h2 (by positivity)
        _ = Cχ * ((n : ℝ) * (Aχ ^ (n - 1) * ((n - 1 : ℕ) : ℝ) ^ (2 * (n - 1)))) := by ring
        _ ≤ Cχ * ((n : ℝ) * (s * Aχ ^ n * (n : ℝ) ^ (2 * n))) :=
            mul_le_mul_of_nonneg_left hA1 hCχ.le
    calc ‖I * θ * ((iteratedDeriv n (cutoff (-w) w η) θ : ℝ) : ℂ) +
          (n : ℂ) * I * ((iteratedDeriv (n - 1) (cutoff (-w) w η) θ : ℝ) : ℂ)‖
        ≤ |θ| * |iteratedDeriv n (cutoff (-w) w η) θ| +
            n * |iteratedDeriv (n - 1) (cutoff (-w) w η) θ| := by
          refine (norm_add_le _ _).trans (le_of_eq ?_)
          simp [Complex.norm_I, Complex.norm_real, Real.norm_eq_abs]
      _ ≤ s * (Cχ * Aχ ^ n * (n : ℝ) ^ (2 * n)) + Cχ * ((n : ℝ) * (s * Aχ ^ n * (n : ℝ) ^ (2 * n))) :=
          add_le_add hfirst hsecond
      _ = (1 + n) * s * Cχ * Aχ ^ n * (n : ℝ) ^ (2 * n) := by ring
  -- vanishing off the support
  have hzero : ∀ θ ∉ Icc (-s) s, iteratedDeriv n (bernG w η) θ = 0 := by
    intro θ hθ
    have hθ' : θ ∉ Icc (-w - η) (w + η) := by
      simpa [s, show -w - η = -(w + η) by ring] using hθ
    rw [iteratedDeriv_bernG hη, iteratedDeriv_cutoff_eq_zero hη hθ',
      iteratedDeriv_cutoff_eq_zero hη hθ']
    simp
  have hcont : Continuous fun θ => ‖iteratedDeriv n (bernG w η) θ‖ :=
    ((bernG_isTest (w := w) hη).iterDeriv n).continuous.norm
  rw [← setIntegral_eq_integral_of_forall_compl_eq_zero (s := Icc (-s) s)
    (fun θ hθ => by rw [hzero θ hθ, norm_zero])]
  calc ∫ θ in Icc (-s) s, ‖iteratedDeriv n (bernG w η) θ‖
      ≤ ∫ θ in Icc (-s) s, (1 + n) * s * Cχ * Aχ ^ n * (n : ℝ) ^ (2 * n) :=
        setIntegral_mono_on (hcont.continuousOn.integrableOn_compact isCompact_Icc)
          (integrableOn_const (by rw [Real.volume_Icc]; exact ENNReal.ofReal_ne_top))
          measurableSet_Icc hpt
    _ = 2 * s ^ 2 * Cχ * (1 + n) * Aχ ^ n * (n : ℝ) ^ (2 * n) := by
        rw [setIntegral_const, Real.volume_real_Icc_of_le (by linarith), smul_eq_mul]; ring

/-- The constant of the Bernstein inequality. -/
def CB : ℝ := B₀ * B₀ * (2 + 1536 * A₀ ^ 2 + 16 * A₀)

/-- **Bernstein inequality** on `[0, ℓ]` for exponential sums with nodes in `[−w, w]`. -/
theorem bernstein (hT : Nazarov1993_Turan) : ∃ c : ℝ, 0 < c ∧
    ∀ (r : ExpSum) (w η ℓ K : ℝ), 0 ≤ w → 0 < η → 0 < ℓ → 1 ≤ r.M →
      (∀ j, r.node j ∈ Icc (-w) w) → 4 * A₀ / η * ((r.M : ℝ) + 1) ^ 2 ≤ ℓ →
      (∀ u ∈ Icc 0 ℓ, ‖r.eval u‖ ≤ K) → ∀ x ∈ Icc 0 ℓ,
        ‖(r.dmod 0).eval x‖ ≤ (2 * c) ^ (r.M - 1) * CB * ((r.M : ℝ) + 1) ^ 3 *
          ((w + η) ^ 2 / η) * K := by
  obtain ⟨c, hc, hml⟩ := ml_bound hT
  refine ⟨c, hc, fun r w η ℓ K hw hη hℓ hM hnodes hℓbig hK x hx => ?_⟩
  set M := r.M
  set s := w + η
  set Cχ := B₀ * B₀
  set Aχ := 2 * (2 / η * A₀)
  have hs : 0 < s := by positivity
  have hCχ : 0 < Cχ := by have := B₀_pos; positivity
  have hA0 : 0 < A₀ := A₀_pos
  have hAχ : 0 < Aχ := by positivity
  have hK0 : 0 ≤ K := (norm_nonneg _).trans (hK 0 ⟨le_rfl, hℓ.le⟩)
  have hG := bernG_isTest (w := w) hη
  rw [← r.Lam_bernG hη hnodes x]
  set R₁ := Aχ * ((M : ℝ) + 1) ^ 2
  have hR₁ : 0 < R₁ := by positivity
  have hmain := hml r hG ℓ K x (1 / η) R₁ hℓ hM hx (by positivity) hR₁ hK
  set A := fun n : ℕ => ∫ θ : ℝ, ‖iteratedDeriv n (bernG w η) θ‖
  have hA : ∀ n, A n ≤ 2 * s ^ 2 * Cχ * (1 + n) * Aχ ^ n * (n : ℝ) ^ (2 * n) :=
    fun n => integral_norm_iteratedDeriv_bernG_le hw hη n
  have hA0' : (∫ θ : ℝ, ‖bernG w η θ‖) = A 0 := by simp [A]
  rw [hA0'] at hmain
  -- the ratio X = Aχ (M+1)² / ℓ ≤ 1
  have hAχeq : Aχ = 4 * A₀ / η := by simp only [Aχ]; ring
  have hX : Aχ * ((M : ℝ) + 1) ^ 2 ≤ ℓ := by rw [hAχeq]; exact hℓbig
  -- bounds on the four terms
  have hT1 : 1 / η * A 0 ≤ 2 * s ^ 2 * Cχ / η := by
    have h := hA 0
    simp only [Nat.cast_zero, add_zero, pow_zero, mul_one, mul_zero] at h
    calc 1 / η * A 0 ≤ 1 / η * (2 * s ^ 2 * Cχ) := mul_le_mul_of_nonneg_left h (by positivity)
      _ = 2 * s ^ 2 * Cχ / η := by ring
  have hT2 : A 2 / (1 / η) ≤ 1536 * A₀ ^ 2 * s ^ 2 * Cχ / η := by
    have := hA 2
    rw [div_div_eq_mul_div, div_one]
    calc A 2 * η ≤ (2 * s ^ 2 * Cχ * (1 + ((2 : ℕ) : ℝ)) * Aχ ^ 2 * ((2 : ℕ) : ℝ) ^ (2 * 2)) * η :=
          mul_le_mul_of_nonneg_right this hη.le
      _ = 1536 * A₀ ^ 2 * s ^ 2 * Cχ / η := by push_cast; rw [hAχeq]; field_simp; ring
  have hpowX : (ℓ ^ (M - 1))⁻¹ * Aχ ^ M * ((M : ℝ) + 1) ^ (2 * M) ≤ Aχ * ((M : ℝ) + 1) ^ 2 := by
    obtain ⟨m, hm⟩ : ∃ m, M = m + 1 := ⟨M - 1, by omega⟩
    rw [hm, Nat.add_sub_cancel]
    have hXle : Aχ * ((m + 1 : ℕ) + 1 : ℝ) ^ 2 / ℓ ≤ 1 := by
      rw [div_le_one hℓ]; push_cast at hX ⊢; rw [hm] at hX; push_cast at hX; linarith
    have hX0 : 0 ≤ Aχ * ((m + 1 : ℕ) + 1 : ℝ) ^ 2 / ℓ := by positivity
    have : (ℓ ^ m)⁻¹ * Aχ ^ (m + 1) * (((m + 1 : ℕ) : ℝ) + 1) ^ (2 * (m + 1)) =
        Aχ * (((m + 1 : ℕ) : ℝ) + 1) ^ 2 * (Aχ * (((m + 1 : ℕ) : ℝ) + 1) ^ 2 / ℓ) ^ m := by
      set X : ℝ := ((m + 1 : ℕ) : ℝ) + 1
      have e1 : X ^ (2 * (m + 1)) = (X ^ 2) ^ (m + 1) := by rw [← pow_mul]
      rw [e1, div_pow]
      have hℓm : ℓ ^ m ≠ 0 := pow_ne_zero _ hℓ.ne'
      field_simp
      ring
    rw [this]
    calc Aχ * (((m + 1 : ℕ) : ℝ) + 1) ^ 2 * (Aχ * (((m + 1 : ℕ) : ℝ) + 1) ^ 2 / ℓ) ^ m
        ≤ Aχ * (((m + 1 : ℕ) : ℝ) + 1) ^ 2 * 1 := by
          gcongr; exact pow_le_one₀ hX0 hXle
      _ = _ := by ring
  have hT34 : (ℓ ^ (M - 1))⁻¹ * (R₁ * A (M - 1) + A (M - 1 + 2) / R₁) ≤
      16 * A₀ * s ^ 2 * Cχ * ((M : ℝ) + 1) ^ 3 / η := by
    have h1 := hA (M - 1)
    have h2 := hA (M - 1 + 2)
    have hM1 : ((M - 1 : ℕ) : ℝ) = (M : ℝ) - 1 := by
      rw [Nat.cast_sub hM]; simp
    have hM2 : ((M - 1 + 2 : ℕ) : ℝ) = (M : ℝ) + 1 := by
      rw [Nat.cast_add, hM1]; push_cast; ring
    have e2 : M - 1 + 2 = M + 1 := by omega
    -- R₁ A_{M−1} ≤ 2 s² Cχ M Aχ^M (M+1)^{2M}
    have hb1 : R₁ * A (M - 1) ≤ 2 * s ^ 2 * Cχ * M * (Aχ ^ M * ((M : ℝ) + 1) ^ (2 * M)) := by
      have hp : ((M - 1 : ℕ) : ℝ) ^ (2 * (M - 1)) * ((M : ℝ) + 1) ^ 2 ≤ ((M : ℝ) + 1) ^ (2 * M) := by
        have : ((M - 1 : ℕ) : ℝ) ^ (2 * (M - 1)) ≤ ((M : ℝ) + 1) ^ (2 * (M - 1)) := by
          gcongr; rw [hM1]; linarith
        calc ((M - 1 : ℕ) : ℝ) ^ (2 * (M - 1)) * ((M : ℝ) + 1) ^ 2
            ≤ ((M : ℝ) + 1) ^ (2 * (M - 1)) * ((M : ℝ) + 1) ^ 2 := by gcongr
          _ = ((M : ℝ) + 1) ^ (2 * M) := by rw [← pow_add]; congr 1; omega
      have hAp : Aχ * Aχ ^ (M - 1) = Aχ ^ M := by
        rw [← pow_succ']; congr 1; omega
      calc R₁ * A (M - 1) ≤ R₁ * (2 * s ^ 2 * Cχ * (1 + ((M - 1 : ℕ) : ℝ)) * Aχ ^ (M - 1) *
            ((M - 1 : ℕ) : ℝ) ^ (2 * (M - 1))) := by gcongr
        _ = 2 * s ^ 2 * Cχ * M * ((Aχ * Aχ ^ (M - 1)) *
            (((M - 1 : ℕ) : ℝ) ^ (2 * (M - 1)) * ((M : ℝ) + 1) ^ 2)) := by
            rw [hM1]; simp only [R₁]; ring
        _ ≤ 2 * s ^ 2 * Cχ * M * (Aχ ^ M * ((M : ℝ) + 1) ^ (2 * M)) := by
            rw [hAp]; gcongr
    -- A_{M+1}/R₁ ≤ 2 s² Cχ (M+2) Aχ^M (M+1)^{2M}
    have hb2 : A (M - 1 + 2) / R₁ ≤ 2 * s ^ 2 * Cχ * ((M : ℝ) + 2) *
        (Aχ ^ M * ((M : ℝ) + 1) ^ (2 * M)) := by
      rw [div_le_iff₀ hR₁]
      calc A (M - 1 + 2) ≤ 2 * s ^ 2 * Cχ * (1 + ((M - 1 + 2 : ℕ) : ℝ)) * Aχ ^ (M - 1 + 2) *
            ((M - 1 + 2 : ℕ) : ℝ) ^ (2 * (M - 1 + 2)) := h2
        _ = 2 * s ^ 2 * Cχ * ((M : ℝ) + 2) * (Aχ ^ M * ((M : ℝ) + 1) ^ (2 * M)) * R₁ := by
            rw [hM2, e2]; simp only [R₁]
            rw [pow_succ, show 2 * (M + 1) = 2 * M + 2 by ring, pow_add]; ring
    calc (ℓ ^ (M - 1))⁻¹ * (R₁ * A (M - 1) + A (M - 1 + 2) / R₁)
        ≤ (ℓ ^ (M - 1))⁻¹ * (2 * s ^ 2 * Cχ * M * (Aχ ^ M * ((M : ℝ) + 1) ^ (2 * M)) +
            2 * s ^ 2 * Cχ * ((M : ℝ) + 2) * (Aχ ^ M * ((M : ℝ) + 1) ^ (2 * M))) := by gcongr
      _ = 2 * s ^ 2 * Cχ * (2 * (M : ℝ) + 2) *
            ((ℓ ^ (M - 1))⁻¹ * Aχ ^ M * ((M : ℝ) + 1) ^ (2 * M)) := by ring
      _ ≤ 2 * s ^ 2 * Cχ * (2 * (M : ℝ) + 2) * (Aχ * ((M : ℝ) + 1) ^ 2) := by gcongr
      _ = 16 * A₀ * s ^ 2 * Cχ * ((M : ℝ) + 1) ^ 3 / η := by rw [hAχeq]; field_simp; ring
  have hbr : 1 / η * A 0 + A 2 / (1 / η) +
      (ℓ ^ (M - 1))⁻¹ * (R₁ * A (M - 1) + A (M - 1 + 2) / R₁) ≤
      CB * ((M : ℝ) + 1) ^ 3 * (s ^ 2 / η) := by
    have hM1 : (1 : ℝ) ≤ ((M : ℝ) + 1) ^ 3 := one_le_pow₀ (by have h0 : (0:ℝ) ≤ M := Nat.cast_nonneg M; linarith)
    calc 1 / η * A 0 + A 2 / (1 / η) + (ℓ ^ (M - 1))⁻¹ * (R₁ * A (M - 1) + A (M - 1 + 2) / R₁)
        ≤ 2 * s ^ 2 * Cχ / η + 1536 * A₀ ^ 2 * s ^ 2 * Cχ / η +
            16 * A₀ * s ^ 2 * Cχ * ((M : ℝ) + 1) ^ 3 / η := by linarith
      _ = (s ^ 2 * Cχ / η) * (2 + 1536 * A₀ ^ 2 + 16 * A₀ * ((M : ℝ) + 1) ^ 3) := by ring
      _ ≤ (s ^ 2 * Cχ / η) * ((2 + 1536 * A₀ ^ 2 + 16 * A₀) * ((M : ℝ) + 1) ^ 3) := by
          gcongr; nlinarith [sq_nonneg A₀]
      _ = CB * ((M : ℝ) + 1) ^ 3 * (s ^ 2 / η) := by simp only [CB, Cχ]; ring
  calc ‖r.Lam x (bernG w η)‖ ≤ c ^ (M - 1) * K * 2 ^ (M - 1) *
        (1 / η * A 0 + A 2 / (1 / η) + (ℓ ^ (M - 1))⁻¹ * (R₁ * A (M - 1) + A (M - 1 + 2) / R₁)) :=
        hmain
    _ ≤ c ^ (M - 1) * K * 2 ^ (M - 1) * (CB * ((M : ℝ) + 1) ^ 3 * (s ^ 2 / η)) := by gcongr
    _ = (2 * c) ^ (M - 1) * CB * ((M : ℝ) + 1) ^ 3 * (s ^ 2 / η) * K := by rw [mul_pow]; ring

end Riemann.Deep.LeanNBRate
