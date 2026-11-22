module

public import Riemann.Deep.LeanLi.LiCoeff

/-!
# LeanLi: Lemma 0 — if RH fails, `|λ_n| ≥ e^{cn}` for infinitely many `n`

A zero `ρ₀` with `Re ρ₀ < 1/2`, `Im ρ₀ > 0` gives `‖1 - 1/ρ₀‖ > 1`. Only finitely many zeros have
`‖1 - 1/ρ‖ ≥ τ > 1`; let `M` be the maximum and `S` the zeros attaining it. Then
`λ_n = A - Mⁿ Σ_S 2m_ρ cos(n arg w_ρ) + O(n² M₂ⁿ)` with `M₂ < M`, and the cosine sum (positive
coefficients) is `≥ c > 0` in absolute value infinitely often (its correlation with one of its own
frequencies has mean `≥ a/2`).
-/

@[expose] public section

open Real Filter Set Finset Topology Complex ComplexConjugate

namespace LeanLi

/-- Off RH, `|λ_n| ≥ e^{cn}` infinitely often (stated by A. Voros, arXiv math/0506326, TeX l.55–60 and
l.284–299, and E. Bombieri, J. C. Lagarias, J. Number Theory 77 (1999), Thm 1(c); proved below). -/
def LiExpGrowth : Prop :=
  ¬ RiemannHypothesis → ∃ c > (0 : ℝ), ∃ᶠ n : ℕ in atTop, Real.exp (c * n) ≤ |lam n|

/-- Telescoping: `2 sin(x/2) Σ_{n<N} cos nx = sin((N - 1/2)x) + sin(x/2)`. -/
lemma two_sin_mul_sum_range_cos (x : ℝ) (N : ℕ) :
    2 * Real.sin (x / 2) * ∑ n ∈ range N, Real.cos (n * x) =
      Real.sin ((N - 1 / 2) * x) + Real.sin (x / 2) := by
  induction N with
  | zero =>
    simp only [Finset.sum_range_zero, mul_zero, Nat.cast_zero, zero_sub]
    rw [show -(1 / 2 : ℝ) * x = -(x / 2) by ring, Real.sin_neg]; ring
  | succ N ih =>
    rw [Finset.sum_range_succ, mul_add, ih]
    have e : 2 * Real.sin (x / 2) * Real.cos (N * x) =
        Real.sin ((N + 1 / 2) * x) - Real.sin ((N - 1 / 2) * x) := by
      rw [show ((N : ℝ) + 1 / 2) * x = N * x + x / 2 by ring,
        show ((N : ℝ) - 1 / 2) * x = N * x - x / 2 by ring, Real.sin_add, Real.sin_sub]
      ring
    rw [e]
    push_cast
    rw [show ((N : ℝ) + 1 - 1 / 2) * x = (N + 1 / 2) * x by ring]
    ring

/-- `Σ_{n<N} cos nx` is bounded below uniformly in `N`. -/
lemma sum_range_cos_ge (x : ℝ) : ∃ B : ℝ, 0 ≤ B ∧ ∀ N : ℕ, -B ≤ ∑ n ∈ range N, Real.cos (n * x) := by
  by_cases hs : Real.sin (x / 2) = 0
  · refine ⟨0, le_rfl, fun N => ?_⟩
    obtain ⟨k, hk⟩ := Real.sin_eq_zero_iff.1 hs
    have h1 : ∀ n : ℕ, Real.cos (n * x) = 1 := by
      intro n
      have hx : x = 2 * ((k : ℝ) * Real.pi) := by linarith
      have : (n : ℝ) * x = ((n * k : ℤ) : ℝ) * (2 * Real.pi) := by push_cast; rw [hx]; ring
      rw [this, Real.cos_int_mul_two_pi]
    simp only [h1, Finset.sum_const, Finset.card_range, nsmul_eq_mul, mul_one, neg_zero]
    exact Nat.cast_nonneg N
  · refine ⟨1 / |Real.sin (x / 2)|, by positivity, fun N => ?_⟩
    have h := two_sin_mul_sum_range_cos x N
    have hb : |2 * Real.sin (x / 2) * ∑ n ∈ range N, Real.cos (n * x)| ≤ 2 := by
      rw [h]
      calc |Real.sin ((N - 1 / 2) * x) + Real.sin (x / 2)| ≤
            |Real.sin ((N - 1 / 2) * x)| + |Real.sin (x / 2)| := abs_add_le _ _
        _ ≤ 1 + 1 := add_le_add (abs_sin_le_one _) (abs_sin_le_one _)
        _ = 2 := by norm_num
    rw [abs_mul, abs_mul, abs_two] at hb
    have hpos : 0 < |Real.sin (x / 2)| := abs_pos.2 hs
    have h2 : |∑ n ∈ range N, Real.cos (n * x)| ≤ 1 / |Real.sin (x / 2)| := by
      rw [le_div_iff₀ hpos]; nlinarith [abs_nonneg (∑ n ∈ range N, Real.cos (n * x))]
    linarith [neg_abs_le (∑ n ∈ range N, Real.cos (n * x))]

/-- **Cosine sums with positive coefficients do not tend to `0`.** -/
lemma osc {ι : Type*} (s : Finset ι) (hs : s.Nonempty) (a φ : ι → ℝ) (ha : ∀ i ∈ s, 0 < a i) :
    ∃ c > (0 : ℝ), ∃ᶠ n : ℕ in atTop, c ≤ |∑ i ∈ s, a i * Real.cos (n * φ i)| := by
  classical
  obtain ⟨i₀, hi₀⟩ := hs
  have ha0 := ha i₀ hi₀
  set f : ℕ → ℝ := fun n => ∑ i ∈ s, a i * Real.cos (n * φ i) with hf
  set A := ∑ i ∈ s, a i with hA
  have hfA : ∀ n, |f n| ≤ A := by
    intro n
    refine (Finset.abs_sum_le_sum_abs _ _).trans (Finset.sum_le_sum fun i hi => ?_)
    rw [abs_mul, abs_of_pos (ha i hi)]
    exact mul_le_of_le_one_right (ha i hi).le (abs_cos_le_one _)
  choose B hB0 hB using fun x : ℝ => sum_range_cos_ge x
  set K := ∑ j ∈ s, a j * (B (φ j - φ i₀) + B (φ j + φ i₀)) / 2 with hK
  -- correlation with the frequency `φ i₀`
  have hlow : ∀ N : ℕ, a i₀ * N / 2 - K ≤ ∑ n ∈ range N, f n * Real.cos (n * φ i₀) := by
    intro N
    have e : ∑ n ∈ range N, f n * Real.cos (n * φ i₀) = ∑ j ∈ s, a j *
        ((∑ n ∈ range N, Real.cos (n * (φ j - φ i₀)) + ∑ n ∈ range N, Real.cos (n * (φ j + φ i₀))) / 2) := by
      simp only [hf, Finset.sum_mul]
      rw [Finset.sum_comm]
      refine Finset.sum_congr rfl fun j _ => ?_
      rw [← Finset.sum_add_distrib, Finset.sum_div, Finset.mul_sum]
      refine Finset.sum_congr rfl fun n _ => ?_
      rw [show (n : ℝ) * (φ j - φ i₀) = n * φ j - n * φ i₀ by ring,
        show (n : ℝ) * (φ j + φ i₀) = n * φ j + n * φ i₀ by ring, Real.cos_sub, Real.cos_add]
      ring
    rw [e]
    have hterm : ∀ j ∈ s, a j * (if j = i₀ then (N : ℝ) / 2 else 0) - a j * (B (φ j - φ i₀) + B (φ j + φ i₀)) / 2 ≤
        a j * ((∑ n ∈ range N, Real.cos (n * (φ j - φ i₀)) + ∑ n ∈ range N, Real.cos (n * (φ j + φ i₀))) / 2) := by
      intro j hj
      have h2 := hB (φ j + φ i₀) N
      have haj := (ha j hj).le
      split_ifs with hji
      · subst hji
        have h1 : ∑ n ∈ range N, Real.cos (n * (φ j - φ j)) = N := by simp
        rw [h1]
        have := hB0 (φ j - φ j)
        nlinarith
      · have h1 := hB (φ j - φ i₀) N
        nlinarith
    have hsum := Finset.sum_le_sum hterm
    rw [Finset.sum_sub_distrib] at hsum
    have e2 : ∑ j ∈ s, a j * (if j = i₀ then (N : ℝ) / 2 else 0) = a i₀ * N / 2 := by
      simp only [mul_ite, mul_zero]
      rw [Finset.sum_ite_eq' s i₀, if_pos hi₀]; ring
    rw [e2] at hsum
    have e3 : ∑ j ∈ s, a j * (B (φ j - φ i₀) + B (φ j + φ i₀)) / 2 = K := rfl
    linarith
  refine ⟨a i₀ / 4, by positivity, ?_⟩
  by_contra hcon
  rw [Filter.not_frequently] at hcon
  obtain ⟨n₀, hn₀⟩ := Filter.eventually_atTop.1 hcon
  set m := ⌈4 * (n₀ * A + K) / a i₀⌉₊ + 1 with hm
  have hup : ∑ n ∈ range (n₀ + m), f n * Real.cos (n * φ i₀) ≤ n₀ * A + m * (a i₀ / 4) := by
    rw [Finset.sum_range_add]
    have h1 : ∑ n ∈ range n₀, f n * Real.cos (n * φ i₀) ≤ n₀ * A := by
      have : ∀ n ∈ range n₀, f n * Real.cos (n * φ i₀) ≤ A := fun n _ => by
        have := hfA n
        have hc := abs_cos_le_one (n * φ i₀)
        calc f n * Real.cos (n * φ i₀) ≤ |f n * Real.cos (n * φ i₀)| := le_abs_self _
          _ = |f n| * |Real.cos (n * φ i₀)| := abs_mul _ _
          _ ≤ A * 1 := mul_le_mul this hc (abs_nonneg _) (le_trans (abs_nonneg _) this)
          _ = A := mul_one _
      have := Finset.sum_le_sum this
      simpa using this
    have h2 : ∑ x ∈ range m, f (n₀ + x) * Real.cos (((n₀ + x : ℕ) : ℝ) * φ i₀) ≤ m * (a i₀ / 4) := by
      have : ∀ x ∈ range m, f (n₀ + x) * Real.cos (((n₀ + x : ℕ) : ℝ) * φ i₀) ≤ a i₀ / 4 := fun x _ => by
        have hlt : |f (n₀ + x)| < a i₀ / 4 := by
          have := hn₀ (n₀ + x) (by omega)
          push Not at this
          exact this
        have hc := abs_cos_le_one (((n₀ + x : ℕ) : ℝ) * φ i₀)
        calc f (n₀ + x) * Real.cos (((n₀ + x : ℕ) : ℝ) * φ i₀) ≤ |f (n₀ + x) * Real.cos (((n₀ + x : ℕ) : ℝ) * φ i₀)| :=
              le_abs_self _
          _ = |f (n₀ + x)| * |Real.cos (((n₀ + x : ℕ) : ℝ) * φ i₀)| := abs_mul _ _
          _ ≤ (a i₀ / 4) * 1 := mul_le_mul hlt.le hc (abs_nonneg _) (by positivity)
          _ = a i₀ / 4 := mul_one _
      have := Finset.sum_le_sum this
      simpa using this
    linarith
  have hl := hlow (n₀ + m)
  have hmK : 4 * (n₀ * A + K) / a i₀ < m := by
    have := Nat.le_ceil (4 * (n₀ * A + K) / a i₀)
    rw [hm]; push_cast; linarith
  rw [div_lt_iff₀ ha0] at hmK
  push_cast at hl
  nlinarith

/-- The Li weight `w_ρ = 1 - 1/ρ` satisfies `|w|² = 1 + (1 - 2 Re ρ)/|ρ|²`. -/
lemma norm_w_sq {ρ : ℂ} (hρ : ρ ≠ 0) : ‖1 - 1 / ρ‖ ^ 2 = 1 + (1 - 2 * ρ.re) / ‖ρ‖ ^ 2 := by
  have h : 1 - 1 / ρ = (ρ - 1) / ρ := by field_simp
  have hn : (0 : ℝ) < ‖ρ‖ ^ 2 := by have := norm_pos_iff.2 hρ; positivity
  rw [h, norm_div, div_pow]
  have e1 : ‖ρ - 1‖ ^ 2 = ‖ρ‖ ^ 2 - 2 * ρ.re + 1 := by
    rw [Complex.sq_norm, Complex.sq_norm, Complex.normSq_apply, Complex.normSq_apply]; simp; ring
  rw [e1]; field_simp; ring

/-- Bound on the term of an upper zero with `‖1 - 1/ρ‖ ≤ K`. -/
lemma abs_upperTerm_le_K {n : ℕ} {K : ℝ} (hK : 1 ≤ K) {ρ : ℂ} (hρ : ρ ∈ NontrivialZeros) (him : 0 < ρ.im)
    (hwK : ‖1 - 1 / ρ‖ ≤ K) :
    |upperTerm n ρ| ≤ 2 * (n + n ^ 2 * K ^ n) * ((mult ρ : ℝ) / ρ.im ^ 2) := by
  have hre1 := re_lt_one_of_mem hρ
  have hρ0 : ρ ≠ 0 := fun h => by rw [h] at him; simp at him
  have hnorm : ρ.im ≤ ‖ρ‖ := (le_abs_self _).trans (Complex.abs_im_le_norm ρ)
  have hnpos : 0 < ‖ρ‖ := norm_pos_iff.2 hρ0
  set w : ℂ := 1 - 1 / ρ with hw
  have h1w : 1 - w = 1 / ρ := by rw [hw]; ring
  have hre_inv : (1 / ρ).re = ρ.re / ‖ρ‖ ^ 2 := by
    rw [one_div, Complex.inv_re, Complex.normSq_eq_norm_sq]
  have hnorm_inv : ‖1 / ρ‖ = 1 / ‖ρ‖ := by rw [norm_div, norm_one]
  have hwle : max 1 ‖w‖ ≤ K := max_le hK hwK
  have hmain := norm_one_sub_pow_sub_le w n
  rw [h1w, hnorm_inv] at hmain
  have hre : |(1 - w ^ n).re| ≤ n * (1 / ‖ρ‖ ^ 2) + n ^ 2 * (1 / ‖ρ‖) ^ 2 * K ^ n := by
    have h2 : |((1 - w ^ n) - n * (1 / ρ)).re| ≤ n ^ 2 * (1 / ‖ρ‖) ^ 2 * K ^ n := by
      refine (Complex.abs_re_le_norm _).trans (hmain.trans ?_)
      gcongr
    have h3 : |((n : ℂ) * (1 / ρ)).re| ≤ n * (1 / ‖ρ‖ ^ 2) := by
      have : ((n : ℂ) * (1 / ρ)).re = n * (ρ.re / ‖ρ‖ ^ 2) := by
        rw [← hre_inv]; simp [Complex.mul_re]
      rw [this, abs_mul, Nat.abs_cast]
      gcongr
      rw [abs_of_pos (by have := re_pos_of_mem hρ; positivity)]
      exact div_le_div_of_nonneg_right hre1.le (by positivity)
    have := abs_add_le ((1 - w ^ n) - n * (1 / ρ)).re ((n : ℂ) * (1 / ρ)).re
    rw [← Complex.add_re, sub_add_cancel] at this
    linarith
  have hγ : 1 / ‖ρ‖ ^ 2 ≤ 1 / ρ.im ^ 2 := by gcongr
  unfold upperTerm
  rw [← hw, abs_mul, abs_of_nonneg (by positivity : (0 : ℝ) ≤ 2 * mult ρ)]
  calc 2 * (mult ρ : ℝ) * |(1 - w ^ n).re|
      ≤ 2 * (mult ρ : ℝ) * (n * (1 / ρ.im ^ 2) + n ^ 2 * (1 / ρ.im ^ 2) * K ^ n) := by
        gcongr
        refine hre.trans ?_
        have : (1 / ‖ρ‖) ^ 2 = 1 / ‖ρ‖ ^ 2 := by ring
        rw [this]
        gcongr
    _ = 2 * (n + n ^ 2 * K ^ n) * ((mult ρ : ℝ) / ρ.im ^ 2) := by ring

/-- A zero with `Re ρ < 1/2` and `Im ρ > 0` when RH fails. -/
lemma exists_zero_left (L : LitCore) (hRH : ¬ RiemannHypothesis) :
    ∃ ρ ∈ NontrivialZeros, ρ.re < 1 / 2 ∧ 0 < ρ.im := by
  have h1 : ∃ ρ ∈ NontrivialZeros, ρ.re ≠ 1 / 2 := by
    by_contra h
    push Not at h
    exact hRH fun s hs htriv hs1 => h s ⟨hs, htriv, hs1⟩
  obtain ⟨σ, hσ, hσre⟩ := h1
  obtain ⟨τ, hτ, hτre⟩ : ∃ τ ∈ NontrivialZeros, τ.re < 1 / 2 := by
    rcases lt_or_gt_of_ne hσre with h | h
    · exact ⟨σ, hσ, h⟩
    · exact ⟨1 - σ, one_sub_mem hσ, by simp; linarith⟩
  have hτim := im_ne_zero L.zeta_unit_interval hτ
  rcases lt_or_gt_of_ne hτim with h | h
  · exact ⟨conj τ, conj_mem hτ, by simpa using hτre, by simp; linarith⟩
  · exact ⟨τ, hτ, hτre, h⟩

set_option maxHeartbeats 4000000 in
/-- **Lemma 0** (proved): if RH fails, `|λ_n| ≥ e^{cn}` for infinitely many `n`. -/
theorem lam_exp_growth (L : LitCore) : LiExpGrowth := by
  classical
  intro hRH
  obtain ⟨ρ₀, hρ₀, hre₀, him₀⟩ := exists_zero_left L hRH
  set wn : ℂ → ℝ := fun ρ => ‖1 - 1 / ρ‖ with hwn
  have hρ₀0 : ρ₀ ≠ 0 := fun h => by rw [h] at him₀; simp at him₀
  have hw₀ : 1 < wn ρ₀ := by
    have h := norm_w_sq hρ₀0
    have hn : 0 < ‖ρ₀‖ ^ 2 := by have := norm_pos_iff.2 hρ₀0; positivity
    have hpos : 0 < (1 - 2 * ρ₀.re) / ‖ρ₀‖ ^ 2 := div_pos (by linarith) hn
    have h2 : 1 < wn ρ₀ ^ 2 := by simp only [hwn]; rw [h]; linarith
    have h3 : 0 ≤ wn ρ₀ := norm_nonneg _
    nlinarith
  set τ := (1 + wn ρ₀) / 2 with hτ
  have hτ1 : 1 < τ := by linarith
  have hτw : τ < wn ρ₀ := by linarith
  have hbd : ∀ ρ ∈ NontrivialZeros, τ ≤ wn ρ → ‖ρ‖ ≤ 1 / (τ ^ 2 - 1) + 1 := by
    intro ρ hρ hτρ
    have hre := re_pos_of_mem hρ
    have hρ0 : ρ ≠ 0 := fun h => by rw [h] at hre; simp at hre
    have h := norm_w_sq hρ0
    have hn : 0 < ‖ρ‖ ^ 2 := by have := norm_pos_iff.2 hρ0; positivity
    have h2 : τ ^ 2 ≤ wn ρ ^ 2 := pow_le_pow_left₀ (by linarith) hτρ 2
    simp only [hwn] at h2
    rw [h] at h2
    have h3 : (1 - 2 * ρ.re) / ‖ρ‖ ^ 2 ≤ 1 / ‖ρ‖ ^ 2 := div_le_div_of_nonneg_right (by linarith) hn.le
    have h5 : 0 < τ ^ 2 - 1 := by nlinarith
    have h4 : τ ^ 2 - 1 ≤ 1 / ‖ρ‖ ^ 2 := by linarith
    have h6 : ‖ρ‖ ^ 2 ≤ 1 / (τ ^ 2 - 1) := by
      rw [le_div_iff₀ h5]; rw [le_div_iff₀ hn] at h4; linarith
    nlinarith [norm_nonneg ρ]
  have hFfin : {ρ : ℂ | ρ ∈ NontrivialZeros ∧ 0 < ρ.im ∧ τ ≤ wn ρ}.Finite :=
    (finite_norm_le (1 / (τ ^ 2 - 1) + 1)).subset fun ρ hρ => ⟨hρ.1, hbd ρ hρ.1 hρ.2.2⟩
  set Fs := hFfin.toFinset with hFs
  have hmemF : ∀ ρ : ℂ, ρ ∈ Fs ↔ ρ ∈ NontrivialZeros ∧ 0 < ρ.im ∧ τ ≤ wn ρ := fun ρ => by
    rw [hFs, Set.Finite.mem_toFinset]; rfl
  have hρ₀F : ρ₀ ∈ Fs := (hmemF ρ₀).2 ⟨hρ₀, him₀, hτw.le⟩
  have hFne : Fs.Nonempty := ⟨ρ₀, hρ₀F⟩
  set M := Fs.sup' hFne wn with hM
  have hMw₀ : wn ρ₀ ≤ M := Finset.le_sup' wn hρ₀F
  have hM1 : 1 < M := lt_of_lt_of_le hw₀ hMw₀
  have hM0 : 0 < M := by linarith
  obtain ⟨ρm, hρmF, hρmM⟩ := Finset.exists_mem_eq_sup' hFne wn
  set g : ℂ → ℝ := fun ρ => if wn ρ < M then wn ρ else τ with hg
  set M₂ := Fs.sup' hFne g with hM₂def
  have hgM : ∀ ρ ∈ Fs, g ρ < M := fun ρ _ => by
    simp only [hg]; split_ifs with h
    · exact h
    · linarith
  have hM₂ : M₂ < M := (Finset.sup'_lt_iff hFne).2 hgM
  have hτM₂ : τ ≤ M₂ := by
    refine le_trans ?_ (Finset.le_sup' g hρ₀F)
    simp only [hg]; split_ifs <;> linarith
  have hM₂1 : 1 ≤ M₂ := by linarith
  -- the extremal zeros
  set S : Finset UpperZeros := (Fs.filter (fun ρ => wn ρ = M)).subtype (fun ρ => ρ ∈ UpperZeros) with hS
  have hSmem : ∀ ρ : UpperZeros, ρ ∈ S ↔ (ρ : ℂ) ∈ Fs ∧ wn ρ = M := fun ρ => by
    rw [hS, Finset.mem_subtype, Finset.mem_filter]
  have hρmU : ρm ∈ UpperZeros := ⟨((hmemF ρm).1 hρmF).1, ((hmemF ρm).1 hρmF).2.1⟩
  have hSne : S.Nonempty := ⟨⟨ρm, hρmU⟩, (hSmem _).2 ⟨hρmF, hρmM.symm⟩⟩
  have hout : ∀ ρ : UpperZeros, ρ ∉ S → wn ρ ≤ M₂ := by
    intro ρ hρ
    by_cases hF : (ρ : ℂ) ∈ Fs
    · have h1 : wn ρ ≤ M := Finset.le_sup' wn hF
      have h2 : wn ρ ≠ M := fun h => hρ ((hSmem ρ).2 ⟨hF, h⟩)
      have h3 : wn ρ < M := lt_of_le_of_ne h1 h2
      have h4 : g ρ = wn ρ := by simp only [hg]; rw [if_pos h3]
      rw [← h4]; exact Finset.le_sup' g hF
    · have h5 : wn ρ < τ := by
        by_contra h; push Not at h
        exact hF ((hmemF ρ).2 ⟨ρ.2.1, ρ.2.2, h⟩)
      linarith
  -- split λ_n
  have hsplit : ∀ n : ℕ, lam n = (∑ ρ ∈ S, upperTerm n ρ) +
      ∑' ρ : ↑((S : Set UpperZeros)ᶜ), upperTerm n ρ := fun n => by
    rw [lam_eq L n, (summable_upperTerm L n).sum_add_tsum_compl]
  set A := ∑ ρ ∈ S, 2 * (mult (ρ : ℂ) : ℝ) with hA
  set φ : UpperZeros → ℝ := fun ρ => Complex.arg (1 - 1 / (ρ : ℂ)) with hφ
  have hhead : ∀ n : ℕ, ∑ ρ ∈ S, upperTerm n ρ =
      A - M ^ n * ∑ ρ ∈ S, 2 * (mult (ρ : ℂ) : ℝ) * Real.cos (n * φ ρ) := by
    intro n
    rw [hA, Finset.mul_sum, ← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl fun ρ hρ => ?_
    have hwM : wn ρ = M := ((hSmem ρ).1 hρ).2
    have hpolar : (1 - 1 / (ρ : ℂ)) = (M : ℂ) * Complex.exp ((φ ρ : ℂ) * I) := by
      rw [← hwM]; exact (Complex.norm_mul_exp_arg_mul_I _).symm
    have hre : ((1 - 1 / (ρ : ℂ)) ^ n).re = M ^ n * Real.cos (n * φ ρ) := by
      rw [hpolar, mul_pow, ← Complex.exp_nat_mul,
        show (n : ℂ) * ((φ ρ : ℂ) * I) = ((n * φ ρ : ℝ) : ℂ) * I by push_cast; ring,
        ← Complex.ofReal_pow, Complex.re_ofReal_mul, Complex.exp_ofReal_mul_I_re]
    unfold upperTerm
    rw [Complex.sub_re, Complex.one_re, hre]
    ring
  -- tail
  set Z := ∑' ρ : UpperZeros, (mult (ρ : ℂ) : ℝ) / (ρ : ℂ).im ^ 2 with hZ
  have hZs := summable_mult_div_sq L
  have hZ0 : 0 ≤ Z := tsum_nonneg fun ρ => by positivity
  have htail : ∀ n : ℕ, |∑' ρ : ↑((S : Set UpperZeros)ᶜ), upperTerm n ρ| ≤ 2 * (n + n ^ 2 * M₂ ^ n) * Z := by
    intro n
    set b : UpperZeros → ℝ := fun ρ => 2 * (n + n ^ 2 * M₂ ^ n) * ((mult (ρ : ℂ) : ℝ) / (ρ : ℂ).im ^ 2) with hb
    have hbs : Summable b := hZs.mul_left _
    have hb0 : ∀ ρ, 0 ≤ b ρ := fun ρ => by simp only [hb]; positivity
    have hpt : ∀ ρ : ↑((S : Set UpperZeros)ᶜ), |upperTerm n ρ| ≤ b ρ :=
      fun ρ => abs_upperTerm_le_K hM₂1 (ρ : UpperZeros).2.1 (ρ : UpperZeros).2.2 (hout ρ ρ.2)
    have hus : Summable fun ρ : ↑((S : Set UpperZeros)ᶜ) => upperTerm n ρ := (summable_upperTerm L n).subtype _
    have hbs' : Summable fun ρ : ↑((S : Set UpperZeros)ᶜ) => b ρ := hbs.subtype _
    have h1 : ∑' ρ : ↑((S : Set UpperZeros)ᶜ), upperTerm n ρ ≤ ∑' ρ : ↑((S : Set UpperZeros)ᶜ), b ρ :=
      hus.tsum_le_tsum (fun ρ => (le_abs_self _).trans (hpt ρ)) hbs'
    have h2 : -∑' ρ : ↑((S : Set UpperZeros)ᶜ), b ρ ≤ ∑' ρ : ↑((S : Set UpperZeros)ᶜ), upperTerm n ρ := by
      rw [← tsum_neg]
      exact hbs'.neg.tsum_le_tsum (fun ρ => by have := hpt ρ; rw [abs_le] at this; linarith) hus
    have h3 : ∑' ρ : ↑((S : Set UpperZeros)ᶜ), b ρ ≤ ∑' ρ : UpperZeros, b ρ :=
      Summable.tsum_subtype_le b _ hb0 hbs
    have h4 : ∑' ρ : UpperZeros, b ρ = 2 * (n + n ^ 2 * M₂ ^ n) * Z := by
      simp only [hb, hZ]; exact tsum_mul_left
    rw [abs_le]; constructor <;> linarith
  -- oscillation
  obtain ⟨c, hc0, hfreq⟩ := osc S hSne (fun ρ => 2 * (mult (ρ : ℂ) : ℝ)) φ
    (fun ρ _ => by have := mult_pos (ρ : UpperZeros).2.1; positivity)
  -- eventual domination
  have hr1 : |1 / M| < 1 := by rw [abs_of_pos (by positivity)]; exact (div_lt_one hM0).2 hM1
  have hr2 : |M₂ / M| < 1 := by
    rw [abs_of_pos (by have : 0 < M₂ := by linarith
                       positivity)]
    exact (div_lt_one hM0).2 hM₂
  have ht0 : Tendsto (fun n : ℕ => (1 / M) ^ n) atTop (𝓝 0) :=
    tendsto_pow_atTop_nhds_zero_of_lt_one (by positivity) (by rw [abs_of_pos (by positivity)] at hr1; exact hr1)
  have ht1 : Tendsto (fun n : ℕ => (n : ℝ) ^ 1 * (1 / M) ^ n) atTop (𝓝 0) :=
    tendsto_pow_const_mul_const_pow_of_abs_lt_one 1 hr1
  have ht2 : Tendsto (fun n : ℕ => (n : ℝ) ^ 2 * (M₂ / M) ^ n) atTop (𝓝 0) :=
    tendsto_pow_const_mul_const_pow_of_abs_lt_one 2 hr2
  have hu : Tendsto (fun n : ℕ => A * (1 / M) ^ n + 2 * Z * ((n : ℝ) ^ 1 * (1 / M) ^ n) +
      2 * Z * ((n : ℝ) ^ 2 * (M₂ / M) ^ n)) atTop (𝓝 0) := by
    have := ((ht0.const_mul A).add (ht1.const_mul (2 * Z))).add (ht2.const_mul (2 * Z))
    simpa using this
  have hev1 : ∀ᶠ n : ℕ in atTop, A + 2 * (n + n ^ 2 * M₂ ^ n) * Z ≤ M ^ n * c / 2 := by
    filter_upwards [hu.eventually (Iio_mem_nhds (show (0 : ℝ) < c / 2 by positivity))] with n hn
    have hMn : 0 < M ^ n := pow_pos hM0 n
    have key : A + 2 * (n + n ^ 2 * M₂ ^ n) * Z = M ^ n * (A * (1 / M) ^ n + 2 * Z * ((n : ℝ) ^ 1 * (1 / M) ^ n) +
        2 * Z * ((n : ℝ) ^ 2 * (M₂ / M) ^ n)) := by
      have hMne : M ≠ 0 := hM0.ne'
      rw [div_pow, div_pow, one_pow]
      field_simp
      ring
    rw [key]
    have := mul_lt_mul_of_pos_left hn hMn
    linarith
  set lM := Real.log M with hlM
  have hlM0 : 0 < lM := Real.log_pos hM1
  have hev2 : ∀ᶠ n : ℕ in atTop, 2 / c ≤ Real.exp (n * (lM / 2)) :=
    (Real.tendsto_exp_atTop.comp (tendsto_natCast_atTop_atTop.atTop_mul_const (by positivity))).eventually_ge_atTop _
  refine ⟨lM / 2, by positivity, (hfreq.and_eventually (hev1.and hev2)).mono fun n ⟨h1, h2, h3⟩ => ?_⟩
  rw [hsplit n, hhead n]
  have ht := htail n
  set f := ∑ ρ ∈ S, 2 * (mult (ρ : ℂ) : ℝ) * Real.cos (n * φ ρ) with hf
  set T := ∑' ρ : ↑((S : Set UpperZeros)ᶜ), upperTerm n ρ with hT
  have hMn : 0 < M ^ n := pow_pos hM0 n
  have hMf : M ^ n * c ≤ |M ^ n * f| := by
    rw [abs_mul, abs_of_pos hMn]; exact mul_le_mul_of_nonneg_left h1 hMn.le
  have hlow : M ^ n * c / 2 ≤ |A - M ^ n * f + T| := by
    have e1 := abs_sub_abs_le_abs_sub (M ^ n * f) (A + T)
    have e2 : |A + T| ≤ A + 2 * (n + n ^ 2 * M₂ ^ n) * Z := by
      have hA0 : 0 ≤ A := Finset.sum_nonneg fun ρ _ => by positivity
      calc |A + T| ≤ |A| + |T| := abs_add_le _ _
        _ = A + |T| := by rw [abs_of_nonneg hA0]
        _ ≤ A + 2 * (n + n ^ 2 * M₂ ^ n) * Z := by linarith
    have e3 : |M ^ n * f - (A + T)| = |A - M ^ n * f + T| := by
      rw [← abs_neg]; ring_nf
    linarith
  have hMexp : M ^ n = Real.exp (n * (lM / 2)) * Real.exp (n * (lM / 2)) := by
    rw [← Real.exp_add, show (n : ℝ) * (lM / 2) + n * (lM / 2) = n * lM by ring, Real.exp_nat_mul,
      hlM, Real.exp_log hM0]
  have hfin : Real.exp (lM / 2 * n) ≤ M ^ n * c / 2 := by
    rw [hMexp, mul_comm (lM / 2)]
    have he := Real.exp_pos (n * (lM / 2))
    have : 2 / c * (c / 2) = 1 := by field_simp
    calc Real.exp (n * (lM / 2)) = Real.exp (n * (lM / 2)) * (2 / c * (c / 2)) := by rw [this, mul_one]
      _ ≤ Real.exp (n * (lM / 2)) * (Real.exp (n * (lM / 2)) * (c / 2)) := by
          apply mul_le_mul_of_nonneg_left _ he.le
          exact mul_le_mul_of_nonneg_right h3 (by positivity)
      _ = Real.exp (n * (lM / 2)) * Real.exp (n * (lM / 2)) * c / 2 := by ring
  linarith

end LeanLi
