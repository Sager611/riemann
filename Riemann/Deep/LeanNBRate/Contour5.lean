import Riemann.Deep.LeanNBRate.Contour4

/-!
# U1 contour argument, part 5: the left strips (trivial-zero residues) and the far-left limit
-/

open Complex Filter Topology Set MeasureTheory intervalIntegral
open scoped Real Interval

noncomputable section

namespace Riemann.Deep.LeanNBRate

variable {H : ℂ → ℂ} {M β D : ℝ} {N : ℕ}

/-- The trivial zero -2(j+1). -/
def tz (j : ℕ) : ℂ := -2 * ((j : ℂ) + 1)

/-- The half-integer line -2j - 3/2. -/
def hl (j : ℕ) : ℝ := -2 * j - 3 / 2

theorem zeta_ne_zero_left {z : ℂ} (hz : z.re < 0) (hk : ∀ k : ℕ, z ≠ tz k) : riemannZeta z ≠ 0 := by
  intro h
  obtain ⟨k, hk'⟩ := zeta_eq_zero_left h hz
  exact hk k (by rw [hk', tz])

theorem continuousAt_Psi_left (hH : HBound H M β D N) (L : ℝ) {s z : ℂ} (hs : 3 < s.re)
    (hz : z.re ≤ -1) (hk : ∀ k : ℕ, z ≠ tz k) : ContinuousAt (Psi H L s) z :=
  continuousAt_Psi hH L (fun h => by rw [h] at hz; norm_num at hz) (zeta_ne_zero_left (by linarith) hk)
    (fun h => by rw [h] at hz; linarith)

theorem differentiableAt_Psi_left (hH : HBound H M β D N) (L : ℝ) {s z : ℂ} (hs : 3 < s.re)
    (hz : z.re ≤ -1) (hk : ∀ k : ℕ, z ≠ tz k) : DifferentiableAt ℂ (Psi H L s) z :=
  differentiableAt_Psi hH L (fun h => by rw [h] at hz; norm_num at hz) (zeta_ne_zero_left (by linarith) hk)
    (fun h => by rw [h] at hz; linarith)

theorem tz_re (j : ℕ) : (tz j).re = -2 * j - 2 := by simp [tz]; ring

theorem line_ne_tz {σ : ℝ} (hσ : ∀ k : ℕ, σ ≠ -2 * k - 2) (y : ℝ) (k : ℕ) : (σ : ℂ) + y * I ≠ tz k := by
  intro h
  have := congrArg Complex.re h
  rw [tz_re] at this
  simp at this
  exact hσ k (by linarith)

theorem hl_ne (j k : ℕ) : hl j ≠ -2 * k - 2 := by
  intro h
  unfold hl at h
  have h2 : (4 * j + 3 : ℝ) = 4 * k + 4 := by linarith
  have h3 : (4 * j + 3 : ℕ) = 4 * k + 4 := by exact_mod_cast h2
  omega

/-- Product bound on the half-integer lines. -/
theorem half_prod_le (j : ℕ) (y : ℝ) :
    1 / 2 ≤ ∏ i ∈ Finset.range (2 * j + 2), ‖((hl j : ℝ) : ℂ) + y * I + i‖ := by
  have hre : ∀ i : ℕ, (((hl j : ℝ) : ℂ) + y * I + i).re = hl j + i := fun i => by simp
  have hfac : ∀ i : ℕ, i ≤ 2 * j → 1 ≤ ‖((hl j : ℝ) : ℂ) + y * I + i‖ := by
    intro i hi
    have h1 := Complex.abs_re_le_norm (((hl j : ℝ) : ℂ) + y * I + i)
    rw [hre] at h1
    have hi' : (i : ℝ) ≤ 2 * j := by exact_mod_cast hi
    have hneg : hl j + i < 0 := by unfold hl; linarith
    rw [abs_of_neg hneg] at h1
    have : 1 ≤ -(hl j + i) := by unfold hl; linarith
    linarith
  have hlast : 1 / 2 ≤ ‖((hl j : ℝ) : ℂ) + y * I + ((2 * j + 1 : ℕ) : ℂ)‖ := by
    have h1 := Complex.abs_re_le_norm (((hl j : ℝ) : ℂ) + y * I + ((2 * j + 1 : ℕ) : ℂ))
    rw [hre] at h1
    have e : hl j + ((2 * j + 1 : ℕ) : ℝ) = -1 / 2 := by unfold hl; push_cast; ring
    rw [e] at h1
    norm_num at h1
    push_cast at h1 ⊢
    linarith
  have hfirst : 1 ≤ ∏ i ∈ Finset.range (2 * j + 1), ‖((hl j : ℝ) : ℂ) + y * I + i‖ := by
    have : ∀ n : ℕ, n ≤ 2 * j + 1 → 1 ≤ ∏ i ∈ Finset.range n, ‖((hl j : ℝ) : ℂ) + y * I + i‖ := by
      intro n
      induction n with
      | zero => intro _; simp
      | succ n ih =>
        intro hn
        rw [Finset.prod_range_succ]
        have h1 := ih (by omega)
        have h2 := hfac n (by omega)
        nlinarith
    exact this _ le_rfl
  rw [show 2 * j + 2 = (2 * j + 1) + 1 by ring, Finset.prod_range_succ]
  nlinarith

/-- Bound on the half-integer line hl j, all y. -/
theorem norm_Psi_hl (hH : HBound H M β D N) {L : ℝ} (hL : 0 ≤ L) {s : ℂ} (hs : 3 < s.re) (j : ℕ) (y : ℝ) :
    ‖Psi H L s ((hl j : ℝ) + y * I)‖ ≤ (M * Real.exp (β * (1 / 2 - hl j)) *
      (12 * Real.sqrt (2 * π) * (2 * π) ^ (-hl j) * (4 * j + 5)) * Real.exp (-(L * (s.re - hl j))) *
      (D / min D 1) ^ N) * ((1 + |y|) ^ 2)⁻¹ := by
  set z : ℂ := ((hl j : ℝ) : ℂ) + y * I with hz
  have hzre : z.re = hl j := by simp [hz]
  have hzim : z.im = y := by simp [hz]
  have hle : z.re ≤ -1 := by rw [hzre]; unfold hl; linarith [(Nat.cast_nonneg j : (0 : ℝ) ≤ j)]
  have hzn : ∀ n : ℕ, z ≠ -n := by
    intro n h
    have := congrArg Complex.re h
    rw [hzre] at this; simp at this
    unfold hl at this
    have h2 : (4 * j + 3 : ℝ) = 2 * n := by linarith
    have h3 : (4 * j + 3 : ℕ) = 2 * n := by exact_mod_cast h2
    omega
  have hm1 : 1 / 2 ≤ z.re + ((2 * j + 2 : ℕ) : ℝ) := by rw [hzre]; unfold hl; push_cast; linarith
  have hm2 : z.re + ((2 * j + 2 : ℕ) : ℝ) ≤ 3 / 2 := by rw [hzre]; unfold hl; push_cast; linarith
  have h := norm_Psi_left hH hL hs hle hzn (2 * j + 2) hm1 hm2 (by norm_num) (half_prod_le j y)
  rw [hzre, hzim] at h
  have hM := hH.M_nonneg
  have hD := hH.D_pos
  have hnz : ‖z‖ + ((2 * j + 2 : ℕ) : ℝ) + 1 ≤ (4 * j + 5) * (1 + |y|) := by
    have h1 : ‖z‖ ≤ |hl j| + |y| := by
      calc ‖z‖ ≤ |z.re| + |z.im| := Complex.norm_le_abs_re_add_abs_im z
        _ = |hl j| + |y| := by rw [hzre, hzim]
    have h2 : |hl j| = 2 * j + 3 / 2 := by
      unfold hl; rw [abs_of_neg (by linarith [(Nat.cast_nonneg j : (0 : ℝ) ≤ j)])]; ring
    push_cast
    nlinarith [abs_nonneg y, (Nat.cast_nonneg j : (0 : ℝ) ≤ j)]
  have hwl := weight_le hD (by linarith [hH.N_ge] : 4 ≤ N) y
  have hy1 : 1 + |y| ≤ (1 + |y|) ^ 2 := by nlinarith [abs_nonneg y]
  calc ‖Psi H L s z‖ ≤ M * Real.exp (β * (1 / 2 - hl j)) * (D / (D + |y|)) ^ N *
        (6 * Real.sqrt (2 * π) * (2 * π) ^ (-hl j) * (‖z‖ + ((2 * j + 2 : ℕ) : ℝ) + 1) / (1 / 2)) *
        Real.exp (-(L * (s.re - hl j))) := h
    _ ≤ M * Real.exp (β * (1 / 2 - hl j)) * (D / (D + |y|)) ^ N *
        (6 * Real.sqrt (2 * π) * (2 * π) ^ (-hl j) * ((4 * j + 5) * (1 + |y|) ^ 2) / (1 / 2)) *
        Real.exp (-(L * (s.re - hl j))) := by
        gcongr
        exact hnz.trans (mul_le_mul_of_nonneg_left hy1 (by positivity))
    _ = M * Real.exp (β * (1 / 2 - hl j)) * (12 * Real.sqrt (2 * π) * (2 * π) ^ (-hl j) * (4 * j + 5)) *
        Real.exp (-(L * (s.re - hl j))) * ((1 + |y|) ^ 2 * (D / (D + |y|)) ^ N) := by ring
    _ ≤ M * Real.exp (β * (1 / 2 - hl j)) * (12 * Real.sqrt (2 * π) * (2 * π) ^ (-hl j) * (4 * j + 5)) *
        Real.exp (-(L * (s.re - hl j))) * ((D / min D 1) ^ N * ((1 + |y|) ^ 2)⁻¹) := by gcongr
    _ = _ := by ring

theorem integrable_Psi_hl (hH : HBound H M β D N) {L : ℝ} (hL : 0 ≤ L) {s : ℂ} (hs : 3 < s.re) (j : ℕ) :
    Integrable fun y : ℝ => Psi H L s ((hl j : ℝ) + y * I) := by
  have hc : Continuous fun y : ℝ => Psi H L s ((hl j : ℝ) + y * I) := by
    rw [continuous_iff_continuousAt]
    intro y
    have hle : (((hl j : ℝ) : ℂ) + y * I).re ≤ -1 := by
      simp; unfold hl; linarith [(Nat.cast_nonneg j : (0 : ℝ) ≤ j)]
    exact ContinuousAt.comp (f := fun y : ℝ => ((hl j : ℝ) : ℂ) + y * I) (x := y)
      (continuousAt_Psi_left hH L hs hle (line_ne_tz (hl_ne j) y))
      (by fun_prop : Continuous fun y : ℝ => ((hl j : ℝ) : ℂ) + y * I).continuousAt
  exact integrable_of_le_inv_sq hc.aestronglyMeasurable (norm_Psi_hl hH hL hs j)

end Riemann.Deep.LeanNBRate

