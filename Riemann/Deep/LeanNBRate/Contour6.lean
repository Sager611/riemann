import Riemann.Deep.LeanNBRate.Contour5

/-!
# U1 contour argument, part 6: residues at the trivial zeros and the telescoped left half-plane
-/

open Complex Filter Topology Set MeasureTheory intervalIntegral
open scoped Real Interval Nat

noncomputable section

namespace Riemann.Deep.LeanNBRate

variable {H : ℂ → ℂ} {M β D : ℝ} {N : ℕ}

/-- The residue of Psi_s at the trivial zero tz j. -/
def Ares (H : ℂ → ℂ) (L : ℝ) (s : ℂ) (j : ℕ) : ℂ := rk H j * kap L (s - tz j) / (s - tz j)

theorem s_sub_tz (s : ℂ) (j : ℕ) : s - tz j = s + 2 * ((j : ℂ) + 1) := by rw [tz]; ring

theorem Rtr_eq (H : ℂ → ℂ) (L : ℝ) (s : ℂ) : Rtr H L s = ∑' j : ℕ, Ares H L s j := by
  unfold Rtr Ares
  congr 1; funext j; rw [s_sub_tz]

theorem tz_ne_one (j : ℕ) : tz j ≠ 1 := fun h => by
  have := congrArg Complex.re h; rw [tz_re] at this; simp at this; linarith [(Nat.cast_nonneg j : (0:ℝ) ≤ j)]

theorem zeta_tz (j : ℕ) : riemannZeta (tz j) = 0 := by
  rw [tz]; exact riemannZeta_neg_two_mul_nat_add_one j

theorem deriv_zeta_tz_ne (j : ℕ) : deriv riemannZeta (tz j) ≠ 0 := by
  have h := norm_deriv_zeta_neg_two_ge j
  have hpos : (0 : ℝ) < ((2 * j + 2).factorial : ℝ) / (6 * (2 * π) ^ (2 * j + 2)) := by
    have : (0 : ℝ) < ((2 * j + 2).factorial : ℝ) := by exact_mod_cast Nat.factorial_pos _
    positivity
  intro h0
  rw [tz] at h0
  rw [h0, norm_zero] at h
  linarith

/-- Local expansion of Psi_s at the trivial zero tz j. -/
theorem Psi_pole_tz (hH : HBound H M β D N) (L : ℝ) {s : ℂ} (hs : 3 < s.re) (j : ℕ) :
    ∃ g : ℂ → ℂ, ContinuousAt g (tz j) ∧
      ∀ᶠ z in 𝓝[≠] (tz j), Psi H L s z = Ares H L s j / (z - tz j) + g z := by
  set p := tz j with hp
  set Q : ℂ → ℂ := fun z => H z * kap L (s - z) / (s - z) with hQ
  set ζ₁ : ℂ → ℂ := dslope riemannZeta p with hζ₁
  have hps : s - p ≠ 0 := by
    intro h; have := congrArg Complex.re h; rw [hp, sub_re, tz_re] at this; simp at this
    linarith [(Nat.cast_nonneg j : (0:ℝ) ≤ j)]
  have hU : {z : ℂ | z ≠ 1} ∈ 𝓝 p := isOpen_ne.mem_nhds (tz_ne_one j)
  have hζd : DifferentiableOn ℂ riemannZeta {z : ℂ | z ≠ 1} := fun z hz =>
    (differentiableAt_riemannZeta hz).differentiableWithinAt
  have hζ₁d : DifferentiableAt ℂ ζ₁ p :=
    ((differentiableOn_dslope hU).mpr hζd).differentiableAt hU
  have hζ₁p : ζ₁ p = deriv riemannZeta p := dslope_same _ _
  have hζ₁ne : ζ₁ p ≠ 0 := by rw [hζ₁p]; exact deriv_zeta_tz_ne j
  have hQd : DifferentiableAt ℂ Q p := by
    simp only [hQ, kap]
    exact ((hH.diff p).mul (by fun_prop)).div ((differentiableAt_const _).sub differentiableAt_id) hps
  have hRd : DifferentiableAt ℂ (fun z => Q z / ζ₁ z) p := hQd.div hζ₁d hζ₁ne
  refine ⟨dslope (fun z => Q z / ζ₁ z) p, continuousAt_dslope_same.mpr hRd, ?_⟩
  have hev : ∀ᶠ z in 𝓝 p, ζ₁ z ≠ 0 := hζ₁d.continuousAt.eventually_ne hζ₁ne
  filter_upwards [nhdsWithin_le_nhds hev, self_mem_nhdsWithin] with z hz1 hz2
  have hzp : z - p ≠ 0 := sub_ne_zero.mpr hz2
  have hζz : riemannZeta z = (z - p) * ζ₁ z := by
    have := sub_smul_dslope riemannZeta p z
    rw [smul_eq_mul, zeta_tz j, sub_zero] at this
    exact this.symm
  rw [dslope_of_ne _ hz2, slope_def_field]
  have hrk : rk H j = H p / deriv riemannZeta p := rfl
  have hd0 : deriv riemannZeta p ≠ 0 := deriv_zeta_tz_ne j
  have hA : Ares H L s j = Q p / ζ₁ p := by
    rw [Ares, hrk, hζ₁p, ← hp]
    simp only [hQ]
    first | (field_simp; done) | (field_simp; ring)
  rw [hA]
  simp only [Psi, Fq, hζz, hQ]
  field_simp
  ring

/-- No trivial zero has real part in [hl (j+1), hl j] except tz j. -/
theorem tz_in_strip {j k : ℕ} (h1 : hl (j + 1) ≤ (tz k).re) (h2 : (tz k).re ≤ hl j) : k = j := by
  rw [tz_re] at h1 h2
  unfold hl at h1 h2
  push_cast at h1 h2
  have a1 : (k : ℝ) < j + 1 := by linarith
  have a2 : (j : ℝ) < k + 1 := by linarith
  have b1 : k < j + 1 := by exact_mod_cast a1
  have b2 : j < k + 1 := by exact_mod_cast a2
  omega

/-- Step F: one strip, one residue. -/
theorem VI_Psi_hl_step (hH : HBound H M β D N) {L : ℝ} (hL : 0 ≤ L) {s : ℂ} (hs : 3 < s.re) (j : ℕ) :
    I * VI (Psi H L s) (hl j) - I * VI (Psi H L s) (hl (j + 1)) = 2 * π * I * Ares H L s j := by
  have hjle : hl j ≤ -1 := by unfold hl; linarith [(Nat.cast_nonneg j : (0:ℝ) ≤ j)]
  have hab : hl (j + 1) < hl j := by unfold hl; push_cast; linarith
  have hnot : ∀ x : ℂ, hl (j + 1) ≤ x.re → x.re ≤ hl j → x ∉ ({tz j} : Finset ℂ) → ∀ k : ℕ, x ≠ tz k := by
    intro x h1 h2 hx k hk
    subst hk
    have := tz_in_strip h1 h2
    subst this
    exact hx (Finset.mem_singleton_self _)
  have h := strip_residues (f := Psi H L s) (a := hl (j + 1)) (b := hl j) hab {tz j} (fun _ => Ares H L s j)
    (by
      intro p hp; rw [Finset.mem_singleton] at hp; subst hp
      rw [tz_re]; unfold hl; push_cast; constructor <;> linarith)
    (fun x h1 h2 hx => continuousAt_Psi_left hH L hs (h2.trans hjle) (hnot x h1 h2 hx))
    (fun x h1 h2 hx => differentiableAt_Psi_left hH L hs (h2.le.trans hjle) (hnot x h1.le h2.le hx))
    (by
      intro p hp; rw [Finset.mem_singleton] at hp; subst hp
      exact Psi_pole_tz hH L hs j)
    (integrable_Psi_hl hH hL hs (j + 1)) (integrable_Psi_hl hH hL hs j)
    (decay_of_le_inv_sq fun z h1 h2 h3 => norm_Psi_left_strip hH hL hs hab.le hjle h1 h2 h3)
  simpa using h

theorem neg_one_ne_tz (y : ℝ) (k : ℕ) : ((-1 : ℝ) : ℂ) + y * I ≠ tz k :=
  line_ne_tz (fun k h => by linarith [(Nat.cast_nonneg k : (0:ℝ) ≤ k)]) y k

theorem integrable_Psi_neg_one (hH : HBound H M β D N) {L : ℝ} (hL : 0 ≤ L) {s : ℂ} (hs : 3 < s.re) :
    Integrable fun y : ℝ => Psi H L s ((-1 : ℝ) + y * I) := by
  have hc : Continuous fun y : ℝ => Psi H L s ((-1 : ℝ) + y * I) := by
    rw [continuous_iff_continuousAt]
    intro y
    exact ContinuousAt.comp (f := fun y : ℝ => ((-1 : ℝ) : ℂ) + y * I) (x := y)
      (continuousAt_Psi_left hH L hs (by simp) (neg_one_ne_tz y))
      (by fun_prop : Continuous fun y : ℝ => ((-1 : ℝ) : ℂ) + y * I).continuousAt
  refine integrable_of_cont_le_inv_sq hc (Y₀ := 1) (C := M * Real.exp (β * (1 / 2 - -1)) *
    (6 * Real.sqrt (2 * π) * (2 * π) ^ (-(-1 : ℝ)) * (|(-1 : ℝ)| + (3 / 2 - -1) + 1 + 1)) * (D / min D 1) ^ N)
    fun y hy => ?_
  have := norm_Psi_left_strip hH hL hs (a := -1) (b := -1) le_rfl le_rfl (z := ((-1 : ℝ) : ℂ) + y * I)
    (by simp) (by simp) (by simpa using hy)
  simpa using this

/-- Step E: the strip [-3/2, -1]. -/
theorem VI_Psi_neg_one (hH : HBound H M β D N) {L : ℝ} (hL : 0 ≤ L) {s : ℂ} (hs : 3 < s.re) :
    I * VI (Psi H L s) (-1) - I * VI (Psi H L s) (hl 0) = 0 := by
  have h0 : hl 0 = -3 / 2 := by unfold hl; simp; norm_num
  have hnot : ∀ x : ℂ, hl 0 ≤ x.re → x.re ≤ -1 → ∀ k : ℕ, x ≠ tz k := by
    intro x h1 _ k hk; rw [hk, tz_re, h0] at h1; linarith [(Nat.cast_nonneg k : (0:ℝ) ≤ k)]
  have h := strip_residues (f := Psi H L s) (a := hl 0) (b := -1) (by rw [h0]; norm_num) ∅ (fun _ => 0)
    (by simp) (fun x h1 h2 _ => continuousAt_Psi_left hH L hs h2 (hnot x h1 h2))
    (fun x h1 h2 _ => differentiableAt_Psi_left hH L hs h2.le (hnot x h1.le h2.le)) (by simp)
    (integrable_Psi_hl hH hL hs 0) (integrable_Psi_neg_one hH hL hs)
    (decay_of_le_inv_sq fun z h1 h2 h3 => norm_Psi_left_strip hH hL hs (by rw [h0]; norm_num) le_rfl h1 h2 h3)
  simpa using h

/-- Telescoping: I VI(hl 0) - I VI(hl K) = 2 pi i sum_{j<K} A_j. -/
theorem VI_Psi_telescope (hH : HBound H M β D N) {L : ℝ} (hL : 0 ≤ L) {s : ℂ} (hs : 3 < s.re) (K : ℕ) :
    I * VI (Psi H L s) (hl 0) - I * VI (Psi H L s) (hl K) = 2 * π * I * ∑ j ∈ Finset.range K, Ares H L s j := by
  induction K with
  | zero => simp
  | succ K ih =>
    rw [Finset.sum_range_succ, mul_add, ← ih, ← VI_Psi_hl_step hH hL hs K]
    ring

end Riemann.Deep.LeanNBRate

