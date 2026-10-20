import Mathlib

/-!
# Hardy's theorem, part 1: the theta kernel

`Hk u = e^{-u/4} θ(e^{-u}) - 2 cosh(u/4)` with `θ(y) = Σ_{n∈ℤ} e^{-π n² y}`.
On the real line `Λ(1/2 + i t) = (1/2) ∫ Hk(v) e^{-i v t / 2} dv` (see `Repr.lean`).
Here: holomorphy on the strip `|Im u| < π/2`, evenness, reality on `ℝ`, and the
uniform bound `‖Hk u‖ ≤ C e^{-|Re u|/4}` on closed sub-strips.
-/

open Complex Real MeasureTheory Set Filter Topology
open scoped ComplexConjugate

namespace Fayed.Hardy

/-- `θ(y) = Σ_{n∈ℤ} exp(-π n² y)` (meaningful for `Re y > 0`). -/
noncomputable def theta (y : ℂ) : ℂ := jacobiTheta₂ 0 (I * y)

/-- Theta part of the kernel. -/
noncomputable def Rk (u : ℂ) : ℂ := cexp (-u / 4) * theta (cexp (-u))

/-- Elementary part of the kernel. -/
noncomputable def Ck (u : ℂ) : ℂ := 2 * Complex.cosh (u / 4)

/-- The kernel. -/
noncomputable def Hk (u : ℂ) : ℂ := Rk u - Ck u

/-- The strip `|Im u| < π/2`. -/
def strip : Set ℂ := {u | |u.im| < π / 2}

lemma isOpen_strip : IsOpen strip :=
  isOpen_lt (continuous_abs.comp continuous_im) continuous_const

lemma re_exp_neg_pos {u : ℂ} (hu : u ∈ strip) : 0 < (cexp (-u)).re := by
  have h : |u.im| < π / 2 := hu
  rw [exp_re, neg_re, neg_im, Real.cos_neg]
  exact mul_pos (Real.exp_pos _) (Real.cos_pos_of_mem_Ioo
    ⟨by linarith [neg_abs_le u.im], by linarith [le_abs_self u.im]⟩)

lemma neg_mem_strip {u : ℂ} (hu : u ∈ strip) : -u ∈ strip := by
  show |(-u).im| < π / 2
  rw [neg_im, abs_neg]; exact hu

lemma differentiableAt_theta {y : ℂ} (hy : 0 < y.re) : DifferentiableAt ℂ theta y := by
  have h1 : 0 < (I * y).im := by simpa using hy
  exact (differentiableAt_jacobiTheta₂_snd 0 h1).comp y (differentiableAt_id.const_mul I)

lemma differentiableAt_Rk {u : ℂ} (hu : u ∈ strip) : DifferentiableAt ℂ Rk u := by
  have h1 : DifferentiableAt ℂ (fun u : ℂ => cexp (-u / 4)) u := by fun_prop
  have h2 : DifferentiableAt ℂ (fun u : ℂ => cexp (-u)) u := by fun_prop
  have h3 := DifferentiableAt.comp (g := theta) u (differentiableAt_theta (re_exp_neg_pos hu)) h2
  exact h1.mul h3

lemma differentiable_Ck : Differentiable ℂ Ck := by
  unfold Ck; fun_prop

lemma differentiableOn_Rk : DifferentiableOn ℂ Rk strip :=
  fun _ hu => (differentiableAt_Rk hu).differentiableWithinAt

lemma differentiableOn_Hk : DifferentiableOn ℂ Hk strip :=
  fun _ hu => ((differentiableAt_Rk hu).sub (differentiable_Ck _)).differentiableWithinAt

/-- `Hk` is real on the real line. -/
lemma conj_Hk_ofReal (x : ℝ) : conj (Hk x) = Hk x := by
  have he : conj (cexp (-(x : ℂ))) = cexp (-(x : ℂ)) := by
    rw [← exp_conj, map_neg, conj_ofReal]
  have he4 : conj (cexp (-(x : ℂ) / 4)) = cexp (-(x : ℂ) / 4) := by
    rw [← exp_conj, map_div₀, map_neg, conj_ofReal, map_ofNat]
  have hth : conj (theta (cexp (-(x : ℂ)))) = theta (cexp (-(x : ℂ))) := by
    simp only [theta, jacobiTheta₂_conj, map_zero, map_mul, conj_I, he]
    congr 1; ring
  have hc : conj (Ck x) = Ck x := by
    simp only [Ck, map_mul, map_ofNat, ← cosh_conj, map_div₀, conj_ofReal]
  simp only [Hk, Rk, map_sub, map_mul, he4, hth, hc]

/-- The theta part is even on the strip (theta functional equation). -/
lemma Rk_neg {u : ℂ} (hu : u ∈ strip) : Rk (-u) = Rk u := by
  have h : |u.im| < π / 2 := hu
  have hlog : Complex.log (cexp (-u)) = -u := by
    apply Complex.log_exp
    · rw [neg_im]; linarith [le_abs_self u.im, Real.pi_pos]
    · rw [neg_im]; linarith [neg_abs_le u.im, Real.pi_pos]
  have hFE := jacobiTheta₂_functional_equation 0 (I * cexp (-u))
  have h1 : -I * (I * cexp (-u)) = cexp (-u) := by
    rw [← mul_assoc, neg_mul, I_mul_I, neg_neg, one_mul]
  have h2 : -1 / (I * cexp (-u)) = I * cexp u := by
    rw [div_eq_iff (mul_ne_zero I_ne_zero (Complex.exp_ne_zero _))]
    rw [show I * cexp u * (I * cexp (-u)) = (I * I) * (cexp u * cexp (-u)) by ring,
      ← Complex.exp_add, add_neg_cancel, Complex.exp_zero, I_mul_I, mul_one]
  have h3 : cexp (-u) ^ (1 / 2 : ℂ) = cexp (-u / 2) := by
    rw [cpow_def_of_ne_zero (Complex.exp_ne_zero _), hlog]; ring_nf
  simp only [Rk, theta]
  rw [hFE, h1, h2, h3, neg_neg]
  simp only [ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true, zero_pow, mul_zero, zero_div,
    Complex.exp_zero, mul_one]
  rw [← mul_assoc, one_div, ← Complex.exp_neg, ← Complex.exp_add]
  congr 2; ring

lemma Ck_neg (u : ℂ) : Ck (-u) = Ck u := by
  simp only [Ck, neg_div, Complex.cosh_neg]

lemma Hk_neg {u : ℂ} (hu : u ∈ strip) : Hk (-u) = Hk u := by
  simp only [Hk, Rk_neg hu, Ck_neg]

lemma Hk_eq_aux (u : ℂ) :
    Hk u = cexp (-u / 4) * (jacobiTheta (I * cexp (-u)) - 1) - cexp (u / 4) := by
  simp only [Hk, Rk, Ck, theta, jacobiTheta_eq_jacobiTheta₂, Complex.two_cosh, neg_div]
  ring

/-- Uniform exponential decay of `Hk` on closed sub-strips, half-plane `Re u ≤ 0`. -/
lemma Hk_bound_left {b : ℝ} (hb0 : 0 ≤ b) (hb : b < π / 2) :
    ∃ C, ∀ u : ℂ, u.re ≤ 0 → |u.im| ≤ b → ‖Hk u‖ ≤ C * rexp (-|u.re| / 4) := by
  set c := Real.cos b with hc_def
  have hc : 0 < c := Real.cos_pos_of_mem_Ioo ⟨by linarith [Real.pi_pos], hb⟩
  have hpc : 0 < π * c := mul_pos Real.pi_pos hc
  have hE : rexp (-π * c) < 1 := by
    have := Real.exp_lt_exp.mpr (show -π * c < 0 by nlinarith)
    simpa using this
  set A := 2 / (1 - rexp (-π * c)) with hA_def
  have hA : 0 < A := div_pos two_pos (by linarith)
  refine ⟨A / (π * c) + 1, fun u hre him => ?_⟩
  set x := |u.re| with hx_def
  have hx : u.re = -x := by rw [hx_def, abs_of_nonpos hre]; ring
  have hx0 : 0 ≤ x := abs_nonneg _
  set s := rexp x with hs_def
  have hs1 : 1 ≤ s := Real.one_le_exp hx0
  -- imaginary part of τ
  have hcos : c ≤ Real.cos u.im := by
    rw [← Real.cos_abs u.im]
    exact Real.cos_le_cos_of_nonneg_of_le_pi (abs_nonneg _) (by linarith [Real.pi_pos]) him
  have hτ : s * c ≤ (I * cexp (-u)).im := by
    rw [I_mul_im, exp_re, neg_re, neg_im, Real.cos_neg, hx, neg_neg]
    exact mul_le_mul_of_nonneg_left hcos (Real.exp_pos _).le
  have hτpos : 0 < (I * cexp (-u)).im := lt_of_lt_of_le (mul_pos (by linarith) hc) hτ
  set Y := (I * cexp (-u)).im
  have hθ := norm_jacobiTheta_sub_one_le hτpos
  have hY1 : c ≤ Y := le_trans (by nlinarith) hτ
  have hEY : rexp (-π * Y) ≤ rexp (-π * c) := Real.exp_le_exp.mpr (by nlinarith [Real.pi_pos])
  have hEY' : rexp (-π * Y) ≤ rexp (-(π * c * s)) :=
    Real.exp_le_exp.mpr (by nlinarith [Real.pi_pos])
  have hA' : 2 / (1 - rexp (-π * Y)) ≤ A :=
    div_le_div_of_nonneg_left (by norm_num) (by linarith) (by linarith)
  have hθ' : ‖jacobiTheta (I * cexp (-u)) - 1‖ ≤ A * rexp (-(π * c * s)) :=
    hθ.trans (mul_le_mul hA' hEY' (Real.exp_pos _).le hA.le)
  -- key elementary inequality: e^{x/2} e^{-π c s} ≤ 1/(π c)
  have hkey : rexp (x / 2) * rexp (-(π * c * s)) ≤ 1 / (π * c) := by
    have h1 : rexp (x / 2) ≤ s := Real.exp_le_exp.mpr (by linarith)
    have h2 : π * c * s ≤ rexp (π * c * s) := by
      linarith [Real.add_one_le_exp (π * c * s)]
    have h3 : s * rexp (-(π * c * s)) ≤ 1 / (π * c) := by
      rw [Real.exp_neg, le_div_iff₀ hpc]
      have hp : 0 < rexp (π * c * s) := Real.exp_pos _
      rw [show s * (rexp (π * c * s))⁻¹ * (π * c) = (π * c * s) / rexp (π * c * s) by
        field_simp]
      rw [div_le_one hp]; exact h2
    exact (mul_le_mul_of_nonneg_right h1 (Real.exp_pos _).le).trans h3
  have hn1 : ‖cexp (-u / 4)‖ = rexp (x / 4) := by
    rw [Complex.norm_exp]; congr 1; simp [hx]
  have hn2 : ‖cexp (u / 4)‖ = rexp (-x / 4) := by
    rw [Complex.norm_exp]; congr 1; simp [hx]
  rw [Hk_eq_aux]
  calc ‖cexp (-u / 4) * (jacobiTheta (I * cexp (-u)) - 1) - cexp (u / 4)‖
      ≤ ‖cexp (-u / 4)‖ * ‖jacobiTheta (I * cexp (-u)) - 1‖ + ‖cexp (u / 4)‖ := by
        refine (norm_sub_le _ _).trans ?_; rw [norm_mul]
    _ ≤ rexp (x / 4) * (A * rexp (-(π * c * s))) + rexp (-x / 4) := by
        rw [hn1, hn2]; gcongr
    _ = A * (rexp (-x / 4) * (rexp (x / 2) * rexp (-(π * c * s)))) + rexp (-x / 4) := by
        rw [show rexp (x / 4) = rexp (-x / 4) * rexp (x / 2) by
          rw [← Real.exp_add]; congr 1; ring]
        ring
    _ ≤ A * (rexp (-x / 4) * (1 / (π * c))) + rexp (-x / 4) := by gcongr
    _ = (A / (π * c) + 1) * rexp (-x / 4) := by ring

/-- Uniform exponential decay of `Hk` on closed sub-strips. -/
lemma Hk_bound {b : ℝ} (hb0 : 0 ≤ b) (hb : b < π / 2) :
    ∃ C, ∀ u : ℂ, |u.im| ≤ b → ‖Hk u‖ ≤ C * rexp (-|u.re| / 4) := by
  obtain ⟨C, hC⟩ := Hk_bound_left hb0 hb
  refine ⟨C, fun u him => ?_⟩
  rcases le_total u.re 0 with h | h
  · exact hC u h him
  · have hu : u ∈ strip := show |u.im| < π / 2 by linarith
    rw [← Hk_neg hu]
    have := hC (-u) (by simp [h]) (by simpa using him)
    simpa using this

end Fayed.Hardy
