/-
Riemann/Deep/Yoshida/Taylor.lean — alternating Taylor bounds for `cos` and `sin` on `x ≥ 0`.

`cosT n x = Σ_{k<n} (−1)^k x^{2k}/(2k)!`, `sinT n x = Σ_{k<n} (−1)^k x^{2k+1}/(2k+1)!`. For `x ≥ 0`:
`cosT (2m+2) x ≤ cos x ≤ cosT (2m+3) x` and `sinT (2m+2) x ≤ sin x ≤ sinT (2m+1) x`
(induction, comparing derivatives on `[0, ∞)`).
-/
import Mathlib

noncomputable section

open Finset

namespace Fayed.Yoshida.Taylor

/-- Taylor partial sum of `cos`. -/
def cosT (n : ℕ) (x : ℝ) : ℝ := ∑ k ∈ range n, (-1 : ℝ) ^ k * x ^ (2 * k) / ((2 * k).factorial : ℝ)

/-- Taylor partial sum of `sin`. -/
def sinT (n : ℕ) (x : ℝ) : ℝ :=
  ∑ k ∈ range n, (-1 : ℝ) ^ k * x ^ (2 * k + 1) / ((2 * k + 1).factorial : ℝ)

lemma hasDerivAt_sinT (n : ℕ) (x : ℝ) : HasDerivAt (sinT n) (cosT n x) x := by
  have h : ∀ k ∈ range n, HasDerivAt
      (fun y : ℝ => (-1 : ℝ) ^ k * y ^ (2 * k + 1) / ((2 * k + 1).factorial : ℝ))
      ((-1 : ℝ) ^ k * x ^ (2 * k) / ((2 * k).factorial : ℝ)) x := by
    intro k _
    have h1 := ((hasDerivAt_pow (2 * k + 1) x).const_mul ((-1 : ℝ) ^ k)).div_const
      ((2 * k + 1).factorial : ℝ)
    convert h1 using 1
    rw [Nat.factorial_succ, Nat.add_sub_cancel]
    push_cast
    have : ((2 * k).factorial : ℝ) ≠ 0 := by positivity
    field_simp
  have hs := HasDerivAt.sum h
  unfold sinT cosT
  convert hs using 1
  funext y
  simp

lemma hasDerivAt_cosT (n : ℕ) (x : ℝ) : HasDerivAt (cosT (n + 1)) (-sinT n x) x := by
  have e : cosT (n + 1) = fun y => 1 + ∑ k ∈ range n,
      (-1 : ℝ) ^ (k + 1) * y ^ (2 * (k + 1)) / ((2 * (k + 1)).factorial : ℝ) := by
    funext y
    unfold cosT
    rw [sum_range_succ']
    simp [add_comm]
  have h : ∀ k ∈ range n, HasDerivAt
      (fun y : ℝ => (-1 : ℝ) ^ (k + 1) * y ^ (2 * (k + 1)) / ((2 * (k + 1)).factorial : ℝ))
      (-((-1 : ℝ) ^ k * x ^ (2 * k + 1) / ((2 * k + 1).factorial : ℝ))) x := by
    intro k _
    have h1 := ((hasDerivAt_pow (2 * (k + 1)) x).const_mul ((-1 : ℝ) ^ (k + 1))).div_const
      ((2 * (k + 1)).factorial : ℝ)
    convert h1 using 1
    have hf : ((2 * (k + 1)).factorial : ℝ) = (2 * k + 2) * ((2 * k + 1).factorial : ℝ) := by
      rw [show 2 * (k + 1) = (2 * k + 1) + 1 by ring, Nat.factorial_succ]
      push_cast
      ring
    rw [hf, show 2 * (k + 1) - 1 = 2 * k + 1 by omega]
    push_cast
    have : ((2 * k + 1).factorial : ℝ) ≠ 0 := by positivity
    field_simp
    ring
  rw [e]
  have hs := (hasDerivAt_const x (1 : ℝ)).add (HasDerivAt.sum h)
  convert hs using 1
  · funext y
    simp
  · unfold sinT
    rw [zero_add, ← sum_neg_distrib]

lemma cosT_zero (n : ℕ) : cosT (n + 1) 0 = 1 := by
  unfold cosT
  rw [sum_range_succ']
  simp

lemma sinT_zero (n : ℕ) : sinT n 0 = 0 := by
  unfold sinT
  simp

/-- Comparison of functions with ordered derivatives on `[0, ∞)`. -/
lemma le_of_deriv_le {f g f' g' : ℝ → ℝ} (hf : ∀ x, HasDerivAt f (f' x) x)
    (hg : ∀ x, HasDerivAt g (g' x) x) (h0 : f 0 ≤ g 0) (h : ∀ x, 0 ≤ x → f' x ≤ g' x)
    {x : ℝ} (hx : 0 ≤ x) : f x ≤ g x := by
  have hm : MonotoneOn (fun y => g y - f y) (Set.Ici 0) := by
    apply monotoneOn_of_hasDerivWithinAt_nonneg (convex_Ici 0) (f' := fun y => g' y - f' y)
    · exact fun y _ => ((hg y).sub (hf y)).continuousAt.continuousWithinAt
    · intro y _
      exact ((hg y).sub (hf y)).hasDerivWithinAt
    · intro y hy
      rw [interior_Ici] at hy
      have := h y (le_of_lt hy)
      linarith
  have := hm (Set.mem_Ici.2 le_rfl) (Set.mem_Ici.2 hx) hx
  simp only at this
  linarith

/-- All four alternating bounds at level `m`. -/
theorem bounds (m : ℕ) : (∀ x : ℝ, 0 ≤ x → Real.sin x ≤ sinT (2 * m + 1) x) ∧
    (∀ x : ℝ, 0 ≤ x → cosT (2 * m + 2) x ≤ Real.cos x) ∧
    (∀ x : ℝ, 0 ≤ x → sinT (2 * m + 2) x ≤ Real.sin x) ∧
    (∀ x : ℝ, 0 ≤ x → Real.cos x ≤ cosT (2 * m + 3) x) := by
  induction m with
  | zero =>
    have hs : ∀ x : ℝ, 0 ≤ x → Real.sin x ≤ sinT 1 x := by
      intro x hx
      simp [sinT]
      exact Real.sin_le hx
    exact ⟨hs, by
      have hc : ∀ x : ℝ, 0 ≤ x → cosT 2 x ≤ Real.cos x := fun x hx =>
        le_of_deriv_le (fun y => hasDerivAt_cosT 1 y) (fun y => Real.hasDerivAt_cos y)
          (by rw [cosT_zero, Real.cos_zero]) (fun y hy => by linarith [hs y hy]) hx
      have hs2 : ∀ x : ℝ, 0 ≤ x → sinT 2 x ≤ Real.sin x := fun x hx =>
        le_of_deriv_le (fun y => hasDerivAt_sinT 2 y) (fun y => Real.hasDerivAt_sin y)
          (by rw [sinT_zero, Real.sin_zero]) hc hx
      have hc3 : ∀ x : ℝ, 0 ≤ x → Real.cos x ≤ cosT 3 x := fun x hx =>
        le_of_deriv_le (fun y => Real.hasDerivAt_cos y) (fun y => hasDerivAt_cosT 2 y)
          (by rw [cosT_zero, Real.cos_zero]) (fun y hy => by linarith [hs2 y hy]) hx
      exact ⟨hc, hs2, hc3⟩⟩
  | succ m ih =>
    obtain ⟨-, -, hs2, hc3⟩ := ih
    have hs : ∀ x : ℝ, 0 ≤ x → Real.sin x ≤ sinT (2 * (m + 1) + 1) x := fun x hx => by
      rw [show 2 * (m + 1) + 1 = 2 * m + 3 by ring]
      exact le_of_deriv_le (fun y => Real.hasDerivAt_sin y) (fun y => hasDerivAt_sinT _ y)
        (by rw [sinT_zero, Real.sin_zero]) hc3 hx
    have hc : ∀ x : ℝ, 0 ≤ x → cosT (2 * (m + 1) + 2) x ≤ Real.cos x := fun x hx =>
      le_of_deriv_le (fun y => hasDerivAt_cosT _ y) (fun y => Real.hasDerivAt_cos y)
        (by rw [cosT_zero, Real.cos_zero]) (fun y hy => by linarith [hs y hy]) hx
    have hs2' : ∀ x : ℝ, 0 ≤ x → sinT (2 * (m + 1) + 2) x ≤ Real.sin x := fun x hx =>
      le_of_deriv_le (fun y => hasDerivAt_sinT _ y) (fun y => Real.hasDerivAt_sin y)
        (by rw [sinT_zero, Real.sin_zero]) hc hx
    have hc3' : ∀ x : ℝ, 0 ≤ x → Real.cos x ≤ cosT (2 * (m + 1) + 3) x := fun x hx =>
      le_of_deriv_le (fun y => Real.hasDerivAt_cos y) (fun y => hasDerivAt_cosT _ y)
        (by rw [cosT_zero, Real.cos_zero]) (fun y hy => by linarith [hs2' y hy]) hx
    exact ⟨hs, hc, hs2', hc3'⟩

theorem cos_ge (m : ℕ) {x : ℝ} (hx : 0 ≤ x) : cosT (2 * m + 2) x ≤ Real.cos x := (bounds m).2.1 x hx

theorem sin_ge (m : ℕ) {x : ℝ} (hx : 0 ≤ x) : sinT (2 * m + 2) x ≤ Real.sin x :=
  (bounds m).2.2.1 x hx

end Fayed.Yoshida.Taylor
