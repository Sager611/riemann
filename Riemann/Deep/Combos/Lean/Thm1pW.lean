module

public import Riemann.Deep.Combos.Lean.Thm1pI
public import Riemann.Deep.Combos.Lean.Thm1pL

@[expose] public section

/-! # Theorem 1′(a): the window bound (RCApps2 l.33, 57) and the last line -/

open Real Finset

namespace RingClass

/-- every window holds ≤ c₀ (2/ε₀) W^{ε₀} zeros, W = 9|D|T -/
lemma RCZeros.window {X : Type*} [Fintype X] (d : RCZeros X) {c₀ c₁ : ℝ} {CI CII Q0 : ℝ → ℝ}
    (h : d.Hyp c₀ c₁ CI CII Q0) (ε₀ : ℝ) (he : 0 < ε₀) (χ : X) (t : ℝ) :
    (∑ ρ ∈ d.win χ t, (d.mult χ ρ : ℝ)) ≤ c₀ * (2 / ε₀) * (9 * d.absD * (|t| + 3)) ^ ε₀ := by
  have hT : 1 ≤ |t| + 3 := by have := abs_nonneg t; linarith only [this]
  have hD : 1 ≤ d.absD := le_trans (by norm_num) d.absD_ge
  have hA := d.A_ge χ
  have hT2 : 0 < (|t| + 3) ^ 2 := pow_pos (lt_of_lt_of_le one_pos hT) 2
  have h1 : Real.log (d.A χ * (|t| + 3) ^ 2) ≤ Real.log (d.absD * (|t| + 3) ^ 2) :=
    Real.log_le_log (mul_pos (lt_of_lt_of_le one_pos hA) hT2)
      (mul_le_mul_of_nonneg_right (d.A_le χ) hT2.le)
  have h2 := window_log d.absD (|t| + 3) ε₀ hD hT he
  have h3 := h.F4 χ t
  have h4 := mul_le_mul_of_nonneg_left (h1.trans h2) h.c₀_nonneg
  have e : c₀ * (2 / ε₀ * (9 * d.absD * (|t| + 3)) ^ ε₀)
      = c₀ * (2 / ε₀) * (9 * d.absD * (|t| + 3)) ^ ε₀ := by ring
  linarith only [h3, h4, e]

/-- the last line: from N ≤ A K P with K ≤ K_tot, P ≤ 81 R -/
lemma fin_step (N A K Ktot P R : ℝ) (hA : 0 ≤ A) (hP : 0 ≤ P) (hK : K ≤ Ktot) (hKt : 0 ≤ Ktot)
    (hN : N ≤ A * K * P) (hPR : P ≤ 81 * R) : N ≤ 81 * (A * Ktot) * R := by
  have g1 : A * K * P ≤ A * Ktot * P :=
    mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_left hK hA) hP
  have g2 : A * Ktot * P ≤ A * Ktot * (81 * R) := mul_le_mul_of_nonneg_left hPR (mul_nonneg hA hKt)
  have e : A * Ktot * (81 * R) = 81 * (A * Ktot) * R := by ring
  linarith only [hN, g1, g2, e]

end RingClass

