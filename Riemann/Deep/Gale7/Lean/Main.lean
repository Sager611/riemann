import Riemann.Deep.Gale7.Lean.N7
import Riemann.Deep.Gale7.Lean.N8

/-!
# Gale7: main theorems

V₇ has 69 points and V₈ has 85 points (distinct vertices have a descriptor with `neD = true`, whose class
values are positive by `decide +kernel`); Gale's conjecture fails in ℝ⁷ and ℝ⁸; every partition of V₇ (V₈)
into 8 (9) parts has a part of diameter > 0.9433983 · diam V₇ (> 0.949947 · diam V₈), while σ₇ < 0.9419652
and σ₈ < 0.9486834.
-/

namespace Gale7

open Vx

theorem pt7_injective : Function.Injective pt7 := by
  intro v w h
  by_contra hne
  have h1 : (0 : ℝ) < sqd (pt7 v) (pt7 w) := by
    rw [sqd7]
    exact Rat.cast_pos.mpr (cls7_pos _ _ _ _ (ne_neD hne))
  rw [h, sqd_self] at h1
  exact lt_irrefl _ h1

theorem pt8_injective : Function.Injective pt8 := by
  intro v w h
  by_contra hne
  have h1 : (0 : ℝ) < sqd (pt8 v) (pt8 w) := by
    rw [sqd8]
    exact Rat.cast_pos.mpr (cls8_pos _ _ _ _ (ne_neD hne))
  rw [h, sqd_self] at h1
  exact lt_irrefl _ h1

theorem card_Vx (p q : ℕ) : Fintype.card (Vx p q) = p + q + 3 * (p * q) := by
  have e : Vx p q ≃ Fin p ⊕ Fin q ⊕ (Fin p × Fin q) ⊕ (Fin p × Fin q) ⊕ (Fin p × Fin q) :=
    { toFun := fun v => match v with
        | P i => .inl i
        | Q j => .inr (.inl j)
        | A i j => .inr (.inr (.inl (i, j)))
        | B i j => .inr (.inr (.inr (.inl (i, j))))
        | C i j => .inr (.inr (.inr (.inr (i, j))))
      invFun := fun s => match s with
        | .inl i => P i
        | .inr (.inl j) => Q j
        | .inr (.inr (.inl (i, j))) => A i j
        | .inr (.inr (.inr (.inl (i, j)))) => B i j
        | .inr (.inr (.inr (.inr (i, j)))) => C i j
      left_inv := fun v => by cases v <;> rfl
      right_inv := fun s => by
        rcases s with i | j | ⟨i, j⟩ | ⟨i, j⟩ | ⟨i, j⟩ <;> rfl }
  rw [Fintype.card_congr e]
  simp only [Fintype.card_sum, Fintype.card_prod, Fintype.card_fin]
  ring

/-- V₇ has 4 + 5 + 3·4·5 = 69 points -/
theorem card_X7 : X7.card = 69 := by
  rw [X7, Finset.card_image_of_injective _ pt7_injective, Finset.card_univ, card_Vx]

/-- V₈ has 5 + 5 + 3·5·5 = 85 points -/
theorem card_X8 : X8.card = 85 := by
  rw [X8, Finset.card_image_of_injective _ pt8_injective, Finset.card_univ, card_Vx]

/-- **Gale's conjecture is false in ℝ⁷.** -/
theorem gale_false_7 : ∃ X : Finset (Fin 7 → ℝ), X.card = 69 ∧ ¬ GaleSplit 7 X :=
  ⟨X7, card_X7, not_galeSplit_X7⟩

/-- **Gale's conjecture is false in ℝ⁸.** -/
theorem gale_false_8 : ∃ X : Finset (Fin 8 → ℝ), X.card = 85 ∧ ¬ GaleSplit 8 X :=
  ⟨X8, card_X8, not_galeSplit_X8⟩

/-- β₈(V₇) > 0.9433983 · diam V₇: every partition into 8 parts has a part with two points farther apart
(and 0.9433983 > 0.9419652 > σ₇, `sigma7_lt`). -/
theorem beta7_gt (c : (Fin 7 → ℝ) → Fin 8) :
    ∃ x ∈ X7, ∃ y ∈ X7, c x = c y ∧ 0.9433983 * diam X7 < eucl x y := by
  obtain ⟨x, hx, y, hy, hc, ht⟩ := X7_long c
  refine ⟨x, hx, y, hy, hc, ?_⟩
  rw [diam_X7]
  have h1 : (0.9433983 : ℝ) * Real.sqrt (D7 : ℝ) < Real.sqrt (t7 : ℝ) := by
    rw [Real.lt_sqrt (by positivity), mul_pow, Real.sq_sqrt (by norm_num [D7])]
    norm_num [D7, t7]
  exact h1.trans_le (Real.sqrt_le_sqrt ht)

/-- β₉(V₈) > 0.949947 · diam V₈ (and 0.949947 > 0.9486834 > σ₈, `sigma8_lt`). -/
theorem beta8_gt (c : (Fin 8 → ℝ) → Fin 9) :
    ∃ x ∈ X8, ∃ y ∈ X8, c x = c y ∧ 0.949947 * diam X8 < eucl x y := by
  obtain ⟨x, hx, y, hy, hc, ht⟩ := X8_long c
  refine ⟨x, hx, y, hy, hc, ?_⟩
  rw [diam_X8]
  have h1 : (0.949947 : ℝ) * Real.sqrt (D8 : ℝ) < Real.sqrt (t8 : ℝ) := by
    rw [Real.lt_sqrt (by positivity), mul_pow, Real.sq_sqrt (by norm_num [D8])]
    norm_num [D8, t8]
  exact h1.trans_le (Real.sqrt_le_sqrt ht)

end Gale7
