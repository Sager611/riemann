import Riemann.Deep.Gale7.Lean.Defs

/-!
# Gale7: the reduction to the colouring lemma (NOTE §3, Theorem)

`not_galeSplit_of_H`: if all edges of `H(p, q)` have squared length ≥ t, all squared distances are ≤ D², and
σₙ²·D² < t, then the point set is not Gale-splittable (no partition into n + 1 parts of diameter ≤ σₙ · diam).
-/

namespace Gale7

/-- Every partition of the point set into `n + 1` parts has a part with two points at squared distance
≥ t (NOTE §3, Theorem; uses `Vx.not_colourable`). -/
theorem long_pair {n p q : ℕ} (hpq : n + 2 ≤ p + q) (pt : Vx p q → Fin n → ℝ) {t : ℝ}
    (hH : ∀ v w, Vx.Adj v w → t ≤ sqd (pt v) (pt w)) (c : (Fin n → ℝ) → Fin (n + 1)) :
    ∃ x ∈ Finset.univ.image pt, ∃ y ∈ Finset.univ.image pt, c x = c y ∧ t ≤ sqd x y := by
  obtain ⟨v, w, hvw, hc⟩ := Vx.not_colourable (K := Fin (n + 1))
    (by rw [Fintype.card_fin]; omega) (fun v => c (pt v))
  exact ⟨pt v, Finset.mem_image_of_mem pt (Finset.mem_univ v), pt w,
    Finset.mem_image_of_mem pt (Finset.mem_univ w), hc, hH v w hvw⟩

theorem diam_image_le {n p q : ℕ} (pt : Vx p q → Fin n → ℝ) {D2 : ℝ}
    (hD : ∀ v w, sqd (pt v) (pt w) ≤ D2) : diam (Finset.univ.image pt) ≤ Real.sqrt D2 := by
  apply diam_le (Real.sqrt_nonneg _)
  intro x hx y hy
  obtain ⟨v, _, rfl⟩ := Finset.mem_image.mp hx
  obtain ⟨w, _, rfl⟩ := Finset.mem_image.mp hy
  exact Real.sqrt_le_sqrt (hD v w)

/-- **Reduction** (NOTE §3, Theorem): the colouring lemma plus three metric facts refute Gale's property. -/
theorem not_galeSplit_of_H {n p q : ℕ} (hpq : n + 2 ≤ p + q) (pt : Vx p q → Fin n → ℝ)
    {t D2 : ℝ} (hD0 : 0 ≤ D2) (hD : ∀ v w, sqd (pt v) (pt w) ≤ D2)
    (hH : ∀ v w, Vx.Adj v w → t ≤ sqd (pt v) (pt w)) (hσ : sigma n ^ 2 * D2 < t) :
    ¬ GaleSplit n (Finset.univ.image pt) := by
  rw [galeSplit_iff]
  rintro ⟨c, hc⟩
  obtain ⟨x, hx, y, hy, hxy, ht⟩ := long_pair hpq pt hH c
  have h1 := hc x hx y hy hxy
  have hs0 : 0 ≤ sigma n := sigma_nonneg n
  have h2 : sigma n * diam (Finset.univ.image pt) ≤ sigma n * Real.sqrt D2 :=
    mul_le_mul_of_nonneg_left (diam_image_le pt hD) hs0
  have h3 : sigma n * Real.sqrt D2 < Real.sqrt t := by
    rw [← Real.sqrt_sq hs0, ← Real.sqrt_mul (sq_nonneg _)]
    exact Real.sqrt_lt_sqrt (mul_nonneg (sq_nonneg _) hD0) hσ
  have h4 : Real.sqrt t ≤ eucl x y := Real.sqrt_le_sqrt ht
  linarith

end Gale7
