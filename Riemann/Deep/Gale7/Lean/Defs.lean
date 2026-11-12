import Riemann.Deep.Gale7.Lean.Fam
import Mathlib.Analysis.Real.Sqrt

/-!
# Gale7: definitions (NOTE §1)

Points of ℝⁿ are coordinate vectors `Fin n → ℝ` with the Euclidean distance
`eucl x y = √(Σₖ (x k - y k)²)` (Mathlib's `dist` on `Fin n → ℝ` is the sup distance, so it is not used).
`diam X` is the largest distance in the finite set `X`; `sigma n` is σₙ of NOTE §1; `GaleSplit n X` says that
`X` splits into `n + 1` parts (colour classes, possibly empty) of diameter ≤ σₙ · diam X.
-/

namespace Gale7

/-- Euclidean distance in ℝⁿ. -/
noncomputable def eucl {n : ℕ} (x y : Fin n → ℝ) : ℝ := Real.sqrt (sqd x y)

/-- Diameter of a finite set: its largest distance (0 for the empty set). -/
noncomputable def diam {n : ℕ} (X : Finset (Fin n → ℝ)) : ℝ :=
  if h : X.Nonempty then (X ×ˢ X).sup' (h.product h) (fun z => eucl z.1 z.2) else 0

/-- σₙ (NOTE §1; verify.py `above_ball`): the part diameter of the regular-simplex partition of a ball of
diameter 1, σₙ² = (n+1)/(n+2) for even n and σₙ² = 1/2 + √((n−1)/(n+3))/2 for odd n. -/
noncomputable def sigma (n : ℕ) : ℝ :=
  if n % 2 = 0 then Real.sqrt (((n : ℝ) + 1) / ((n : ℝ) + 2))
  else Real.sqrt (1 / 2 + Real.sqrt (((n : ℝ) - 1) / ((n : ℝ) + 3)) / 2)

/-- Gale's property for a finite `X ⊂ ℝⁿ`: a partition into `n + 1` parts (the classes of `c`, empty parts
allowed) such that every part has diameter at most σₙ · diam X. -/
def GaleSplit (n : ℕ) (X : Finset (Fin n → ℝ)) : Prop :=
  ∃ c : (Fin n → ℝ) → Fin (n + 1), ∀ k, diam (X.filter fun x => c x = k) ≤ sigma n * diam X

theorem eucl_nonneg {n : ℕ} (x y : Fin n → ℝ) : 0 ≤ eucl x y := Real.sqrt_nonneg _

theorem le_diam {n : ℕ} {X : Finset (Fin n → ℝ)} {x y : Fin n → ℝ} (hx : x ∈ X) (hy : y ∈ X) :
    eucl x y ≤ diam X := by
  have hX : X.Nonempty := ⟨x, hx⟩
  have hm : (x, y) ∈ X ×ˢ X := Finset.mem_product.mpr ⟨hx, hy⟩
  rw [diam, dif_pos hX]
  exact Finset.le_sup' (fun z : (Fin n → ℝ) × (Fin n → ℝ) => eucl z.1 z.2) hm

theorem diam_le {n : ℕ} {X : Finset (Fin n → ℝ)} {M : ℝ} (hM : 0 ≤ M)
    (h : ∀ x ∈ X, ∀ y ∈ X, eucl x y ≤ M) : diam X ≤ M := by
  unfold diam
  split_ifs with hX
  · apply Finset.sup'_le
    rintro ⟨x, y⟩ hxy
    rw [Finset.mem_product] at hxy
    exact h x hxy.1 y hxy.2
  · exact hM

theorem diam_nonneg {n : ℕ} (X : Finset (Fin n → ℝ)) : 0 ≤ diam X := by
  by_cases hX : X.Nonempty
  · obtain ⟨x, hx⟩ := hX
    exact (eucl_nonneg x x).trans (le_diam hx hx)
  · rw [diam, dif_neg hX]

/-- The diameter is attained. -/
theorem exists_eq_diam {n : ℕ} {X : Finset (Fin n → ℝ)} (hX : X.Nonempty) :
    ∃ x ∈ X, ∃ y ∈ X, eucl x y = diam X := by
  rw [diam, dif_pos hX]
  obtain ⟨z, hz, hz'⟩ := Finset.exists_mem_eq_sup' (hX.product hX) (fun z => eucl z.1 z.2)
  rw [Finset.mem_product] at hz
  exact ⟨z.1, hz.1, z.2, hz.2, hz'.symm⟩

theorem sigma_nonneg (n : ℕ) : 0 ≤ sigma n := by
  unfold sigma
  split_ifs <;> exact Real.sqrt_nonneg _

/-- `GaleSplit` in pair form: same-colour points are at distance ≤ σₙ · diam X. -/
theorem galeSplit_iff (n : ℕ) (X : Finset (Fin n → ℝ)) :
    GaleSplit n X ↔ ∃ c : (Fin n → ℝ) → Fin (n + 1),
      ∀ x ∈ X, ∀ y ∈ X, c x = c y → eucl x y ≤ sigma n * diam X := by
  constructor
  · rintro ⟨c, hc⟩
    refine ⟨c, fun x hx y hy hxy => ?_⟩
    have hx' : x ∈ X.filter (fun z => c z = c x) := Finset.mem_filter.mpr ⟨hx, rfl⟩
    have hy' : y ∈ X.filter (fun z => c z = c x) := Finset.mem_filter.mpr ⟨hy, hxy.symm⟩
    exact (le_diam hx' hy').trans (hc (c x))
  · rintro ⟨c, hc⟩
    refine ⟨c, fun k => diam_le (mul_nonneg (sigma_nonneg n) (diam_nonneg X)) ?_⟩
    intro x hx y hy
    rw [Finset.mem_filter] at hx hy
    exact hc x hx.1 y hy.1 (hx.2.trans hy.2.symm)

theorem sigma_odd_sq {n : ℕ} (hn : n % 2 = 1) (h1 : 1 ≤ n) :
    sigma n ^ 2 = 1 / 2 + Real.sqrt (((n : ℝ) - 1) / ((n : ℝ) + 3)) / 2 := by
  have hn' : ¬ n % 2 = 0 := by omega
  have h1' : (1 : ℝ) ≤ n := by exact_mod_cast h1
  rw [sigma, if_neg hn', Real.sq_sqrt]
  have : 0 ≤ Real.sqrt (((n : ℝ) - 1) / ((n : ℝ) + 3)) := Real.sqrt_nonneg _
  linarith

theorem sigma_even_sq {n : ℕ} (hn : n % 2 = 0) :
    sigma n ^ 2 = ((n : ℝ) + 1) / ((n : ℝ) + 2) := by
  have h0 : (0 : ℝ) ≤ n := Nat.cast_nonneg n
  rw [sigma, if_pos hn, Real.sq_sqrt]
  exact div_nonneg (by linarith) (by linarith)

end Gale7
