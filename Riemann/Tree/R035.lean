import Mathlib

/-!
# Round 35 — transferring function-field mechanisms to ζ over ℚ

| mechanism | missing over ℚ | best substitute | verdict |
|-----------|----------------|-----------------|---------|
| Weil: Castelnuovo–Severi on C×C | Frobenius graph on "Spec Z × Spec Z" | Arakelov intersection theory | DEAD |
| Deligne: tensor powers + monodromy | Frobenius with pure-weight cohomology | Weil explicit-formula criterion | CIRCULAR |
| Stepanov–Bombieri | Frobenius, point counts over extensions | explicit formula | CIRCULAR |
| Faltings–Hriljac Hodge index | the surface Spec Z × Spec Z | Yuan–Zhang 2013 adelic Hodge index | partial only |
| standard conjectures | global Frobenius | Connes trace formula, Gillet–Soulé | CIRCULAR |
| Lefschetz trace + Frobenius flow | cohomology realising the zeros | Connes 1999, Deninger | CIRCULAR |
| F₁ / arithmetic site | positive intersection pairing | Connes–Consani Riemann–Roch for Spec Z | partial only |
| rationality (Dwork) | finite-dimensional cohomology | regularized determinants | partial only |

The purely analytic step — Deligne's tensor-power trick — is machine-checked below. Over ℚ it lacks
its input: a family of L-functions whose "eigenvalues" are the k-th powers of ζ's.
-/

namespace Fayed.R35

/-- Tensor-power trick: a bound x^k ≤ C·y^k for every k forces x ≤ y. -/
theorem tensor_power_trick {x y C : ℝ} (hy : 0 < y) (h : ∀ k : ℕ, x ^ k ≤ C * y ^ k) :
    x ≤ y := by
  by_contra hxy
  rw [not_le] at hxy
  have h1 : 1 < x / y := (one_lt_div hy).mpr hxy
  obtain ⟨k, hk⟩ := ((tendsto_pow_atTop_atTop_of_one_lt h1).eventually_gt_atTop C).exists
  have hk' : x ^ k / y ^ k ≤ C := (div_le_iff₀ (pow_pos hy k)).mpr (h k)
  rw [div_pow] at hk
  linarith

end Fayed.R35

#print axioms Fayed.R35.tensor_power_trick
