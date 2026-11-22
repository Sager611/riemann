import Mathlib.Analysis.InnerProductSpace.PiL2
import Riemann.Deep.LeanGale.Row

/-!
# LeanGale: Gale's conjecture fails in every dimension n ≥ 3 (DISCOVERIES #2, discoveries/gale-r7-r8.md)

Objects (Mathlib): ℝⁿ = `EuclideanSpace ℝ (Fin n)`, `Metric.diam`, `Bornology.IsBounded`.
* `beta k V` = β_k(V): the least possible largest part diameter over partitions of V into k parts
  (NOTE §1; Voronov arXiv 2610.28828, l. 59–64, there divided by diam V).
* `sigma n` = σₙ (Gale7, the closed forms of NOTE §1 and Voronov l. 96–107).
* `GaleConjecture n`: every bounded V ⊂ ℝⁿ has β_{n+1}(V) ≤ σₙ · diam V (Gale 1953 with the simplex value σₙ;
  `GaleConjectureBall n`: with β_{n+1} of the ball of diameter 1, Voronov l. 90–92).
Results: `claim1_n7`, `claim1_n8` (the sets points7.txt, points8.txt, exact ratios), `claim2` (the family for
5 ≤ n ≤ 1000), `gale_fails_le_1000` (n = 3…1000, no hypotheses), `gale_fails` (every n ≥ 3, with Hinrichs 2002
for n > 1000), `gale_ball_fails` (the ball form, with the simplex-partition bound as a further input).
-/

namespace LeanGale

open Gale7

/-- ℝⁿ -/
abbrev E (n : ℕ) := EuclideanSpace ℝ (Fin n)

/-- β_k(V): the least possible largest part diameter over partitions of V into k parts (the classes
`V ∩ c ⁻¹' {i}` of a colouring `c`; empty parts allowed) -/
noncomputable def beta (k : ℕ) {n : ℕ} (V : Set (E n)) : ℝ :=
  ⨅ c : E n → Fin k, ⨆ i : Fin k, Metric.diam (V ∩ c ⁻¹' {i})

/-- **Gale's conjecture in ℝⁿ** (with the regular-simplex value σₙ): every bounded set V ⊂ ℝⁿ satisfies
β_{n+1}(V) ≤ σₙ · diam V. -/
def GaleConjecture (n : ℕ) : Prop :=
  ∀ V : Set (E n), Bornology.IsBounded V → beta (n + 1) V ≤ sigma n * Metric.diam V

/-- Gale's conjecture in the form "the ball is extremal" (Voronov l. 90–92): β_{n+1}(V) ≤ β_{n+1}(B) · diam V
for every bounded V ⊂ ℝⁿ, B the closed ball of diameter 1. -/
def GaleConjectureBall (n : ℕ) : Prop :=
  ∀ V : Set (E n), Bornology.IsBounded V →
    beta (n + 1) V ≤ beta (n + 1) (Metric.closedBall (0 : E n) (1 / 2)) * Metric.diam V

/-! ## Coordinates -/

/-- a coordinate vector as a point of ℝⁿ -/
def toE {n : ℕ} (x : Fin n → ℝ) : E n := WithLp.toLp 2 x

theorem toE_injective {n : ℕ} : Function.Injective (toE (n := n)) := WithLp.toLp_injective 2

theorem dist_toE {n : ℕ} (x y : Fin n → ℝ) : dist (toE x) (toE y) = eucl x y := by
  rw [EuclideanSpace.dist_eq, eucl, sqd]
  congr 1
  refine Finset.sum_congr rfl fun i _ => ?_
  simp only [toE, WithLp.ofLp_toLp, Real.dist_eq, sq_abs]

/-- a finite set of coordinate vectors as a subset of ℝⁿ -/
def XE {n : ℕ} (X : Finset (Fin n → ℝ)) : Set (E n) := toE '' (X : Set (Fin n → ℝ))

theorem isBounded_XE {n : ℕ} (X : Finset (Fin n → ℝ)) : Bornology.IsBounded (XE X) :=
  (X.finite_toSet.image _).isBounded

theorem diam_XE {n : ℕ} (X : Finset (Fin n → ℝ)) : Metric.diam (XE X) = diam X := by
  apply le_antisymm
  · apply Metric.diam_le_of_forall_dist_le (diam_nonneg X)
    rintro _ ⟨x, hx, rfl⟩ _ ⟨y, hy, rfl⟩
    rw [dist_toE]
    exact le_diam hx hy
  · by_cases hX : X.Nonempty
    · obtain ⟨x, hx, y, hy, hxy⟩ := exists_eq_diam hX
      rw [← hxy, ← dist_toE]
      exact Metric.dist_le_diam_of_mem (isBounded_XE X) ⟨x, hx, rfl⟩ ⟨y, hy, rfl⟩
    · rw [diam, dif_neg hX]
      exact Metric.diam_nonneg

/-! ## β from lower and upper bounds -/

theorem iSup_part_nonneg {k n : ℕ} (V : Set (E n)) (c : E n → Fin k) :
    0 ≤ ⨆ i : Fin k, Metric.diam (V ∩ c ⁻¹' {i}) :=
  Real.iSup_nonneg fun _ => Metric.diam_nonneg

theorem bddBelow_parts {k n : ℕ} (V : Set (E n)) :
    BddBelow (Set.range fun c : E n → Fin k => ⨆ i : Fin k, Metric.diam (V ∩ c ⁻¹' {i})) :=
  ⟨0, by rintro _ ⟨c, rfl⟩; exact iSup_part_nonneg V c⟩

theorem le_iSup_part {k n : ℕ} (V : Set (E n)) (c : E n → Fin k) (i : Fin k) :
    Metric.diam (V ∩ c ⁻¹' {i}) ≤ ⨆ i : Fin k, Metric.diam (V ∩ c ⁻¹' {i}) :=
  le_ciSup (f := fun i : Fin k => Metric.diam (V ∩ c ⁻¹' {i})) (Set.finite_range _).bddAbove i

/-- β_k(V) ≤ diam V -/
theorem beta_le_diam {k n : ℕ} (hk : 0 < k) {V : Set (E n)} (hV : Bornology.IsBounded V) :
    beta k V ≤ Metric.diam V := by
  have : Nonempty (Fin k) := ⟨⟨0, hk⟩⟩
  refine (ciInf_le (bddBelow_parts V) (fun _ => ⟨0, hk⟩)).trans (ciSup_le fun i => ?_)
  exact Metric.diam_mono Set.inter_subset_left hV

/-- a lower bound: every colouring has a same-colour pair at squared distance ≥ t -/
theorem sqrt_le_beta {k n : ℕ} (hk : 0 < k) {X : Finset (Fin n → ℝ)} {t : ℝ}
    (hlow : ∀ c : (Fin n → ℝ) → Fin k, ∃ x ∈ X, ∃ y ∈ X, c x = c y ∧ t ≤ sqd x y) :
    Real.sqrt t ≤ beta k (XE X) := by
  have : Nonempty (E n → Fin k) := ⟨fun _ => ⟨0, hk⟩⟩
  apply le_ciInf
  intro c
  obtain ⟨x, hx, y, hy, hc, ht⟩ := hlow (fun z => c (toE z))
  have hxm : toE x ∈ XE X ∩ c ⁻¹' {c (toE x)} := ⟨⟨x, hx, rfl⟩, rfl⟩
  have hym : toE y ∈ XE X ∩ c ⁻¹' {c (toE x)} := ⟨⟨y, hy, rfl⟩, hc.symm⟩
  calc Real.sqrt t ≤ eucl x y := Real.sqrt_le_sqrt ht
    _ = dist (toE x) (toE y) := (dist_toE x y).symm
    _ ≤ Metric.diam (XE X ∩ c ⁻¹' {c (toE x)}) :=
        Metric.dist_le_diam_of_mem ((isBounded_XE X).subset Set.inter_subset_left) hxm hym
    _ ≤ _ := le_iSup_part _ c _

/-- an upper bound: a colouring whose same-colour pairs are at squared distance ≤ t -/
theorem beta_le_sqrt {k n : ℕ} (hk : 0 < k) {X : Finset (Fin n → ℝ)} {t : ℝ}
    (hup : ∃ c : (Fin n → ℝ) → Fin k, ∀ x ∈ X, ∀ y ∈ X, c x = c y → sqd x y ≤ t) :
    beta k (XE X) ≤ Real.sqrt t := by
  have : Nonempty (Fin k) := ⟨⟨0, hk⟩⟩
  obtain ⟨c, hc⟩ := hup
  refine (ciInf_le (bddBelow_parts (XE X)) (fun z : E n => c (WithLp.ofLp z))).trans (ciSup_le fun i => ?_)
  apply Metric.diam_le_of_forall_dist_le (Real.sqrt_nonneg _)
  rintro _ ⟨⟨x, hx, rfl⟩, hxi⟩ _ ⟨⟨y, hy, rfl⟩, hyi⟩
  simp only [Set.mem_preimage, Set.mem_singleton_iff, toE, WithLp.ofLp_toLp] at hxi hyi
  rw [dist_toE]
  exact Real.sqrt_le_sqrt (hc x hx y hy (hxi.trans hyi.symm))

/-- σₙ < 1 -/
theorem sigma_lt_one (n : ℕ) : sigma n < 1 := by
  unfold sigma
  split_ifs
  · rw [Real.sqrt_lt' one_pos, one_pow, div_lt_one (by positivity)]
    linarith
  · rw [Real.sqrt_lt' one_pos, one_pow]
    have h1 : Real.sqrt (((n : ℝ) - 1) / ((n : ℝ) + 3)) < 1 := by
      rw [Real.sqrt_lt' one_pos, one_pow, div_lt_one (by positivity)]
      linarith
    linarith

/-! ## Counterexamples from explicit finite sets -/

/-- **The criterion**: a finite set whose colourings with n + 1 colours all have a same-colour pair at
squared distance ≥ t, with all squared distances ≤ D² and σₙ² D² < t, satisfies σₙ · diam < β_{n+1}. -/
theorem sigma_mul_diam_lt_beta {n : ℕ} {X : Finset (Fin n → ℝ)} {t D2 : ℝ}
    (hD : ∀ x ∈ X, ∀ y ∈ X, sqd x y ≤ D2)
    (hlow : ∀ c : (Fin n → ℝ) → Fin (n + 1), ∃ x ∈ X, ∃ y ∈ X, c x = c y ∧ t ≤ sqd x y)
    (hσ : sigma n ^ 2 * D2 < t) : sigma n * Metric.diam (XE X) < beta (n + 1) (XE X) := by
  have hD0 : 0 ≤ D2 := by
    obtain ⟨x, hx, _, _, _, _⟩ := hlow fun _ => 0
    exact (sqd_nonneg x x).trans (hD x hx x hx)
  have hs0 : 0 ≤ sigma n := sigma_nonneg n
  have h1 : Metric.diam (XE X) ≤ Real.sqrt D2 := by
    rw [diam_XE]
    exact diam_le (Real.sqrt_nonneg _) fun x hx y hy => Real.sqrt_le_sqrt (hD x hx y hy)
  have h2 : sigma n * Real.sqrt D2 < Real.sqrt t := by
    rw [← Real.sqrt_sq hs0, ← Real.sqrt_mul (sq_nonneg _)]
    exact Real.sqrt_lt_sqrt (mul_nonneg (sq_nonneg _) hD0) hσ
  calc sigma n * Metric.diam (XE X) ≤ sigma n * Real.sqrt D2 := mul_le_mul_of_nonneg_left h1 hs0
    _ < Real.sqrt t := h2
    _ ≤ beta (n + 1) (XE X) := sqrt_le_beta (Nat.succ_pos n) hlow

theorem not_gale_of {n : ℕ} {V : Set (E n)} (hV : Bornology.IsBounded V)
    (h : sigma n * Metric.diam V < beta (n + 1) V) : ¬ GaleConjecture n :=
  fun hG => absurd (hG V hV) (not_le.mpr h)

/-- every partition into k parts has a part of diameter > r -/
theorem exists_part_gt {k n : ℕ} (hk : 0 < k) {V : Set (E n)} {r : ℝ} (h : r < beta k V) (c : E n → Fin k) :
    ∃ i, r < Metric.diam (V ∩ c ⁻¹' {i}) := by
  have : Nonempty (Fin k) := ⟨⟨0, hk⟩⟩
  exact exists_lt_of_lt_ciSup (h.trans_le (ciInf_le (bddBelow_parts V) c))

/-- σ test for odd n: σₙ² · D < t from integers -/
theorem sigma_sq_mul_lt' {n : ℕ} {t D : ℕ} (hn : 1 ≤ n) (h : sigmaOK n t D = true) :
    sigma n ^ 2 * ((D : ℝ) / 10 ^ 16) < (t : ℝ) / 10 ^ 16 := by
  have := sigma_sq_mul_lt hn h
  rw [← mul_div_assoc]
  exact (div_lt_div_iff_of_pos_right (by positivity)).mpr this

end LeanGale
