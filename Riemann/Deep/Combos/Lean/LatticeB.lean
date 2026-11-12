import Riemann.Deep.Combos.Lean.LatticeA

/-! # Lemma 1, the row decomposition (NOTE l.71–74) -/

open Finset

namespace RingClass

/-- the covering set {(±1, 0)} ∪ ⋃_{y ∈ rows V} row y -/
noncomputable def Tset (a b : ℤ) (p L V : ℝ) : Finset (ℤ × ℤ) :=
  {(1, 0), (-1, 0)} ∪ (rows V).biUnion (row a b p L)

lemma card_Tset_le (a b : ℤ) (p L V : ℝ) (hL0 : 0 ≤ L) (hV0 : 0 ≤ V) :
    ((Tset a b p L V).card : ℝ) ≤ 2 + 2 * V * (L + 1) := by
  have h1 : (Tset a b p L V).card ≤ ({(1, 0), (-1, 0)} : Finset (ℤ × ℤ)).card
      + ((rows V).biUnion (row a b p L)).card := Finset.card_union_le _ _
  have h2 : ({(1, 0), (-1, 0)} : Finset (ℤ × ℤ)).card ≤ 2 := Finset.card_le_two
  have h3 : ((rows V).biUnion (row a b p L)).card ≤ ∑ y ∈ rows V, (row a b p L y).card :=
    Finset.card_biUnion_le
  have h4 : ∀ y ∈ rows V, ((row a b p L y).card : ℝ) ≤ L + 1 := fun y _ => by
    calc ((row a b p L y).card : ℝ) ≤ (Icc ⌈lft a b p y⌉ ⌊lft a b p y + L⌋).card := by
          exact_mod_cast Finset.card_image_le
      _ ≤ L + 1 := card_Icc_le _ _ hL0
  have h5 := card_rows_le V hV0
  have h6 : (∑ y ∈ rows V, ((row a b p L y).card : ℝ)) ≤ (rows V).card * (L + 1) := by
    calc _ ≤ ∑ y ∈ rows V, (L + 1) := Finset.sum_le_sum h4
      _ = _ := by rw [Finset.sum_const, nsmul_eq_mul]
  have h8 := h1.trans (Nat.add_le_add h2 h3)
  have h9 : (((Tset a b p L V).card : ℕ) : ℝ) ≤ ((2 + ∑ y ∈ rows V, (row a b p L y).card : ℕ) : ℝ) := by
    exact_mod_cast h8
  have h7 : ((Tset a b p L V).card : ℝ) ≤ 2 + ∑ y ∈ rows V, ((row a b p L y).card : ℝ) := by
    rw [Nat.cast_add, Nat.cast_sum] at h9; exact_mod_cast h9
  have h10 : ((rows V).card : ℝ) * (L + 1) ≤ 2 * V * (L + 1) :=
    mul_le_mul_of_nonneg_right h5 (by linarith)
  linarith

/-- the row y = 0: primitive points are (±1, 0) -/
lemma row_zero_mem (x : ℤ) (hg : Int.gcd x 0 = 1) :
    (x, (0 : ℤ)) ∈ ({(1, 0), (-1, 0)} : Finset (ℤ × ℤ)) := by
  have h1 : x.natAbs = 1 := by rw [Int.gcd_zero_right] at hg; exact hg
  rcases Int.natAbs_eq_iff.mp h1 with h | h
  · rw [h]; exact Finset.mem_insert_self _ _
  · rw [h]; exact Finset.mem_insert_of_mem (Finset.mem_singleton_self _)

lemma PQset_subset_Tset (a b c : ℤ) (ha : 0 < a) (hb : |b| ≤ a) (hac : a ≤ c) (N : ℝ)
    (p q : ℝ) (hp0 : 0 ≤ p) (hq0 : 0 < q) (hp : a * N ≤ p ^ 2)
    (hq : q ^ 2 ≤ 4 * a * c - b ^ 2) :
    PQset a b c N ⊆ ↑(Tset a b p (2 * p / a) (2 * p / q)) := by
  have haR : (0 : ℝ) < a := by exact_mod_cast ha
  have hΔ3 := disc_facts (a : ℝ) b c haR (by exact_mod_cast hb) (by exact_mod_cast hac)
  have hΔpos : (0 : ℝ) < 4 * a * c - b ^ 2 := by
    have : (0 : ℝ) < 3 * (a : ℝ) ^ 2 := by positivity
    linarith
  have hV0 : 0 ≤ 2 * p / q := by positivity
  rintro ⟨x, y⟩ ⟨hg, hQ⟩
  simp only at hg hQ
  have key : (4 * a * (Qf a b c x y : ℝ)) =
      (2 * a * x + b * y) ^ 2 + (4 * a * c - b ^ 2) * y ^ 2 := by
    simp only [Qf]; push_cast; ring
  have h4 : (4 * a) * (Qf a b c x y : ℝ) ≤ (4 * a) * N :=
    mul_le_mul_of_nonneg_left hQ (by positivity)
  have hDy : 0 ≤ (4 * (a : ℝ) * c - b ^ 2) * (y : ℝ) ^ 2 := mul_nonneg hΔpos.le (sq_nonneg _)
  have hsq : 0 ≤ ((2 * a * x + b * y : ℝ)) ^ 2 := sq_nonneg _
  have hqy : q ^ 2 * (y : ℝ) ^ 2 ≤ (4 * (a : ℝ) * c - b ^ 2) * (y : ℝ) ^ 2 :=
    mul_le_mul_of_nonneg_right hq (sq_nonneg _)
  unfold Tset
  rw [Finset.coe_union, Set.mem_union]
  by_cases hy : y = 0
  · left
    rw [hy] at hg ⊢
    exact row_zero_mem x hg
  · right
    rw [Finset.mem_coe, Finset.mem_biUnion]
    refine ⟨y, ?_, ?_⟩
    · rw [rows, Finset.mem_erase]
      refine ⟨hy, ?_⟩
      rw [Finset.mem_Icc, Int.ceil_le, Int.le_floor]
      have hy2 : (y : ℝ) ^ 2 ≤ (2 * p / q) ^ 2 := by
        rw [div_pow, le_div_iff₀ (by positivity), mul_pow]
        linarith
      have h := abs_le_of_sq_le_sq' hy2 hV0
      constructor <;> linarith [h.1, h.2]
    · rw [row, Finset.mem_image]
      refine ⟨x, ?_, rfl⟩
      rw [Finset.mem_Icc, Int.ceil_le, Int.le_floor]
      have hs : (2 * a * x + b * y : ℝ) ^ 2 ≤ (2 * p) ^ 2 := by
        rw [mul_pow]; linarith
      have h := abs_le_of_sq_le_sq' hs (by positivity)
      have e : lft a b p y + 2 * p / a = (-(b : ℝ) * y + 2 * p) / (2 * a) := by
        rw [lft]; field_simp; ring
      constructor
      · rw [lft, div_le_iff₀ (by positivity)]; linarith [h.1]
      · rw [e, le_div_iff₀ (by positivity)]; linarith [h.2]

/-- The row decomposition of NOTE l.71–74: finiteness and P_Q(N) ≤ 2 + 2V(L + 1) with
V = 2p/q, L = 2p/a, for any p ≥ 0, q > 0 with p² ≥ aN, q² ≤ |D| = 4ac − b². -/
theorem PQset_subset_card (a b c : ℤ) (ha : 0 < a) (hb : |b| ≤ a) (hac : a ≤ c) (N : ℝ)
    (p q : ℝ) (hp0 : 0 ≤ p) (hq0 : 0 < q) (hp : a * N ≤ p ^ 2)
    (hq : q ^ 2 ≤ 4 * a * c - b ^ 2) :
    (PQset a b c N).Finite ∧
      ((PQset a b c N).ncard : ℝ) ≤ 2 + 2 * (2 * p / q) * (2 * p / a + 1) := by
  have haR : (0 : ℝ) < a := by exact_mod_cast ha
  have hsub := PQset_subset_Tset a b c ha hb hac N p q hp0 hq0 hp hq
  have hc := card_Tset_le a b p (2 * p / a) (2 * p / q) (by positivity) (by positivity)
  refine ⟨(Tset a b p (2 * p / a) (2 * p / q)).finite_toSet.subset hsub, ?_⟩
  have h1 : (PQset a b c N).ncard ≤ ((Tset a b p (2 * p / a) (2 * p / q) : Finset (ℤ × ℤ)) :
      Set (ℤ × ℤ)).ncard := Set.ncard_le_ncard hsub (Finset.finite_toSet _)
  rw [Set.ncard_coe_finset] at h1
  have h2 : ((PQset a b c N).ncard : ℝ) ≤ (Tset a b p (2 * p / a) (2 * p / q)).card := by
    exact_mod_cast h1
  linarith

end RingClass

