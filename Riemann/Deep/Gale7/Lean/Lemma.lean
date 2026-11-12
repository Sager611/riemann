import Mathlib.Combinatorics.Enumerative.DoubleCounting
import Mathlib.Data.Fintype.Prod
import Mathlib.Data.Fintype.Sigma
import Mathlib.Data.Fintype.Sum
import Mathlib.Tactic.DeriveFintype
import Mathlib.Tactic.Linarith

/-!
# Gale7: the colouring lemma (NOTE §3)

`H(p, q)`: vertices `P i`, `Q j`, `A i j`, `B i j`, `C i j` (`i : Fin p` rows, `j : Fin q` columns), edges as
in NOTE §3. `Vx.card_le`: every proper colouring of `H(p, q)` uses at least `p + q` colours; hence
`Vx.lemma_H`: no proper colouring with `p + q - 1` colours. The proof follows NOTE §3 (classes `X(i)`,
`Y(j)`, `XY(i,j)`, `L`; `R₀`, `C₀`; the line cover of the `L` classes; the final `r₀ = 1` / `c₀ = 1` cases).
-/

namespace Gale7

/-- Vertices of `H(p, q)` (NOTE §3). -/
inductive Vx (p q : ℕ) : Type where
  | P : Fin p → Vx p q
  | Q : Fin q → Vx p q
  | A : Fin p → Fin q → Vx p q
  | B : Fin p → Fin q → Vx p q
  | C : Fin p → Fin q → Vx p q
  deriving DecidableEq, Fintype

namespace Vx

variable {p q : ℕ}

/-- Edges of `H(p, q)` (NOTE §3), for both orders of the end points:
`PᵢP_k, A_ijA_kl, P_kA_ij, P_kC_ij, A_ijC_kl` (i ≠ k); `Q_jQ_l, B_ijB_kl, Q_lB_ij, Q_lC_ij, B_ijC_kl` (j ≠ l);
`C_ijC_kl` (i ≠ k and j ≠ l); `PᵢB_ij`, `Q_jA_ij` (all i, j). -/
def Adj : Vx p q → Vx p q → Prop
  | P i, P k => i ≠ k
  | P _, Q _ => False
  | P k, A i _ => k ≠ i
  | P k, B i _ => k = i
  | P k, C i _ => k ≠ i
  | Q _, P _ => False
  | Q j, Q l => j ≠ l
  | Q l, A _ j => l = j
  | Q l, B _ j => l ≠ j
  | Q l, C _ j => l ≠ j
  | A i _, P k => k ≠ i
  | A _ j, Q l => l = j
  | A i _, A k _ => i ≠ k
  | A _ _, B _ _ => False
  | A i _, C k _ => i ≠ k
  | B i _, P k => k = i
  | B _ j, Q l => l ≠ j
  | B _ _, A _ _ => False
  | B _ j, B _ l => j ≠ l
  | B _ j, C _ l => j ≠ l
  | C i _, P k => k ≠ i
  | C _ j, Q l => l ≠ j
  | C k _, A i _ => i ≠ k
  | C _ l, B _ j => j ≠ l
  | C i j, C k l => i ≠ k ∧ j ≠ l

theorem adj_symm {x y : Vx p q} (h : Adj x y) : Adj y x := by
  cases x <;> cases y <;>
    first
    | exact h
    | exact Ne.symm h
    | exact ⟨Ne.symm h.1, Ne.symm h.2⟩

theorem adj_irrefl (x : Vx p q) : ¬ Adj x x := by
  cases x <;> intro h
  · exact h rfl
  · exact h rfl
  · exact h rfl
  · exact h rfl
  · exact h.1 rfl

section Defs

variable {K : Type*} (col : Vx p q → K)

/-- The colour `κ` meets the row block `X_i = {P_i} ∪ {A_ij : j}`. -/
def InX (κ : K) (i : Fin p) : Prop := col (P i) = κ ∨ ∃ j, col (A i j) = κ

/-- The colour `κ` meets the column block `Y_j = {Q_j} ∪ {B_ij : i}`. -/
def InY (κ : K) (j : Fin q) : Prop := col (Q j) = κ ∨ ∃ i, col (B i j) = κ

/-- All `C`-vertices of colour `κ` lie in one row. -/
def RowT (κ : K) : Prop := ∀ i j i' j', col (C i j) = κ → col (C i' j') = κ → i = i'

variable [Fintype K]
open Classical

/-- classes `X(i)`: meet a row block, no column block -/
noncomputable def Kx : Finset K := Finset.univ.filter fun κ => (∃ i, InX col κ i) ∧ ¬ ∃ j, InY col κ j
/-- classes `Y(j)` -/
noncomputable def Ky : Finset K := Finset.univ.filter fun κ => (¬ ∃ i, InX col κ i) ∧ ∃ j, InY col κ j
/-- classes `XY(i,j)` -/
noncomputable def Kxy : Finset K := Finset.univ.filter fun κ => (∃ i, InX col κ i) ∧ ∃ j, InY col κ j
/-- classes `L` (no block; also the unused colours) -/
noncomputable def Kl : Finset K := Finset.univ.filter fun κ => (¬ ∃ i, InX col κ i) ∧ ¬ ∃ j, InY col κ j
/-- `R₀`: rows without an `X(i)` class -/
noncomputable def R0 : Finset (Fin p) := Finset.univ.filter fun i => ∀ κ ∈ Kx col, ¬ InX col κ i
/-- `C₀`: columns without a `Y(j)` class -/
noncomputable def C0 : Finset (Fin q) := Finset.univ.filter fun j => ∀ κ ∈ Ky col, ¬ InY col κ j
/-- the row of an `L` class whose `C`-vertices lie in one row -/
noncomputable def rowsL (κ : K) : Finset (Fin p) :=
  if RowT col κ then Finset.univ.filter fun i => ∃ j, col (C i j) = κ else ∅
/-- the column of any other `L` class -/
noncomputable def colsL (κ : K) : Finset (Fin q) :=
  if RowT col κ then ∅ else Finset.univ.filter fun j => ∃ i, col (C i j) = κ
/-- rows met by row-type `L` classes -/
noncomputable def R1 : Finset (Fin p) := (Kl col).biUnion (rowsL col)
/-- columns met by column-type `L` classes -/
noncomputable def C1 : Finset (Fin q) := (Kl col).biUnion (colsL col)

theorem card_split :
    (Kx col).card + (Ky col).card + (Kxy col).card + (Kl col).card ≤ Fintype.card K := by
  have d1 : Disjoint (Kx col) (Ky col) := by
    rw [Finset.disjoint_left]
    intro κ h1 h2
    simp only [Kx, Ky, Finset.mem_filter, Finset.mem_univ, true_and] at h1 h2
    exact h2.1 h1.1
  have d2 : Disjoint (Kx col ∪ Ky col) (Kxy col) := by
    rw [Finset.disjoint_left]
    intro κ h1 h2
    simp only [Kx, Ky, Kxy, Finset.mem_union, Finset.mem_filter, Finset.mem_univ, true_and] at h1 h2
    rcases h1 with h1 | h1
    · exact h1.2 h2.2
    · exact h1.1 h2.1
  have d3 : Disjoint (Kx col ∪ Ky col ∪ Kxy col) (Kl col) := by
    rw [Finset.disjoint_left]
    intro κ h1 h2
    simp only [Kx, Ky, Kxy, Kl, Finset.mem_union, Finset.mem_filter, Finset.mem_univ,
      true_and] at h1 h2
    rcases h1 with (h1 | h1) | h1
    · exact h2.1 h1.1
    · exact h2.2 h1.2
    · exact h2.1 h1.1
  rw [← Finset.card_union_of_disjoint d1, ← Finset.card_union_of_disjoint d2,
    ← Finset.card_union_of_disjoint d3]
  exact Finset.card_le_univ _

theorem mem_Kxy_of_R0 {i : Fin p} (hi : i ∈ R0 col) {κ : K} (h : InX col κ i) : κ ∈ Kxy col := by
  simp only [R0, Finset.mem_filter, Finset.mem_univ, true_and] at hi
  simp only [Kxy, Finset.mem_filter, Finset.mem_univ, true_and]
  refine ⟨⟨i, h⟩, ?_⟩
  by_contra hY
  exact hi κ (by simp only [Kx, Finset.mem_filter, Finset.mem_univ, true_and]; exact ⟨⟨i, h⟩, hY⟩) h

theorem mem_Kxy_of_C0 {j : Fin q} (hj : j ∈ C0 col) {κ : K} (h : InY col κ j) : κ ∈ Kxy col := by
  simp only [C0, Finset.mem_filter, Finset.mem_univ, true_and] at hj
  simp only [Kxy, Finset.mem_filter, Finset.mem_univ, true_and]
  refine ⟨?_, ⟨j, h⟩⟩
  by_contra hX
  exact hj κ (by simp only [Ky, Finset.mem_filter, Finset.mem_univ, true_and]; exact ⟨hX, ⟨j, h⟩⟩) h

theorem cover {i : Fin p} {j : Fin q} (h : col (C i j) ∈ Kl col) : i ∈ R1 col ∨ j ∈ C1 col := by
  by_cases hr : RowT col (col (C i j))
  · left
    simp only [R1, Finset.mem_biUnion]
    refine ⟨col (C i j), h, ?_⟩
    rw [rowsL, if_pos hr]
    simp only [Finset.mem_filter, Finset.mem_univ, true_and]
    exact ⟨j, rfl⟩
  · right
    simp only [C1, Finset.mem_biUnion]
    refine ⟨col (C i j), h, ?_⟩
    rw [colsL, if_neg hr]
    simp only [Finset.mem_filter, Finset.mem_univ, true_and]
    exact ⟨i, rfl⟩

end Defs

theorem nat_aux {a b z : ℕ} (h2a : 2 ≤ a) (h2b : 2 ≤ b) (hS : a * b ≤ z) (hz : z + 1 ≤ a + b) :
    False := by
  obtain ⟨a', rfl⟩ := Nat.exists_eq_add_of_le h2a
  obtain ⟨b', rfl⟩ := Nat.exists_eq_add_of_le h2b
  rw [Nat.add_mul, Nat.mul_add, Nat.mul_add] at hS
  generalize a' * b' = w at hS
  omega

section Proof

variable {K : Type*} {col : Vx p q → K} (hcol : ∀ x y, Adj x y → col x ≠ col y)
include hcol

/-- (a) a colour class meets at most one row block. -/
theorem inX_unique {κ : K} {i k : Fin p} (hi : InX col κ i) (hk : InX col κ k) : i = k := by
  by_contra hne
  rcases hi with hi | ⟨j, hi⟩ <;> rcases hk with hk | ⟨l, hk⟩
  · exact hcol (P i) (P k) hne (hi.trans hk.symm)
  · exact hcol (P i) (A k l) hne (hi.trans hk.symm)
  · exact hcol (A i j) (P k) (Ne.symm hne) (hi.trans hk.symm)
  · exact hcol (A i j) (A k l) hne (hi.trans hk.symm)

/-- (a) a colour class meets at most one column block. -/
theorem inY_unique {κ : K} {j l : Fin q} (hj : InY col κ j) (hl : InY col κ l) : j = l := by
  by_contra hne
  rcases hj with hj | ⟨i, hj⟩ <;> rcases hl with hl | ⟨k, hl⟩
  · exact hcol (Q j) (Q l) hne (hj.trans hl.symm)
  · exact hcol (Q j) (B k l) hne (hj.trans hl.symm)
  · exact hcol (B i j) (Q l) (Ne.symm hne) (hj.trans hl.symm)
  · exact hcol (B i j) (B k l) hne (hj.trans hl.symm)

/-- (b) `C_ij` is adjacent to all of `X_k`, `k ≠ i`. -/
theorem inX_of_C {κ : K} {i k : Fin p} {j : Fin q} (hc : col (C i j) = κ) (hk : InX col κ k) :
    k = i := by
  by_contra hne
  rcases hk with hk | ⟨l, hk⟩
  · exact hcol (P k) (C i j) hne (hk.trans hc.symm)
  · exact hcol (A k l) (C i j) hne (hk.trans hc.symm)

/-- (b) `C_ij` is adjacent to all of `Y_l`, `l ≠ j`. -/
theorem inY_of_C {κ : K} {i : Fin p} {j l : Fin q} (hc : col (C i j) = κ) (hl : InY col κ l) :
    l = j := by
  by_contra hne
  rcases hl with hl | ⟨k, hl⟩
  · exact hcol (Q l) (C i j) hne (hl.trans hc.symm)
  · exact hcol (B k l) (C i j) hne (hl.trans hc.symm)

/-- (c) the `C`-vertices of a class lie in one row or in one column. -/
theorem C_line (κ : K) :
    RowT col κ ∨ ∀ i j i' j', col (C i j) = κ → col (C i' j') = κ → j = j' := by
  have key : ∀ i j i' j', col (C i j) = κ → col (C i' j') = κ → i = i' ∨ j = j' := by
    intro i j i' j' h1 h2
    by_contra h
    push_neg at h
    exact hcol (C i j) (C i' j') h (h1.trans h2.symm)
  by_cases hr : RowT col κ
  · exact Or.inl hr
  · right
    simp only [RowT] at hr
    push_neg at hr
    obtain ⟨i1, j1, i2, j2, h1, h2, hi⟩ := hr
    have e12 : j1 = j2 := (key _ _ _ _ h1 h2).resolve_left hi
    subst e12
    have hcolj : ∀ i j, col (C i j) = κ → j = j1 := by
      intro i j h
      by_contra hj
      have a1 : i = i1 := (key _ _ _ _ h h1).resolve_right hj
      have a2 : i = i2 := (key _ _ _ _ h h2).resolve_right hj
      exact hi (a1.symm.trans a2)
    intro i j i' j' h h'
    rw [hcolj i j h, hcolj i' j' h']

variable [Fintype K]

theorem p_le : p ≤ (Kx col).card + (R0 col).card := by
  classical
  have h1 : (Finset.univ : Finset (Fin p)).card ≤ (Finset.univ \ R0 col).card + (R0 col).card :=
    Finset.card_le_card_sdiff_add_card
  have h2 : (Finset.univ \ R0 col).card ≤ (Kx col).card := by
    apply Finset.card_le_card_of_forall_subsingleton (fun i κ => InX col κ i)
    · intro i hi
      simp only [R0, Finset.mem_sdiff, Finset.mem_univ, true_and, Finset.mem_filter] at hi
      push_neg at hi
      exact hi
    · intro κ _ i hi k hk
      obtain ⟨_, hi⟩ := hi
      obtain ⟨_, hk⟩ := hk
      exact inX_unique hcol hi hk
  rw [Finset.card_univ, Fintype.card_fin] at h1
  omega

theorem q_le : q ≤ (Ky col).card + (C0 col).card := by
  classical
  have h1 : (Finset.univ : Finset (Fin q)).card ≤ (Finset.univ \ C0 col).card + (C0 col).card :=
    Finset.card_le_card_sdiff_add_card
  have h2 : (Finset.univ \ C0 col).card ≤ (Ky col).card := by
    apply Finset.card_le_card_of_forall_subsingleton (fun j κ => InY col κ j)
    · intro j hj
      simp only [C0, Finset.mem_sdiff, Finset.mem_univ, true_and, Finset.mem_filter] at hj
      push_neg at hj
      exact hj
    · intro κ _ j hj l hl
      obtain ⟨_, hj⟩ := hj
      obtain ⟨_, hl⟩ := hl
      exact inY_unique hcol hj hl
  rw [Finset.card_univ, Fintype.card_fin] at h1
  omega

theorem R0_le : (R0 col).card ≤ (Kxy col).card := by
  apply Finset.card_le_card_of_injOn (fun i => col (P i))
  · intro i hi
    exact Finset.mem_coe.mpr (mem_Kxy_of_R0 col (Finset.mem_coe.mp hi) (Or.inl rfl))
  · intro i _ k _ h
    exact inX_unique hcol (Or.inl rfl) (Or.inl h.symm)

theorem C0_le : (C0 col).card ≤ (Kxy col).card := by
  apply Finset.card_le_card_of_injOn (fun j => col (Q j))
  · intro j hj
    exact Finset.mem_coe.mpr (mem_Kxy_of_C0 col (Finset.mem_coe.mp hj) (Or.inl rfl))
  · intro j _ l _ h
    exact inY_unique hcol (Or.inl rfl) (Or.inl h.symm)

/-- For `i ∈ R₀`, `j ∈ C₀` the vertex `C_ij` lies in an `XY` class or an `L` class (NOTE, by (b)). -/
theorem cell {i : Fin p} {j : Fin q} (hi : i ∈ R0 col) (hj : j ∈ C0 col) :
    col (C i j) ∈ Kxy col ∨ col (C i j) ∈ Kl col := by
  classical
  by_cases hx : ∃ k, InX col (col (C i j)) k
  · obtain ⟨k, hk⟩ := hx
    have hki : k = i := inX_of_C hcol (rfl : col (C i j) = _) hk
    rw [hki] at hk
    exact Or.inl (mem_Kxy_of_R0 col hi hk)
  · by_cases hy : ∃ l, InY col (col (C i j)) l
    · obtain ⟨l, hl⟩ := hy
      have hlj : l = j := inY_of_C hcol (rfl : col (C i j) = _) hl
      rw [hlj] at hl
      exact Or.inl (mem_Kxy_of_C0 col hj hl)
    · right
      simp only [Kl, Finset.mem_filter, Finset.mem_univ, true_and]
      exact ⟨hx, hy⟩

theorem R1C1_le : (R1 col).card + (C1 col).card ≤ (Kl col).card := by
  classical
  have h1 : ∀ κ, (rowsL col κ).card + (colsL col κ).card ≤ 1 := by
    intro κ
    by_cases hr : RowT col κ
    · rw [rowsL, colsL, if_pos hr, if_pos hr, Finset.card_empty, add_zero, Finset.card_le_one]
      intro a ha b hb
      simp only [Finset.mem_filter, Finset.mem_univ, true_and] at ha hb
      obtain ⟨j, hj⟩ := ha
      obtain ⟨j', hj'⟩ := hb
      exact hr a j b j' hj hj'
    · rw [rowsL, colsL, if_neg hr, if_neg hr, Finset.card_empty, zero_add, Finset.card_le_one]
      intro a ha b hb
      simp only [Finset.mem_filter, Finset.mem_univ, true_and] at ha hb
      obtain ⟨i, hi⟩ := ha
      obtain ⟨i', hi'⟩ := hb
      exact ((C_line hcol κ).resolve_left hr) i a i' b hi hi'
  have e1 : (R1 col).card ≤ ∑ κ ∈ Kl col, (rowsL col κ).card := by
    unfold R1; exact Finset.card_biUnion_le
  have e2 : (C1 col).card ≤ ∑ κ ∈ Kl col, (colsL col κ).card := by
    unfold C1; exact Finset.card_biUnion_le
  have e3 : ∑ κ ∈ Kl col, ((rowsL col κ).card + (colsL col κ).card) ≤ ∑ _κ ∈ Kl col, 1 :=
    Finset.sum_le_sum fun κ _ => h1 κ
  rw [Finset.sum_add_distrib, ← Finset.card_eq_sum_ones] at e3
  omega

theorem mem_S_Kxy {i : Fin p} {j : Fin q} (hi : i ∈ R0 col \ R1 col) (hj : j ∈ C0 col \ C1 col) :
    col (C i j) ∈ Kxy col := by
  rw [Finset.mem_sdiff] at hi hj
  rcases cell hcol hi.1 hj.1 with h | h
  · exact h
  · rcases cover col h with h' | h'
    · exact absurd h' hi.2
    · exact absurd h' hj.2

/-- The cells of `(R₀ \ R₁) × (C₀ \ C₁)` need distinct `XY` classes. -/
theorem S_le : ((R0 col \ R1 col) ×ˢ (C0 col \ C1 col)).card ≤ (Kxy col).card := by
  apply Finset.card_le_card_of_injOn (fun c => col (C c.1 c.2))
  · intro c hc
    have hc' := Finset.mem_product.mp (Finset.mem_coe.mp hc)
    exact Finset.mem_coe.mpr (mem_S_Kxy hcol hc'.1 hc'.2)
  · intro c hc c' _ h
    have hc' := Finset.mem_product.mp (Finset.mem_coe.mp hc)
    have hm := mem_S_Kxy hcol hc'.1 hc'.2
    simp only [Kxy, Finset.mem_filter, Finset.mem_univ, true_and] at hm
    obtain ⟨⟨k, hk⟩, ⟨l, hl⟩⟩ := hm
    have h' : col (C c'.1 c'.2) = col (C c.1 c.2) := h.symm
    have e1 : k = c.1 := inX_of_C hcol (rfl : col (C c.1 c.2) = _) hk
    have e2 : k = c'.1 := inX_of_C hcol h' hk
    have e3 : l = c.2 := inY_of_C hcol (rfl : col (C c.1 c.2) = _) hl
    have e4 : l = c'.2 := inY_of_C hcol h' hl
    exact Prod.ext (e1.symm.trans e2) (e3.symm.trans e4)

/-- Case `|R₀ \ R₁| = 1` (NOTE: `r₀ = 1`): `P_{i*}` and `B_{i*j}` would share a class. -/
theorem case_a (h1 : (R0 col \ R1 col).card = 1)
    (hz : (Kxy col).card ≤ (C0 col \ C1 col).card) : False := by
  classical
  obtain ⟨i0, hi0⟩ := Finset.card_eq_one.mp h1
  have hi0R : i0 ∈ R0 col \ R1 col := by rw [hi0]; exact Finset.mem_singleton_self _
  have hi0R0 : i0 ∈ R0 col := (Finset.mem_sdiff.mp hi0R).1
  have hTsub : (C0 col \ C1 col).image (fun j => col (C i0 j)) ⊆ Kxy col := by
    intro κ hκ
    rw [Finset.mem_image] at hκ
    obtain ⟨j, hj, rfl⟩ := hκ
    exact mem_S_Kxy hcol hi0R hj
  have hTcard : ((C0 col \ C1 col).image (fun j => col (C i0 j))).card = (C0 col \ C1 col).card := by
    apply Finset.card_image_of_injOn
    intro j hj j' _ h
    have hm := mem_S_Kxy hcol hi0R (Finset.mem_coe.mp hj)
    simp only [Kxy, Finset.mem_filter, Finset.mem_univ, true_and] at hm
    obtain ⟨_, ⟨l, hl⟩⟩ := hm
    have h' : col (C i0 j') = col (C i0 j) := h.symm
    exact (inY_of_C hcol (rfl : col (C i0 j) = _) hl).symm.trans (inY_of_C hcol h' hl)
  have hTeq := Finset.eq_of_subset_of_card_le hTsub (by omega)
  have hP : col (P i0) ∈ (C0 col \ C1 col).image (fun j => col (C i0 j)) := by
    rw [hTeq]; exact mem_Kxy_of_R0 col hi0R0 (Or.inl rfl)
  rw [Finset.mem_image] at hP
  obtain ⟨j0, hj0, hPj0⟩ := hP
  have hj0C0 : j0 ∈ C0 col := (Finset.mem_sdiff.mp hj0).1
  have hB : col (B i0 j0) ∈ (C0 col \ C1 col).image (fun j => col (C i0 j)) := by
    rw [hTeq]; exact mem_Kxy_of_C0 col hj0C0 (Or.inr ⟨i0, rfl⟩)
  rw [Finset.mem_image] at hB
  obtain ⟨j1, _, hBj1⟩ := hB
  have e : j0 = j1 := inY_of_C hcol hBj1 (Or.inr ⟨i0, rfl⟩)
  subst e
  exact hcol (P i0) (B i0 j0) rfl (hPj0.symm.trans hBj1)

/-- Case `|C₀ \ C₁| = 1` (NOTE: `c₀ = 1`): `Q_{j*}` and `A_{ij*}` would share a class. -/
theorem case_b (h1 : (C0 col \ C1 col).card = 1)
    (hz : (Kxy col).card ≤ (R0 col \ R1 col).card) : False := by
  classical
  obtain ⟨j0, hj0⟩ := Finset.card_eq_one.mp h1
  have hj0C : j0 ∈ C0 col \ C1 col := by rw [hj0]; exact Finset.mem_singleton_self _
  have hj0C0 : j0 ∈ C0 col := (Finset.mem_sdiff.mp hj0C).1
  have hTsub : (R0 col \ R1 col).image (fun i => col (C i j0)) ⊆ Kxy col := by
    intro κ hκ
    rw [Finset.mem_image] at hκ
    obtain ⟨i, hi, rfl⟩ := hκ
    exact mem_S_Kxy hcol hi hj0C
  have hTcard : ((R0 col \ R1 col).image (fun i => col (C i j0))).card = (R0 col \ R1 col).card := by
    apply Finset.card_image_of_injOn
    intro i hi i' _ h
    have hm := mem_S_Kxy hcol (Finset.mem_coe.mp hi) hj0C
    simp only [Kxy, Finset.mem_filter, Finset.mem_univ, true_and] at hm
    obtain ⟨⟨k, hk⟩, _⟩ := hm
    have h' : col (C i' j0) = col (C i j0) := h.symm
    exact (inX_of_C hcol (rfl : col (C i j0) = _) hk).symm.trans (inX_of_C hcol h' hk)
  have hTeq := Finset.eq_of_subset_of_card_le hTsub (by omega)
  have hQ : col (Q j0) ∈ (R0 col \ R1 col).image (fun i => col (C i j0)) := by
    rw [hTeq]; exact mem_Kxy_of_C0 col hj0C0 (Or.inl rfl)
  rw [Finset.mem_image] at hQ
  obtain ⟨i0, hi0, hQi0⟩ := hQ
  have hi0R0 : i0 ∈ R0 col := (Finset.mem_sdiff.mp hi0).1
  have hA : col (A i0 j0) ∈ (R0 col \ R1 col).image (fun i => col (C i j0)) := by
    rw [hTeq]; exact mem_Kxy_of_R0 col hi0R0 (Or.inr ⟨j0, rfl⟩)
  rw [Finset.mem_image] at hA
  obtain ⟨i1, _, hAi1⟩ := hA
  have e : i0 = i1 := inX_of_C hcol hAi1 (Or.inr ⟨j0, rfl⟩)
  subst e
  exact hcol (Q j0) (A i0 j0) rfl (hQi0.symm.trans hAi1)

/-- **Colouring lemma** (NOTE §3): a proper colouring of `H(p, q)` uses at least `p + q` colours. -/
theorem card_le : p + q ≤ Fintype.card K := by
  by_contra hlt
  push_neg at hlt
  have hs := card_split col
  have hp := p_le hcol
  have hq := q_le hcol
  have hr := R0_le hcol
  have hc := C0_le hcol
  have hl := R1C1_le hcol
  have ha := Finset.card_le_card_sdiff_add_card (s := R0 col) (t := R1 col)
  have hb := Finset.card_le_card_sdiff_add_card (s := C0 col) (t := C1 col)
  have hz : (Kxy col).card + 1 ≤ (R0 col \ R1 col).card + (C0 col \ C1 col).card := by omega
  have ha1 : 1 ≤ (R0 col \ R1 col).card := by omega
  have hb1 : 1 ≤ (C0 col \ C1 col).card := by omega
  have hab : (R0 col \ R1 col).card = 1 ∨ (C0 col \ C1 col).card = 1 := by
    by_contra h
    push_neg at h
    have h2a : 2 ≤ (R0 col \ R1 col).card := by omega
    have h2b : 2 ≤ (C0 col \ C1 col).card := by omega
    have hS := S_le hcol
    rw [Finset.card_product] at hS
    exact nat_aux h2a h2b hS hz
  rcases hab with h1 | h1
  · exact case_a hcol h1 (by omega)
  · exact case_b hcol h1 (by omega)

end Proof

/-- NOTE §3 Lemma, any colour set: fewer than `p + q` colours force a monochromatic edge. -/
theorem not_colourable {K : Type*} [Fintype K] (hK : Fintype.card K < p + q) (col : Vx p q → K) :
    ∃ x y, Adj x y ∧ col x = col y := by
  by_contra h
  push_neg at h
  have := card_le h
  omega

/-- **NOTE §3 Lemma**: `H(p, q)` has no proper colouring with `p + q - 1` colours. -/
theorem lemma_H (hpq : 1 ≤ p + q) (col : Vx p q → Fin (p + q - 1)) :
    ∃ x y, Adj x y ∧ col x = col y :=
  not_colourable (by rw [Fintype.card_fin]; omega) col

end Vx

end Gale7
