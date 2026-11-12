import Mathlib.Data.ZMod.Defs
import Mathlib.Data.Set.Card
import Mathlib.Algebra.BigOperators.Fin
import Mathlib.Data.ENat.Lattice

/-!
# Wang–Sun, arXiv 2210.12044 (Electron. J. Combin. 33(2) (2026) #P2.64), Conjecture 1.1, inequality (C)

Stated for every additive (abelian) group `G` with `|G| > 1` and refuted in `G = ℤ/7`, `n = 4`.
Sources and fidelity: `../LEAN.md`.
-/

namespace SunLean

/-- `p(G)`: the minimum of the orders of the nonzero elements of `G` (`+∞` if every nonzero element has
infinite order); the order of `g` is the least `m > 0` with `m • g = 0`. -/
noncomputable def pG (G : Type) [AddCommGroup G] : ℕ∞ :=
  ⨅ (g : G) (_ : g ≠ 0) (m : ℕ) (_ : 0 < m) (_ : m • g = 0), (m : ℕ∞)

/-- `L(A_1, …, A_n) = {a_1 + ⋯ + a_n : a_i ∈ A_i, a_i ≠ a_(i+1) for 0 < i < n}` (indices from 0 here). -/
def Lset {G : Type} [AddCommGroup G] (n : ℕ) (A : Fin n → Finset G) : Set G :=
  {s | ∃ a : Fin n → G, (∀ i, a i ∈ A i) ∧ (∀ i j : Fin n, j.val = i.val + 1 → a i ≠ a j) ∧
    s = ∑ i, a i}

/-- `C(A_1, …, A_n)`: as `L`, and also `a_n ≠ a_1`. -/
def Cset {G : Type} [AddCommGroup G] (n : ℕ) (A : Fin n → Finset G) : Set G :=
  {s | ∃ a : Fin n → G, (∀ i, a i ∈ A i) ∧ (∀ i j : Fin n, j.val = i.val + 1 → a i ≠ a j) ∧
    (∀ i j : Fin n, i.val + 1 = n → j.val = 0 → a i ≠ a j) ∧ s = ∑ i, a i}

/-- H. Wang and Z.-W. Sun, arXiv 2210.12044 v1, Conjecture 1.1 (TeX l. 174-184): `G` an additive group
with `|G| > 1`, `A_1, …, A_n` (`n > 1`) finite subsets with `|A_i| > 1`; then (L)
`|L(A_1, …, A_n)| ≥ min{p(G), |A_1| + ⋯ + |A_n| - 2n + 1 + {n}_2}` and (C)
`|C(A_1, …, A_n)| ≥ min{p(G), |A_1| + ⋯ + |A_n| - 2n + (-1)^n (1 + {n}_2)}`, `{n}_2 = n mod 2`.
(A negative integer bound is read as `0`: `toNat`; both sides are then trivially comparable.) -/
def Conj_WS : Prop :=
  ∀ (G : Type) [AddCommGroup G] [Nontrivial G] (n : ℕ), 1 < n → ∀ A : Fin n → Finset G,
    (∀ i, 1 < (A i).card) →
      min (pG G) (((∑ i, ((A i).card : ℤ)) - 2 * n + 1 + (n % 2 : ℕ)).toNat : ℕ∞) ≤ (Lset n A).encard ∧
      min (pG G) (((∑ i, ((A i).card : ℤ)) - 2 * n + (-1) ^ n * (1 + (n % 2 : ℕ))).toNat : ℕ∞) ≤
        (Cset n A).encard

/-- Inequality (C) alone. -/
def Conj_WS_C : Prop :=
  ∀ (G : Type) [AddCommGroup G] [Nontrivial G] (n : ℕ), 1 < n → ∀ A : Fin n → Finset G,
    (∀ i, 1 < (A i).card) →
      min (pG G) (((∑ i, ((A i).card : ℤ)) - 2 * n + (-1) ^ n * (1 + (n % 2 : ℕ))).toNat : ℕ∞) ≤
        (Cset n A).encard

theorem pG_zmod7 : (7 : ℕ∞) ≤ pG (ZMod 7) := by
  unfold pG
  refine le_iInf fun g => le_iInf fun hg => le_iInf fun m => le_iInf fun hm => le_iInf fun hmg => ?_
  by_contra hlt
  have hm7 : m < 7 := by exact_mod_cast lt_of_not_ge hlt
  have key : ∀ g : ZMod 7, g ≠ 0 → ∀ m < 7, 0 < m → m • g ≠ 0 := by decide
  exact key g hg m hm7 hm hmg

instance nontrivial_zmod7 : Nontrivial (ZMod 7) := ⟨⟨0, 1, by decide⟩⟩

/-- The counterexample: `(A_1, A_2, A_3, A_4) = ({0, 1}, {0, 1, 2, 4, 5}, {2, 4, 5}, {0, 1, 4, 5})`. -/
def WSA : Fin 4 → Finset (ZMod 7) := ![{0, 1}, {0, 1, 2, 4, 5}, {2, 4, 5}, {0, 1, 4, 5}]

/-- `G = ℤ/7`, `n = 4`, `A = WSA`: `S = 14`, bound `min{7, 14 - 8 + 1} = 7`, but `2 ∉ C`, so `|C| ≤ 6`. -/
theorem not_conj_WS_C : ¬ Conj_WS_C := by
  intro h
  have hC := h (ZMod 7) 4 (by decide) WSA (by decide)
  have hsub : Cset 4 WSA ⊆ (({0, 1, 3, 4, 5, 6} : Finset (ZMod 7)) : Set (ZMod 7)) := by
    rintro s ⟨a, hA, hadj, hcyc, rfl⟩
    rw [Fin.sum_univ_four]
    have key : ∀ x ∈ WSA 0, ∀ y ∈ WSA 1, ∀ z ∈ WSA 2, ∀ w ∈ WSA 3, x ≠ y → y ≠ z → z ≠ w →
        w ≠ x → x + y + z + w ∈ ({0, 1, 3, 4, 5, 6} : Finset (ZMod 7)) := by decide
    exact key _ (hA 0) _ (hA 1) _ (hA 2) _ (hA 3) (hadj 0 1 rfl) (hadj 1 2 rfl) (hadj 2 3 rfl)
      (hcyc 3 0 rfl rfl)
  have h6 := ((le_min pG_zmod7 (by decide)).trans hC).trans
    ((Set.encard_le_encard hsub).trans_eq (Set.encard_coe_eq_coe_finsetCard _))
  revert h6
  decide

theorem not_conj_WS : ¬ Conj_WS := fun h =>
  not_conj_WS_C fun G _ _ n hn A hA => (h G n hn A hA).2

end SunLean

#print axioms SunLean.not_conj_WS_C
#print axioms SunLean.not_conj_WS
