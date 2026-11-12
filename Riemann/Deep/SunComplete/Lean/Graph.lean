import Mathlib.Combinatorics.SimpleGraph.Basic
import Mathlib.Logic.Equiv.Fin.Rotate
import Mathlib.Data.Fintype.Card

/-!
# sunlean · Graph: Hamiltonian cycles, independent sets, vertex cuts, Chvátal–Erdős (statement)

`HasHamCycle G`: a cyclic numbering of all vertices with consecutive vertices adjacent.
`ChvatalErdos`: the Chvátal–Erdős theorem (κ ≥ α, n ≥ 3 ⇒ Hamiltonian), a named hypothesis.
-/

namespace SunLean

variable {V : Type} [Fintype V]

/-- A Hamiltonian cycle of `G`: a numbering `a 0, …, a (n-1)` of all `n` vertices (a bijection
`Fin n ≃ V`) with `a i` adjacent to `a (i + 1 mod n)` (`finRotate`) for every `i`. For `n ≤ 2` this
degenerates (it holds for the empty graph and for `K₂`); it is used only for `n ≥ 3`. -/
def HasHamCycle (G : SimpleGraph V) : Prop :=
  ∃ a : Fin (Fintype.card V) ≃ V, ∀ i, G.Adj (a i) (a (finRotate _ i))

/-- `I` is an independent set of `G`. -/
def IndepSet (G : SimpleGraph V) (I : Finset V) : Prop :=
  ∀ x ∈ I, ∀ y ∈ I, x ≠ y → ¬ G.Adj x y

/-- `G − S` (delete the vertices of `S`) is nonempty and connected: every split of `V ∖ S` into
two nonempty parts `A`, `(V ∖ S) ∖ A` has an edge across. -/
def ConnDel [DecidableEq V] (G : SimpleGraph V) (S : Finset V) : Prop :=
  (Finset.univ \ S).Nonempty ∧
    ∀ A : Finset V, A ⊆ Finset.univ \ S → A.Nonempty → ((Finset.univ \ S) \ A).Nonempty →
      ∃ x ∈ A, ∃ y ∈ (Finset.univ \ S) \ A, G.Adj x y

/-- **Chvátal–Erdős** (V. Chvátal, P. Erdős, Discrete Math. 2 (1972) 111–113): a graph on `n ≥ 3`
vertices whose connectivity `κ` is at least its independence number `α` is Hamiltonian.
Stated with an integer `k`: `α ≤ k` and `G − S` connected for every `|S| < k`. For `k < n` this
says `κ ≥ k ≥ α`; `k = n` forces `G` complete (Hamiltonian as `n ≥ 3`); for `k > n` the deletion
premise fails at `S = univ`, so that case is vacuous. -/
def ChvatalErdos : Prop :=
  ∀ (V : Type) [Fintype V] [DecidableEq V] (G : SimpleGraph V) (k : ℕ),
    3 ≤ Fintype.card V →
    (∀ I : Finset V, IndepSet G I → I.card ≤ k) →
    (∀ S : Finset V, S.card < k → ConnDel G S) →
    HasHamCycle G

/-- For `n ≥ 2`, `finRotate n` has no fixed point. -/
theorem finRotate_ne : ∀ {n : ℕ}, 2 ≤ n → ∀ i : Fin n, finRotate n i ≠ i
  | m + 2, _, i => by
    rw [finRotate_apply]
    intro h
    have h1 : (1 : Fin (m + 2)) = 0 := by
      have := congrArg (fun x => x - i) h
      simp at this
    have h2 := congrArg Fin.val h1
    simp at h2

/-- For `n ≥ 3`: `rot i ≠ rot⁻¹ i`. -/
theorem finRotate_ne_symm : ∀ {n : ℕ}, 3 ≤ n → ∀ i : Fin n, finRotate n i ≠ (finRotate n).symm i
  | m + 3, _, i => by
    intro h
    have h' := congrArg (finRotate (m + 3)) h
    rw [Equiv.apply_symm_apply, finRotate_apply, finRotate_apply] at h'
    have h1 : (1 + 1 : Fin (m + 3)) = 0 := by
      have := congrArg (fun x => x - i) h'
      simp only [add_assoc, add_sub_cancel_left, sub_self] at this
      exact this
    have h2 : ((1 + 1 : Fin (m + 3)) : ℕ) = 2 := by
      rw [Fin.val_add, Fin.val_one]
      exact Nat.mod_eq_of_lt (by omega)
    rw [h1, Fin.val_zero] at h2
    omega

end SunLean
