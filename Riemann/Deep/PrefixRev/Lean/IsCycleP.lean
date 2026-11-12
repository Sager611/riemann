import Riemann.Deep.PrefixRev.Lean.CycleC
import Mathlib.GroupTheory.Perm.Closure

/-!
# P is an n-cycle in Mathlib's sense; a criterion for G = S_n

On the line with gcd(m+1, k+1) = 1: P.IsCycle and P.support = univ (Theorem B.1, n-cycle).
Then G = S_n as soon as G contains the transposition (x, P x) for one x
(Mathlib: Equiv.Perm.closure_cycle_adjacent_swap).
-/

namespace PrefixRev

open Equiv

theorem P_mem (n m k : ℕ) : P n m k ∈ G n m k :=
  (G n m k).mul_mem ((G n m k).mul_mem (rk_mem n m k) (rn_mem n m k)) (rm_mem n m k)

theorem P_ne {n m k : ℕ} (hline : n = m + k + 1) (hkm : k < m) (x : Fin n) : P n m k x ≠ x := by
  intro h
  have := congrArg Fin.val h
  rw [P_val hline hkm x] at this
  split_ifs at this <;> omega

theorem P_support {n m k : ℕ} (hline : n = m + k + 1) (hkm : k < m) :
    (P n m k).support = Finset.univ :=
  Finset.eq_univ_iff_forall.mpr fun x => Perm.mem_support.mpr (P_ne hline hkm x)

/-- Theorem B.1 (n-cycle, Mathlib form). -/
theorem P_isCycle {n m k : ℕ} (hline : n = m + k + 1) (hkm : k < m)
    (hg : Nat.Coprime (k + 1) (m + 1)) : (P n m k).IsCycle := by
  refine ⟨⟨0, by omega⟩, P_ne hline hkm _, fun y _ => ?_⟩
  obtain ⟨t, ht⟩ := P_cycle hline hkm hg ⟨0, by omega⟩ y
  exact ⟨t, by rw [zpow_natCast]; exact ht⟩

/-- Criterion: if G contains the transposition (x, P x) for some x, then G = S_n. -/
theorem G_eq_top_of_swap {n m k : ℕ} (hline : n = m + k + 1) (hkm : k < m)
    (hg : Nat.Coprime (k + 1) (m + 1)) (x : Fin n) (hs : swap x (P n m k x) ∈ G n m k) :
    G n m k = ⊤ := by
  rw [eq_top_iff, ← Perm.closure_cycle_adjacent_swap (P_isCycle hline hkm hg) (P_support hline hkm) x,
    Subgroup.closure_le]
  intro g hg'
  simp only [Set.mem_insert_iff, Set.mem_singleton_iff] at hg'
  rcases hg' with rfl | rfl
  · exact P_mem n m k
  · exact hs

end PrefixRev
