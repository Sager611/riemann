import Riemann.Deep.SunPrim.Lean.ListAux

/-!
# Sun, arXiv 1309.1679, Conjecture 3.10(ii) (real numbers)

Stated for an arbitrary field `R` (`R = ℝ` in the source; `Mathlib.Data.Real.Basic` is not imported,
see `../LEAN.md`) and refuted for every field of characteristic zero, hence for `ℝ`. The finite check
runs over integer lists (`decide`); the cast `ℤ → R` is injective.
-/

namespace SunLean

/-- Z.-W. Sun, arXiv 1309.1679 v9, Conjecture 3.10(ii) (TeX label `Conj3.12`, l. 668-673), over a
field `R` (the source: `R = ℝ`): for any set `A` of `n > 3` distinct nonzero numbers there is a
circular permutation `(a_1, …, a_n)` of all elements of `A` whose `n` adjacent differences
`a_1 - a_2, …, a_n - a_1` are pairwise distinct and whose `n` adjacent products `a_1 a_2, …, a_n a_1`
are pairwise distinct, unless `|A| = 4` and `A = {±s, ±t}`. -/
def Conj_3_10_ii (R : Type) [Field R] : Prop :=
  ∀ A : Finset R, 3 < A.card → (∀ x ∈ A, x ≠ 0) →
    ¬ (A.card = 4 ∧ ∃ s t : R, (A : Set R) = {s, -s, t, -t}) →
    ∃ l : List R, l.Nodup ∧ (∀ x, x ∈ l ↔ x ∈ A) ∧
      (List.zipWith (· - ·) l (l.rotate 1)).Nodup ∧ (List.zipWith (· * ·) l (l.rotate 1)).Nodup

/-- `A = {2, 3, 4, 6}` (`2 · 6 = 3 · 4`), in every field of characteristic zero, e.g. `ℝ`. -/
theorem not_conj_3_10_ii (R : Type) [Field R] [CharZero R] : ¬ Conj_3_10_ii R := by
  classical
  intro h
  have hLn : (([2, 3, 4, 6] : List ℤ).map (Int.cast : ℤ → R)).Nodup :=
    (List.nodup_map_iff Int.cast_injective).2 (by decide)
  have hcard : 3 < (([2, 3, 4, 6] : List ℤ).map (Int.cast : ℤ → R)).toFinset.card := by
    rw [List.toFinset_card_of_nodup hLn]; simp
  have hnz : ∀ x ∈ (([2, 3, 4, 6] : List ℤ).map (Int.cast : ℤ → R)).toFinset, x ≠ 0 := by
    intro x hx
    obtain ⟨k, hk, rfl⟩ := List.mem_map.1 (List.mem_toFinset.1 hx)
    have key0 : ∀ k ∈ ([2, 3, 4, 6] : List ℤ), k ≠ 0 := by decide
    exact_mod_cast key0 k hk
  have hexc : ¬ ((([2, 3, 4, 6] : List ℤ).map (Int.cast : ℤ → R)).toFinset.card = 4 ∧ ∃ s t : R,
      ((([2, 3, 4, 6] : List ℤ).map (Int.cast : ℤ → R)).toFinset : Set R) = {s, -s, t, -t}) := by
    rintro ⟨-, s, t, hst⟩
    have hs : s ∈ ((([2, 3, 4, 6] : List ℤ).map (Int.cast : ℤ → R)).toFinset : Set R) := by
      rw [hst]; simp
    have hs' : -s ∈ ((([2, 3, 4, 6] : List ℤ).map (Int.cast : ℤ → R)).toFinset : Set R) := by
      rw [hst]; simp
    simp only [Finset.mem_coe, List.mem_toFinset] at hs hs'
    obtain ⟨k, hk, rfl⟩ := List.mem_map.1 hs
    obtain ⟨k', hk', hk'e⟩ := List.mem_map.1 hs'
    have hkk : k' = -k := by exact_mod_cast hk'e
    have key2 : ∀ k ∈ ([2, 3, 4, 6] : List ℤ), -k ∉ ([2, 3, 4, 6] : List ℤ) := by decide
    exact key2 k hk (hkk ▸ hk')
  obtain ⟨l, hl, hmem, hd, hp⟩ :=
    h (([2, 3, 4, 6] : List ℤ).map (Int.cast : ℤ → R)).toFinset hcard hnz hexc
  have hperm : l.Perm (([2, 3, 4, 6] : List ℤ).map (Int.cast : ℤ → R)) :=
    (List.perm_ext_iff_of_nodup hl hLn).2 (fun x => (hmem x).trans List.mem_toFinset)
  obtain ⟨l₀, h₀, rfl⟩ := exists_of_perm_map _ hperm
  rw [zipWith_rotate_map _ (fun x y : R => x - y) (fun x y : ℤ => x - y)
    (fun x y => (Int.cast_sub x y).symm), List.nodup_map_iff Int.cast_injective] at hd
  rw [zipWith_rotate_map _ (fun x y : R => x * y) (fun x y : ℤ => x * y)
    (fun x y => (Int.cast_mul x y).symm), List.nodup_map_iff Int.cast_injective] at hp
  have key : ∀ l₀ ∈ ([2, 3, 4, 6] : List ℤ).permutations',
      ¬ ((List.zipWith (· - ·) l₀ (l₀.rotate 1)).Nodup ∧
        (List.zipWith (· * ·) l₀ (l₀.rotate 1)).Nodup) := by decide
  exact key l₀ h₀ ⟨hd, hp⟩

end SunLean

#print axioms SunLean.not_conj_3_10_ii
