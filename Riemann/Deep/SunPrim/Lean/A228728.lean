import Riemann.Deep.SunPrim.Lean.ListAux

/-!
# Sun, OEIS A228728, comment (real numbers)

Stated for an arbitrary linearly ordered field `R` (`R = ℝ` in the source) and refuted for every such
field, hence for `ℝ`. The finite check runs over integer lists (`decide`).
-/

namespace SunLean

/-- Z.-W. Sun, OEIS A228728, comment ("Conjecture: …"), over an ordered field `R` (the source:
`R = ℝ`): for `a_1 < … < a_n`, if some permutation `b` of the `a_i` has the `n` cyclic distances
`|b_1 - b_2|, …, |b_n - b_1|` pairwise distinct, then some permutation `c` with `c_1 = a_1` and
`c_n = a_n` has the `n` cyclic distances `|c_1 - c_2|, …, |c_n - c_1|` pairwise distinct. -/
def Conj_A228728 (R : Type) [Field R] [LinearOrder R] : Prop :=
  ∀ a : List R, a.Pairwise (· < ·) →
    (∃ b : List R, b.Perm a ∧ (List.zipWith (fun x y => |x - y|) b (b.rotate 1)).Nodup) →
    ∃ c : List R, c.Perm a ∧ c.head? = a.head? ∧ c.getLast? = a.getLast? ∧
      (List.zipWith (fun x y => |x - y|) c (c.rotate 1)).Nodup

/-- `a = (0, 1, 2, 3, 4, 5, 8)`, `b = (0, 3, 4, 2, 8, 1, 5)` (distances `3, 1, 2, 6, 7, 4, 5`); no `c`
with ends `0, 8` works. In every linearly ordered field, e.g. `ℝ`. -/
theorem not_conj_A228728 (R : Type) [Field R] [LinearOrder R] [IsStrictOrderedRing R] :
    ¬ Conj_A228728 R := by
  intro h
  have hg : ∀ x y : ℤ, |(x : R) - (y : R)| = ((|x - y| : ℤ) : R) := fun x y => by
    rw [Int.cast_abs, Int.cast_sub]
  obtain ⟨c, hc, hh, hl, hd⟩ := h (([0, 1, 2, 3, 4, 5, 8] : List ℤ).map (Int.cast : ℤ → R))
    ((List.pairwise_map).2 ((by decide : ([0, 1, 2, 3, 4, 5, 8] : List ℤ).Pairwise (· < ·)).imp
      (fun hxy => Int.cast_lt.2 hxy)))
    ⟨([0, 3, 4, 2, 8, 1, 5] : List ℤ).map (Int.cast : ℤ → R),
      List.Perm.map _ (by decide), by
        rw [zipWith_rotate_map _ (fun x y : R => |x - y|) (fun x y : ℤ => |x - y|) hg,
          List.nodup_map_iff Int.cast_injective]
        decide⟩
  obtain ⟨c₀, h₀, rfl⟩ := exists_of_perm_map _ hc
  rw [List.head?_map, List.head?_map] at hh
  rw [List.getLast?_map, List.getLast?_map] at hl
  have hh₀ : c₀.head? = some 0 := Option.map_injective Int.cast_injective hh
  have hl₀ : c₀.getLast? = some 8 := Option.map_injective Int.cast_injective hl
  rw [zipWith_rotate_map _ (fun x y : R => |x - y|) (fun x y : ℤ => |x - y|) hg,
    List.nodup_map_iff Int.cast_injective] at hd
  have hc₀ : c₀.Perm [0, 1, 2, 3, 4, 5, 8] := List.mem_permutations'.1 h₀
  -- `c₀ = 0 :: (m ++ [8])` with `m` a permutation of `[1, 2, 3, 4, 5]`
  obtain ⟨m, rfl⟩ : ∃ m, c₀ = 0 :: (m ++ [8]) := by
    rcases c₀ with _ | ⟨x, t⟩
    · simp at hh₀
    · simp only [List.head?_cons, Option.some.injEq] at hh₀
      subst hh₀
      rcases List.eq_nil_or_concat t with rfl | ⟨L, b, rfl⟩
      · simp at hl₀
      · simp only [List.concat_eq_append] at hl₀ ⊢
        rw [← List.cons_append, List.getLast?_concat, Option.some.injEq] at hl₀
        subst hl₀
        exact ⟨L, rfl⟩
  have hm8 : (m ++ [8]).Perm ([1, 2, 3, 4, 5] ++ [8]) := (List.perm_cons 0).1 hc₀
  have hm : m.Perm [1, 2, 3, 4, 5] := (List.perm_append_right_iff [8]).1 hm8
  have key : ∀ m ∈ ([1, 2, 3, 4, 5] : List ℤ).permutations',
      ¬ (List.zipWith (fun x y => |x - y|) (0 :: (m ++ [8])) ((0 :: (m ++ [8])).rotate 1)).Nodup := by
    decide +kernel
  exact key m (List.mem_permutations'.2 hm) hd

end SunLean

#print axioms SunLean.not_conj_A228728
