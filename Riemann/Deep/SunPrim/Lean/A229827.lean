import Riemann.Deep.SunPrim.Lean.F8Top
import Riemann.Deep.SunPrim.Lean.F9Top
import Riemann.Deep.SunPrim.Lean.CircAux
import Mathlib.Algebra.Polynomial.Degree.Defs
import Mathlib.Algebra.Polynomial.Eval.Defs

/-!
# Sun, OEIS A229827, Conjecture (i), general part (finite fields)

Stated for every finite field `F` and every polynomial `P ∈ F[X]` (Mathlib `Polynomial`) of degree
`< q - 1` not of the form `c - x`; refuted in `F₈` with `P = x²` and in `F₉` with `P = x³`. `F₈`, `F₉` are
explicit tables on `Fin 8`, `Fin 9` whose field axioms are checked by `decide` (`F8.lean`, `F9.lean`;
every field with `q` elements is isomorphic to these). `GaloisField 2 3` is not used: its import does
not fit the 290 MB cap and its elements are not computable. A circular order is rotated to start
at 0 (`CircAux.lean`); the search below 0 is `F8Top.lean`, `F9Top.lean`.
-/

namespace SunLean

open Polynomial

/-- Z.-W. Sun, OEIS A229827, comment "Conjecture: (i) …" (text format l. 9-11): for any finite field
`F(q)` with `q` elements and a polynomial `P(x)` over `F(q)` of degree smaller than `q - 1`, if `P(x)`
is not of the form `c - x`, then there is a circular permutation `a_1, …, a_q` of all the elements of
`F(q)` with `{P(a_1) + a_2, …, P(a_(q-1)) + a_q, P(a_q) + a_1}` equal to `F(q)`. -/
def Conj_A229827_i : Prop :=
  ∀ (F : Type) [Field F] [Fintype F] (P : F[X]),
    P.degree < ((Fintype.card F - 1 : ℕ) : WithBot ℕ) → (∀ c : F, P ≠ C c - X) →
    ∃ a : List F, a.Nodup ∧ (∀ x, x ∈ a) ∧
      ∀ y : F, y ∈ List.zipWith (fun x z => P.eval x + z) a (a.rotate 1)

/-- F₈, `P(x) = x^2` (Frobenius, not of the form `c - x`): no circular order works. -/
theorem not_conj_A229827_i_F8 : ¬ Conj_A229827_i := by
  intro h
  have hcard : Fintype.card F8 = 8 := F8.card_eq
  have hdeg : (X ^ 2 : F8[X]).degree < ((Fintype.card F8 - 1 : ℕ) : WithBot ℕ) := by
    rw [degree_X_pow, hcard]; decide
  have key : ∀ c : F8, ∃ x : F8, x ^ 2 ≠ c - x := by decide +kernel
  have hnot : ∀ c : F8, (X ^ 2 : F8[X]) ≠ C c - X := by
    intro c hc
    obtain ⟨x, hx⟩ := key c
    apply hx
    simpa using congrArg (eval x) hc
  obtain ⟨a, ha, hall, hcov⟩ := h F8 (X ^ 2) hdeg hnot
  simp only [eval_X_pow] at hcov
  have huniv : a.toFinset = Finset.univ :=
    Finset.eq_univ_iff_forall.2 (fun x => List.mem_toFinset.2 (hall x))
  have hlen : a.length = 8 := by
    rw [← List.toFinset_card_of_nodup ha, huniv, Finset.card_univ, hcard]
  obtain ⟨k, hk⟩ := exists_rotate_head (hall 0)
  have hb : (a.rotate k).Nodup := List.nodup_rotate.2 ha
  have hlenb : (a.rotate k).length = 8 := by rw [List.length_rotate, hlen]
  have hcovb : ∀ y : F8, y ∈ List.zipWith (fun x z => x ^ 2 + z) (a.rotate k) ((a.rotate k).rotate 1) := by
    intro y; rw [sums_rotate]; exact List.mem_rotate.2 (hcov y)
  have hS := nodup_of_cover _ (by simp only [List.length_zipWith, List.length_rotate, hlen, min_self, hcard]) hcovb
  have hfun : (fun x z : F8 => x ^ 2 + z) = F8.sT := by
    funext x z; exact (F8.sT_eq x z).symm
  rw [hfun] at hS
  obtain ⟨x0, x1, x2, x3, x4, x5, x6, x7, hx⟩ := eq_of_length_eight _ hlenb
  rw [hx] at hb hS hk
  simp only [List.head?_cons, Option.some.injEq] at hk
  subst hk
  have g1 : x1 ∉ [0] := (List.nodup_cons.1 (nodup_rev_take hb 2)).1
  have g2 : x2 ∉ [x1, 0] := (List.nodup_cons.1 (nodup_rev_take hb 3)).1
  have g3 : x3 ∉ [x2, x1, 0] := (List.nodup_cons.1 (nodup_rev_take hb 4)).1
  have g4 : x4 ∉ [x3, x2, x1, 0] := (List.nodup_cons.1 (nodup_rev_take hb 5)).1
  have g5 : x5 ∉ [x4, x3, x2, x1, 0] := (List.nodup_cons.1 (nodup_rev_take hb 6)).1
  have g6 : x6 ∉ [x5, x4, x3, x2, x1, 0] := (List.nodup_cons.1 (nodup_rev_take hb 7)).1
  have g7 : x7 ∉ [x6, x5, x4, x3, x2, x1, 0] := (List.nodup_cons.1 (nodup_rev_take hb 8)).1
  have s2 : F8.sT x1 x2 ∉ [F8.sT 0 x1] := (List.nodup_cons.1 (nodup_rev_take hS 2)).1
  have s3 : F8.sT x2 x3 ∉ [F8.sT x1 x2, F8.sT 0 x1] := (List.nodup_cons.1 (nodup_rev_take hS 3)).1
  have s4 : F8.sT x3 x4 ∉ [F8.sT x2 x3, F8.sT x1 x2, F8.sT 0 x1] := (List.nodup_cons.1 (nodup_rev_take hS 4)).1
  have s5 : F8.sT x4 x5 ∉ [F8.sT x3 x4, F8.sT x2 x3, F8.sT x1 x2, F8.sT 0 x1] := (List.nodup_cons.1 (nodup_rev_take hS 5)).1
  have s6 : F8.sT x5 x6 ∉ [F8.sT x4 x5, F8.sT x3 x4, F8.sT x2 x3, F8.sT x1 x2, F8.sT 0 x1] := (List.nodup_cons.1 (nodup_rev_take hS 6)).1
  have s7 : F8.sT x6 x7 ∉ [F8.sT x5 x6, F8.sT x4 x5, F8.sT x3 x4, F8.sT x2 x3, F8.sT x1 x2, F8.sT 0 x1] := (List.nodup_cons.1 (nodup_rev_take hS 7)).1
  have sc : F8.sT x7 0 ∉ [F8.sT x6 x7, F8.sT x5 x6, F8.sT x4 x5, F8.sT x3 x4, F8.sT x2 x3, F8.sT x1 x2, F8.sT 0 x1] := (List.nodup_cons.1 (nodup_rev_take hS 8)).1
  exact F8.no_order_zero x1 (F8.mem_LL x1) g1 x2 (F8.mem_LL x2) g2 s2 x3 (F8.mem_LL x3) g3 s3 x4 (F8.mem_LL x4) g4 s4 x5 (F8.mem_LL x5) g5 s5 x6 (F8.mem_LL x6) g6 s6 x7 (F8.mem_LL x7) g7 s7 sc

/-- F₉, `P(x) = x^3` (Frobenius, not of the form `c - x`): no circular order works. -/
theorem not_conj_A229827_i_F9 : ¬ Conj_A229827_i := by
  intro h
  have hcard : Fintype.card F9 = 9 := F9.card_eq
  have hdeg : (X ^ 3 : F9[X]).degree < ((Fintype.card F9 - 1 : ℕ) : WithBot ℕ) := by
    rw [degree_X_pow, hcard]; decide
  have key : ∀ c : F9, ∃ x : F9, x ^ 3 ≠ c - x := by decide +kernel
  have hnot : ∀ c : F9, (X ^ 3 : F9[X]) ≠ C c - X := by
    intro c hc
    obtain ⟨x, hx⟩ := key c
    apply hx
    simpa using congrArg (eval x) hc
  obtain ⟨a, ha, hall, hcov⟩ := h F9 (X ^ 3) hdeg hnot
  simp only [eval_X_pow] at hcov
  have huniv : a.toFinset = Finset.univ :=
    Finset.eq_univ_iff_forall.2 (fun x => List.mem_toFinset.2 (hall x))
  have hlen : a.length = 9 := by
    rw [← List.toFinset_card_of_nodup ha, huniv, Finset.card_univ, hcard]
  obtain ⟨k, hk⟩ := exists_rotate_head (hall 0)
  have hb : (a.rotate k).Nodup := List.nodup_rotate.2 ha
  have hlenb : (a.rotate k).length = 9 := by rw [List.length_rotate, hlen]
  have hcovb : ∀ y : F9, y ∈ List.zipWith (fun x z => x ^ 3 + z) (a.rotate k) ((a.rotate k).rotate 1) := by
    intro y; rw [sums_rotate]; exact List.mem_rotate.2 (hcov y)
  have hS := nodup_of_cover _ (by simp only [List.length_zipWith, List.length_rotate, hlen, min_self, hcard]) hcovb
  have hfun : (fun x z : F9 => x ^ 3 + z) = F9.sT := by
    funext x z; exact (F9.sT_eq x z).symm
  rw [hfun] at hS
  obtain ⟨x0, x1, x2, x3, x4, x5, x6, x7, x8, hx⟩ := eq_of_length_nine _ hlenb
  rw [hx] at hb hS hk
  simp only [List.head?_cons, Option.some.injEq] at hk
  subst hk
  have g1 : x1 ∉ [0] := (List.nodup_cons.1 (nodup_rev_take hb 2)).1
  have g2 : x2 ∉ [x1, 0] := (List.nodup_cons.1 (nodup_rev_take hb 3)).1
  have g3 : x3 ∉ [x2, x1, 0] := (List.nodup_cons.1 (nodup_rev_take hb 4)).1
  have g4 : x4 ∉ [x3, x2, x1, 0] := (List.nodup_cons.1 (nodup_rev_take hb 5)).1
  have g5 : x5 ∉ [x4, x3, x2, x1, 0] := (List.nodup_cons.1 (nodup_rev_take hb 6)).1
  have g6 : x6 ∉ [x5, x4, x3, x2, x1, 0] := (List.nodup_cons.1 (nodup_rev_take hb 7)).1
  have g7 : x7 ∉ [x6, x5, x4, x3, x2, x1, 0] := (List.nodup_cons.1 (nodup_rev_take hb 8)).1
  have g8 : x8 ∉ [x7, x6, x5, x4, x3, x2, x1, 0] := (List.nodup_cons.1 (nodup_rev_take hb 9)).1
  have s2 : F9.sT x1 x2 ∉ [F9.sT 0 x1] := (List.nodup_cons.1 (nodup_rev_take hS 2)).1
  have s3 : F9.sT x2 x3 ∉ [F9.sT x1 x2, F9.sT 0 x1] := (List.nodup_cons.1 (nodup_rev_take hS 3)).1
  have s4 : F9.sT x3 x4 ∉ [F9.sT x2 x3, F9.sT x1 x2, F9.sT 0 x1] := (List.nodup_cons.1 (nodup_rev_take hS 4)).1
  have s5 : F9.sT x4 x5 ∉ [F9.sT x3 x4, F9.sT x2 x3, F9.sT x1 x2, F9.sT 0 x1] := (List.nodup_cons.1 (nodup_rev_take hS 5)).1
  have s6 : F9.sT x5 x6 ∉ [F9.sT x4 x5, F9.sT x3 x4, F9.sT x2 x3, F9.sT x1 x2, F9.sT 0 x1] := (List.nodup_cons.1 (nodup_rev_take hS 6)).1
  have s7 : F9.sT x6 x7 ∉ [F9.sT x5 x6, F9.sT x4 x5, F9.sT x3 x4, F9.sT x2 x3, F9.sT x1 x2, F9.sT 0 x1] := (List.nodup_cons.1 (nodup_rev_take hS 7)).1
  have s8 : F9.sT x7 x8 ∉ [F9.sT x6 x7, F9.sT x5 x6, F9.sT x4 x5, F9.sT x3 x4, F9.sT x2 x3, F9.sT x1 x2, F9.sT 0 x1] := (List.nodup_cons.1 (nodup_rev_take hS 8)).1
  have sc : F9.sT x8 0 ∉ [F9.sT x7 x8, F9.sT x6 x7, F9.sT x5 x6, F9.sT x4 x5, F9.sT x3 x4, F9.sT x2 x3, F9.sT x1 x2, F9.sT 0 x1] := (List.nodup_cons.1 (nodup_rev_take hS 9)).1
  exact F9.no_order_zero x1 (F9.mem_LL x1) g1 x2 (F9.mem_LL x2) g2 s2 x3 (F9.mem_LL x3) g3 s3 x4 (F9.mem_LL x4) g4 s4 x5 (F9.mem_LL x5) g5 s5 x6 (F9.mem_LL x6) g6 s6 x7 (F9.mem_LL x7) g7 s7 x8 (F9.mem_LL x8) g8 s8 sc

theorem not_conj_A229827_i : ¬ Conj_A229827_i := not_conj_A229827_i_F8

end SunLean

#print axioms SunLean.not_conj_A229827_i_F8
#print axioms SunLean.not_conj_A229827_i_F9
