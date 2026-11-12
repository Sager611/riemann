import Riemann.Deep.SunComplete.Lean.Graph
import Mathlib.GroupTheory.OrderOfElement
import Mathlib.Algebra.Field.Defs

/-!
# sunlean · Defs: primitive elements, the Sun–Hou property, the graph `H_q`

`F` is a finite field (`Field F`, `Fintype F`), `q = Fintype.card F`.
-/

namespace SunLean

section
variable {F : Type} [Field F] [Fintype F]

/-- `x` is a primitive element of the finite field `F` (`q = card F`): a generator of the cyclic
group `Fˣ` of order `q − 1`, i.e. an element of multiplicative order `q − 1`. -/
def IsPrim (x : F) : Prop := orderOf x = Fintype.card F - 1

variable (F) in
/-- **The Sun–Hou property** (Sun–Hou, Conjecture 3.5(i)): the elements of `F` can be numbered
`a_1, …, a_q` (here `a 0, …, a (q-1)`, a bijection `Fin q ≃ F`) so that all `q` cyclic adjacent
sums `a_1 + a_2, …, a_{q-1} + a_q, a_q + a_1` are primitive. -/
def SunHou : Prop := ∃ a : Fin (Fintype.card F) ≃ F, ∀ i, IsPrim (a i + a (finRotate _ i))

variable (F) in
/-- The graph `H_q` (NOTE §1): vertex set `F`, `x ~ y` iff `x ≠ y` and `x + y` is primitive. -/
def sumGraph : SimpleGraph F where
  Adj x y := x ≠ y ∧ IsPrim (x + y)
  symm := ⟨fun x y h => ⟨h.1.symm, by rw [add_comm]; exact h.2⟩⟩
  loopless := ⟨fun x h => h.1 rfl⟩

/-- For `q ≥ 2` the Sun–Hou property is exactly a Hamiltonian cycle of `H_q`. -/
theorem sunHou_iff (h2 : 2 ≤ Fintype.card F) : SunHou F ↔ HasHamCycle (sumGraph F) := by
  constructor
  · rintro ⟨a, ha⟩
    exact ⟨a, fun i => ⟨fun h => finRotate_ne h2 i (a.injective h).symm, ha i⟩⟩
  · rintro ⟨a, ha⟩
    exact ⟨a, fun i => (ha i).2⟩

/-- `0` is not primitive (`q ≥ 2`). -/
theorem not_isPrim_zero (h : 2 ≤ Fintype.card F) : ¬ IsPrim (0 : F) := by
  intro h0
  have h1 := pow_orderOf_eq_one (0 : F)
  unfold IsPrim at h0
  rw [h0, zero_pow (by omega)] at h1
  exact zero_ne_one h1

/-- In a Sun–Hou numbering (`q ≥ 3`) every element `v` has two distinct neighbours `y ≠ z`
(`≠ v`) with `v + y`, `v + z` primitive. -/
theorem two_nbrs (h3 : 3 ≤ Fintype.card F) (hS : SunHou F) (v : F) :
    ∃ y z : F, y ≠ z ∧ y ≠ v ∧ z ≠ v ∧ IsPrim (v + y) ∧ IsPrim (v + z) := by
  obtain ⟨a, ha⟩ := hS
  refine ⟨a (finRotate _ (a.symm v)), a ((finRotate _).symm (a.symm v)), ?_, ?_, ?_, ?_, ?_⟩
  · exact fun h => finRotate_ne_symm h3 _ (a.injective h)
  · intro h
    exact finRotate_ne (by omega) _ (a.injective (h.trans (a.apply_symm_apply v).symm))
  · intro h
    have h' : (finRotate _).symm (a.symm v) = a.symm v :=
      a.injective (h.trans (a.apply_symm_apply v).symm)
    have h'' := congrArg (finRotate _) h'
    rw [Equiv.apply_symm_apply] at h''
    exact finRotate_ne (by omega) _ h''.symm
  · have := ha (a.symm v)
    rwa [a.apply_symm_apply] at this
  · have := ha ((finRotate _).symm (a.symm v))
    rw [Equiv.apply_symm_apply, a.apply_symm_apply, add_comm] at this
    exact this

variable (F) in
/-- The set `P` of primitive elements. -/
noncomputable def primSet : Finset F := by
  classical exact Finset.univ.filter (fun x : F => IsPrim x)

theorem mem_primSet {x : F} : x ∈ primSet F ↔ IsPrim x := by
  unfold primSet
  simp

noncomputable instance [DecidableEq F] : DecidableRel (sumGraph F).Adj := fun x y =>
  decidable_of_iff (x ≠ y ∧ x + y ∈ primSet F) (by rw [mem_primSet]; rfl)

end

end SunLean
