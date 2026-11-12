import Riemann.Deep.SunComplete.Lean.Small

/-!
# sunlean · Cycles: `q = 4` works; certificate cycles for `q = 11` and `q = 13` (`decide`)
-/

namespace SunLean

/-- The field `F_4 = F_2[ω]/(ω² + ω + 1)`: `o = 0`, `i = 1`, `w = ω`, `w2 = ω² = ω + 1`. -/
inductive F4 | o | i | w | w2
  deriving DecidableEq

namespace F4

instance : Fintype F4 where
  elems := {o, i, w, w2}
  complete := fun x => by cases x <;> decide

/-- Addition: `F_2`-coordinates in the basis `(1, ω)`. -/
def add : F4 → F4 → F4
  | o, x => x
  | x, o => x
  | i, i => o
  | i, w => w2
  | i, w2 => w
  | w, i => w2
  | w, w => o
  | w, w2 => i
  | w2, i => w
  | w2, w => i
  | w2, w2 => o

/-- Multiplication: `ω³ = 1`. -/
def mul : F4 → F4 → F4
  | o, _ => o
  | _, o => o
  | i, x => x
  | x, i => x
  | w, w => w2
  | w, w2 => i
  | w2, w => i
  | w2, w2 => w

def inv : F4 → F4
  | o => o
  | i => i
  | w => w2
  | w2 => w

instance : Zero F4 := ⟨o⟩
instance : One F4 := ⟨i⟩
instance : Add F4 := ⟨add⟩
instance : Mul F4 := ⟨mul⟩
instance : Neg F4 := ⟨fun x => x⟩
instance : Inv F4 := ⟨inv⟩

instance : Field F4 where
  add_assoc := by decide
  zero_add := by decide
  add_zero := by decide
  add_comm := by decide
  neg_add_cancel := by decide
  mul_assoc := by decide
  one_mul := by decide
  mul_one := by decide
  left_distrib := by decide
  right_distrib := by decide
  zero_mul := by decide
  mul_zero := by decide
  mul_comm := by decide
  exists_pair_ne := ⟨o, i, by decide⟩
  mul_inv_cancel := by decide
  inv_zero := by decide
  nsmul := nsmulRec
  zsmul := zsmulRec
  nnqsmul := _
  nnqsmul_def := fun _ _ => rfl
  qsmul := _
  qsmul_def := fun _ _ => rfl

theorem card : Fintype.card F4 = 4 := rfl

theorem isPrim_iff (x : F4) : IsPrim x ↔ x ^ 3 = 1 ∧ ∀ m, m < 3 → 0 < m → x ^ m ≠ 1 := by
  unfold IsPrim
  rw [card]
  exact orderOf_eq_iff (by decide)

/-- The order `(0, ω, 1, ω²)` (NOTE §0). -/
def cyc : Fin 4 → F4 := ![o, w, i, w2]
def pos : F4 → Fin 4
  | o => 0
  | w => 1
  | i => 2
  | w2 => 3

theorem pos_cyc : ∀ x : Fin 4, pos (cyc x) = x := by decide
theorem cyc_pos : ∀ y : F4, cyc (pos y) = y := by decide

def e : Fin 4 ≃ F4 := ⟨cyc, pos, pos_cyc, cyc_pos⟩

end F4

/-- **`q = 4` works**: `(0, ω, 1, ω²)` has primitive cyclic adjacent sums. -/
theorem sunHou_F4 : SunHou F4 := by
  refine sunHou_of_cycle F4.card F4.e ?_
  simp only [finRotate_apply, F4.isPrim_iff]
  decide

/-- `q = 11`: She's construction (Lemma 8) with `x = 6` (certA.json `"11": she, x = 6`),
`(a_1, …, a_11) = (4, 2, 6, 0, 8, 9, 10, 7, 1, 5, 3)`. -/
def cyc11 : Fin 11 → ZMod 11 := ![4, 2, 6, 0, 8, 9, 10, 7, 1, 5, 3]
def pos11 : ZMod 11 → Fin 11 := ![3, 8, 1, 10, 0, 9, 2, 7, 4, 5, 6]
theorem pos_cyc11 : ∀ x : Fin 11, pos11 (cyc11 x) = x := by decide
theorem cyc_pos11 : ∀ y : ZMod 11, cyc11 (pos11 y) = y := by decide
def e11 : Fin 11 ≃ ZMod 11 := ⟨cyc11, pos11, pos_cyc11, cyc_pos11⟩

theorem sunHou_zmod11 : SunHou (ZMod 11) := by
  refine sunHou_of_cycle (ZMod.card 11) e11 ?_
  simp only [finRotate_apply, isPrim_zmod_iff]
  decide

/-- `q = 13`: the explicit cycle of certA.json. -/
def cyc13 : Fin 13 → ZMod 13 := ![1, 10, 5, 2, 9, 11, 4, 3, 8, 12, 7, 0, 6]
def pos13 : ZMod 13 → Fin 13 := ![11, 0, 3, 7, 6, 2, 12, 10, 8, 4, 1, 5, 9]
theorem pos_cyc13 : ∀ x : Fin 13, pos13 (cyc13 x) = x := by decide
theorem cyc_pos13 : ∀ y : ZMod 13, cyc13 (pos13 y) = y := by decide
def e13 : Fin 13 ≃ ZMod 13 := ⟨cyc13, pos13, pos_cyc13, cyc_pos13⟩

theorem sunHou_zmod13 : SunHou (ZMod 13) := by
  refine sunHou_of_cycle (ZMod.card 13) e13 ?_
  simp only [finRotate_apply, isPrim_zmod_iff]
  decide

/-- Every field with `11` or `13` elements has the Sun–Hou property. -/
theorem sunHou_11_13 (F : Type) [Field F] [Fintype F]
    (h : Fintype.card F = 11 ∨ Fintype.card F = 13) : SunHou F := by
  rcases h with h | h
  · exact sunHou_of_ringEquiv (ZMod.ringEquivOfPrime F (by decide) h) sunHou_zmod11
  · exact sunHou_of_ringEquiv (ZMod.ringEquivOfPrime F (by decide) h) sunHou_zmod13

end SunLean
