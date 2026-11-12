import Riemann.Deep.SunComplete.Lean.Small
import Mathlib.Tactic.Ring

/-!
# sunlean · Cycles9: the certificate cycle for `q = 9` (certA.json), kernel-checked by `decide`

`F_9 = F_3[x]/(x² + x + 2)` (certA.json `"f": [2, 1, 1]`); `a + b·x` is coded `a + 3b` there.
-/

namespace SunLean

/-- `a + b·x` in `F_3[x]/(x² + x + 2)`, so `x² = 1 + 2x`. -/
structure F9 where
  a : ZMod 3
  b : ZMod 3
  deriving DecidableEq

namespace F9

instance : Fintype F9 :=
  Fintype.ofEquiv (ZMod 3 × ZMod 3) ⟨fun p => ⟨p.1, p.2⟩, fun x => (x.a, x.b), fun _ => rfl, fun _ => rfl⟩

def add (x y : F9) : F9 := ⟨x.a + y.a, x.b + y.b⟩
/-- `(a + bx)(c + dx) = (ac + bd) + (ad + bc + 2bd)x`. -/
def mul (x y : F9) : F9 := ⟨x.a * y.a + x.b * y.b, x.a * y.b + x.b * y.a + 2 * x.b * y.b⟩
/-- `x⁻¹ = x⁷` (`x⁸ = 1` for `x ≠ 0`; `0⁷ = 0`). -/
def inv (x : F9) : F9 := mul x (mul x (mul x (mul x (mul x (mul x x)))))

instance : Zero F9 := ⟨⟨0, 0⟩⟩
instance : One F9 := ⟨⟨1, 0⟩⟩
instance : Add F9 := ⟨add⟩
instance : Mul F9 := ⟨mul⟩
instance : Neg F9 := ⟨fun x => ⟨-x.a, -x.b⟩⟩
instance : Inv F9 := ⟨inv⟩

theorem ext' {x y : F9} (h1 : x.a = y.a) (h2 : x.b = y.b) : x = y := by
  cases x; cases y; simp_all

@[simp] theorem add_def (x y : F9) : x + y = add x y := rfl
@[simp] theorem mul_def (x y : F9) : x * y = mul x y := rfl
@[simp] theorem add_a (x y : F9) : (add x y).a = x.a + y.a := rfl
@[simp] theorem add_b (x y : F9) : (add x y).b = x.b + y.b := rfl
@[simp] theorem mul_a (x y : F9) : (mul x y).a = x.a * y.a + x.b * y.b := rfl
@[simp] theorem mul_b (x y : F9) : (mul x y).b = x.a * y.b + x.b * y.a + 2 * x.b * y.b := rfl

instance : Field F9 where
  add_assoc := fun x y z => by
    apply ext' <;> simp only [add_def, add_a, add_b] <;> ring
  zero_add := by decide
  add_zero := by decide
  add_comm := by decide
  neg_add_cancel := by decide
  mul_assoc := fun x y z => by
    apply ext' <;> simp only [mul_def, mul_a, mul_b] <;> ring
  one_mul := by decide
  mul_one := by decide
  left_distrib := fun x y z => by
    apply ext' <;> simp only [add_def, mul_def, add_a, add_b, mul_a, mul_b] <;> ring
  right_distrib := fun x y z => by
    apply ext' <;> simp only [add_def, mul_def, add_a, add_b, mul_a, mul_b] <;> ring
  zero_mul := by decide
  mul_zero := by decide
  mul_comm := by decide
  exists_pair_ne := ⟨0, 1, by decide⟩
  mul_inv_cancel := by decide
  inv_zero := by decide
  nsmul := nsmulRec
  zsmul := zsmulRec
  nnqsmul := _
  nnqsmul_def := fun _ _ => rfl
  qsmul := _
  qsmul_def := fun _ _ => rfl

theorem card : Fintype.card F9 = 9 := rfl

theorem isPrim_iff (x : F9) : IsPrim x ↔ x ^ 8 = 1 ∧ ∀ m, m < 8 → 0 < m → x ^ m ≠ 1 := by
  unfold IsPrim
  rw [card]
  exact orderOf_eq_iff (by decide)

/-- certA.json `"9"`: `cycle = [3, 0, 4, 5, 2, 6, 7, 8, 1]` in the code `a + 3b`. -/
def cyc : Fin 9 → F9 :=
  ![⟨0, 1⟩, ⟨0, 0⟩, ⟨1, 1⟩, ⟨2, 1⟩, ⟨2, 0⟩, ⟨0, 2⟩, ⟨1, 2⟩, ⟨2, 2⟩, ⟨1, 0⟩]
def pos (x : F9) : Fin 9 := (![![1, 0, 5], ![8, 2, 6], ![4, 3, 7]] : Fin 3 → Fin 3 → Fin 9) x.a x.b

theorem pos_cyc : ∀ i : Fin 9, pos (cyc i) = i := by decide
theorem cyc_pos : ∀ x : F9, cyc (pos x) = x := by decide

def e : Fin 9 ≃ F9 := ⟨cyc, pos, pos_cyc, cyc_pos⟩

end F9

/-- **`q = 9`**: the certA.json cycle has primitive cyclic adjacent sums in `F_9`. -/
theorem sunHou_F9 : SunHou F9 := by
  refine sunHou_of_cycle F9.card F9.e ?_
  simp only [finRotate_apply, F9.isPrim_iff]
  decide

end SunLean
