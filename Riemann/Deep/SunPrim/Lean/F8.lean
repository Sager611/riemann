import Riemann.Deep.SunPrim.Lean.F8Add
import Riemann.Deep.SunPrim.Lean.F8Mul
import Riemann.Deep.SunPrim.Lean.F8DistL
import Riemann.Deep.SunPrim.Lean.F8DistR

/-! F₈ = F₂[t]/(t³ + t + 1), element b₀ + b₁ t + b₂ t² ↦ b₀ + 2 b₁ + 4 b₂ is a field with 8 elements. -/

namespace SunLean

namespace F8

instance instField : Field F8 where
  add := (· + ·)
  add_assoc := add_assoc'
  zero := 0
  zero_add := zero_add'
  add_zero := add_zero'
  nsmul := nsmulRec
  neg := Neg.neg
  zsmul := zsmulRec
  neg_add_cancel := neg_add_cancel'
  add_comm := add_comm'
  mul := (· * ·)
  left_distrib := left_distrib'
  right_distrib := right_distrib'
  zero_mul := zero_mul'
  mul_zero := mul_zero'
  mul_assoc := mul_assoc'
  one := 1
  one_mul := one_mul'
  mul_one := mul_one'
  mul_comm := mul_comm'
  inv := Inv.inv
  exists_pair_ne := ⟨0, 1, zero_ne_one'⟩
  mul_inv_cancel := mul_inv_cancel'
  inv_zero := inv_zero'
  nnqsmul := _
  nnqsmul_def := fun _ _ => rfl
  qsmul := _
  qsmul_def := fun _ _ => rfl

theorem card_eq : Fintype.card F8 = 8 := by decide

end F8

end SunLean

#print axioms SunLean.F8.instField
