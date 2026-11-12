import Mathlib.Algebra.Field.Defs
import Mathlib.Data.Fintype.Card
import Mathlib.Data.List.Nodup

/-! The field F₈ = F₂[t]/(t³ + t + 1), element b₀ + b₁ t + b₂ t² ↦ b₀ + 2 b₁ + 4 b₂ as packed tables on `Fin 8` (definitions only). -/

namespace SunLean

/-- F₈ = F₂[t]/(t³ + t + 1), element b₀ + b₁ t + b₂ t² ↦ b₀ + 2 b₁ + 4 b₂. -/
def F8 : Type := Fin 8

namespace F8

instance : DecidableEq F8 := inferInstanceAs (DecidableEq (Fin 8))
instance : Fintype F8 := inferInstanceAs (Fintype (Fin 8))

/-- An element given by its code. -/
def ofFin (i : Fin 8) : F8 := i

/-- Tables packed into one natural number, 4 bits per entry (entry `(a, b)` at bit `4 (8 a + b)`). -/
def ADD : ℕ := 514631494928646303341455936939146235248427799039887685160188751890747765264
def MUL : ℕ := 23691910682714665197238709536041582193718400246807875409196429481426745819136
def NEG : ℕ := 1985229328
def INV : ℕ := 1126655248

def addT (a b : Fin 8) : Fin 8 := ⟨(ADD >>> (4 * (8 * a.val + b.val))) % 16 % 8, Nat.mod_lt _ (by decide)⟩
def mulT (a b : Fin 8) : Fin 8 := ⟨(MUL >>> (4 * (8 * a.val + b.val))) % 16 % 8, Nat.mod_lt _ (by decide)⟩
def negT (a : Fin 8) : Fin 8 := ⟨(NEG >>> (4 * a.val)) % 16 % 8, Nat.mod_lt _ (by decide)⟩
def invT (a : Fin 8) : Fin 8 := ⟨(INV >>> (4 * a.val)) % 16 % 8, Nat.mod_lt _ (by decide)⟩

instance : Zero F8 := ⟨(0 : Fin 8)⟩
instance : One F8 := ⟨(1 : Fin 8)⟩
instance : Add F8 := ⟨addT⟩
instance : Mul F8 := ⟨mulT⟩
instance : Neg F8 := ⟨negT⟩
instance : Inv F8 := ⟨invT⟩

end F8

end SunLean
