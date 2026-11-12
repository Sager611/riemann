import Mathlib.Algebra.Field.Defs
import Mathlib.Data.Fintype.Card
import Mathlib.Data.List.Nodup

/-! The field F₉ = F₃[i]/(i² + 1), element a + b i ↦ a + 3 b as packed tables on `Fin 9` (definitions only). -/

namespace SunLean

/-- F₉ = F₃[i]/(i² + 1), element a + b i ↦ a + 3 b. -/
def F9 : Type := Fin 9

namespace F9

instance : DecidableEq F9 := inferInstanceAs (DecidableEq (Fin 9))
instance : Fintype F9 := inferInstanceAs (Fintype (Fin 9))

/-- An element given by its code. -/
def ofFin (i : Fin 9) : F9 := i

/-- Tables packed into one natural number, 4 bits per entry (entry `(a, b)` at bit `4 (9 a + b)`). -/
def ADD : ℕ := 8986690708024339884259668388453865091868182506470258411551265922896781674643265104972283578167824
def MUL : ℕ := 12992295474418726299168584010125499800446975643965808714274855111902658395364003326119775607717888
def NEG : ℕ := 18580267296
def INV : ℕ := 32267133456

def addT (a b : Fin 9) : Fin 9 := ⟨(ADD >>> (4 * (9 * a.val + b.val))) % 16 % 9, Nat.mod_lt _ (by decide)⟩
def mulT (a b : Fin 9) : Fin 9 := ⟨(MUL >>> (4 * (9 * a.val + b.val))) % 16 % 9, Nat.mod_lt _ (by decide)⟩
def negT (a : Fin 9) : Fin 9 := ⟨(NEG >>> (4 * a.val)) % 16 % 9, Nat.mod_lt _ (by decide)⟩
def invT (a : Fin 9) : Fin 9 := ⟨(INV >>> (4 * a.val)) % 16 % 9, Nat.mod_lt _ (by decide)⟩

instance : Zero F9 := ⟨(0 : Fin 9)⟩
instance : One F9 := ⟨(1 : Fin 9)⟩
instance : Add F9 := ⟨addT⟩
instance : Mul F9 := ⟨mulT⟩
instance : Neg F9 := ⟨negT⟩
instance : Inv F9 := ⟨invT⟩

end F9

end SunLean
