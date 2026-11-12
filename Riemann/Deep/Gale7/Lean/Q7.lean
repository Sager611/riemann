import Riemann.Deep.Gale7.Lean.Desc

/-!
# Gale7: the rational class checks for V7 (kernel evaluation, `decide +kernel`)

Parameters of family_params.txt row n = 7 (decimals as exact rationals), split 3 + 4. Over all 5·5·2·2
descriptors: every class value ≤ D7; every edge class of `H` (`adjD`) ≥ t7; every class of distinct vertices
(`neD`) > 0.
-/

namespace Gale7

def rA7 : ℚ := 612372435696 / 1000000000000
def rB7 : ℚ := 632455532034 / 1000000000000
def a7 : ℚ := 606069108095 / 1000000000000
def b7 : ℚ := 90511059877 / 1000000000000
def c7 : ℚ := 105433773038 / 1000000000000
def d7 : ℚ := 612172119187 / 1000000000000
def e7 : ℚ := 381914257728 / 1000000000000
def f7 : ℚ := 486126417449 / 1000000000000
/-- the smallest edge class (`A_ij C_kl`, i ≠ k, j ≠ l) -/
def t7 : ℚ := 356000182722431308912917 / 400000000000000000000000
/-- the largest class: D² -/
def D7 : ℚ := 100000000000102500044289 / 100000000000000000000000

theorem cls7_le : ∀ o o' : Orb, ∀ e₁ e₂ : Bool,
    clsQ 3 4 rA7 rB7 a7 b7 c7 d7 e7 f7 o o' e₁ e₂ ≤ D7 := by
  decide +kernel

theorem cls7_ge : ∀ o o' : Orb, ∀ e₁ e₂ : Bool, adjD o o' e₁ e₂ = true →
    t7 ≤ clsQ 3 4 rA7 rB7 a7 b7 c7 d7 e7 f7 o o' e₁ e₂ := by
  decide +kernel

theorem cls7_pos : ∀ o o' : Orb, ∀ e₁ e₂ : Bool, neD o o' e₁ e₂ = true →
    0 < clsQ 3 4 rA7 rB7 a7 b7 c7 d7 e7 f7 o o' e₁ e₂ := by
  decide +kernel

end Gale7
