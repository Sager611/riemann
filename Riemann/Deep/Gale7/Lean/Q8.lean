import Riemann.Deep.Gale7.Lean.Desc

/-!
# Gale7: the rational class checks for V8 (kernel evaluation, `decide +kernel`)

Parameters of family_params.txt row n = 8 (decimals as exact rationals), split 4 + 4. Over all 5·5·2·2
descriptors: every class value ≤ D8; every edge class of `H` (`adjD`) ≥ t8; every class of distinct vertices
(`neD`) > 0.
-/

namespace Gale7

def rA8 : ℚ := 632455532034 / 1000000000000
def rB8 : ℚ := 632455532034 / 1000000000000
def a8 : ℚ := 627253113456 / 1000000000000
def b8 : ℚ := 80953886008 / 1000000000000
def c8 : ℚ := 96579073850 / 1000000000000
def d8 : ℚ := 609022174896 / 1000000000000
def e8 : ℚ := 404230564733 / 1000000000000
def f8 : ℚ := 481355888451 / 1000000000000
/-- the smallest edge class (`A_ij C_kl`, i ≠ k, j ≠ l) -/
def t8 : ℚ := 90239946609080831680751 / 100000000000000000000000
/-- the largest class: D² -/
def D8 : ℚ := 100000000000102500044289 / 100000000000000000000000

theorem cls8_le : ∀ o o' : Orb, ∀ e₁ e₂ : Bool,
    clsQ 4 4 rA8 rB8 a8 b8 c8 d8 e8 f8 o o' e₁ e₂ ≤ D8 := by
  decide +kernel

theorem cls8_ge : ∀ o o' : Orb, ∀ e₁ e₂ : Bool, adjD o o' e₁ e₂ = true →
    t8 ≤ clsQ 4 4 rA8 rB8 a8 b8 c8 d8 e8 f8 o o' e₁ e₂ := by
  decide +kernel

theorem cls8_pos : ∀ o o' : Orb, ∀ e₁ e₂ : Bool, neD o o' e₁ e₂ = true →
    0 < clsQ 4 4 rA8 rB8 a8 b8 c8 d8 e8 f8 o o' e₁ e₂ := by
  decide +kernel

end Gale7
