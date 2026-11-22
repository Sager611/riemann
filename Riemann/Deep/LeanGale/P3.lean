import Riemann.Deep.LeanGale.Pts
import Riemann.Deep.LeanGale.Data.P3

/-!
# LeanGale: the 170-point set of Riemann/Deep/Gale/points.txt in ℝ^3

Kernel checks (`decide +kernel`): the point data (`ok`), the spread masks of G(≥ dc²) (`sp_ok`) and the
refutation certificate (`run_ok`). Hence: 170 points, every partition into 4 parts has a part with two points
at squared distance ≥ dc² (`lower`), some partition into 4 parts has all same-part squared distances ≤ dc²
(`upper`); D² = 10000000211897483/10¹⁶, dc² = 7932058815501876/10¹⁶. (Our R³ record set, discoveries/gale-r3.md.)
-/

namespace LeanGale.P3

open LeanGale Cert Gale7

/-- Riemann/Deep/Gale/points.txt -/
def S : PtSet := ⟨170, 3, pts, D2, dc2, iD, jD, col⟩

theorem ok : S.ok 4 = true := by decide +kernel

theorem sp_ok : SpreadOK sp 170 S.E = true := by decide +kernel

theorem run_ok : run sp 170 4 toks [initS 4 170] = true := by decide +kernel

/-- the 170-point set -/
noncomputable def X : Finset (Fin 3 → ℝ) := S.X

theorem card_X : X.card = 170 := PtSet.card_X ok

theorem diam_X : diam X = Real.sqrt ((10000000211897483 : ℝ) / 10 ^ 16) := by
  have := PtSet.diam_X (S := S) ok
  norm_num [S, D2] at this ⊢
  exact this

theorem sqd_le (x y : Fin 3 → ℝ) (hx : x ∈ X) (hy : y ∈ X) :
    sqd x y ≤ (10000000211897483 : ℝ) / 10 ^ 16 := by
  have := PtSet.sqd_le ok hx hy
  norm_num [S, D2] at this ⊢
  exact this

theorem lower (c : (Fin 3 → ℝ) → Fin 4) :
    ∃ x ∈ X, ∃ y ∈ X, c x = c y ∧ (7932058815501876 : ℝ) / 10 ^ 16 ≤ sqd x y := by
  have := PtSet.lower_of_cert (S := S) ok (by norm_num) (by norm_num) sp_ok run_ok c
  norm_num [S, dc2] at this ⊢
  exact this

theorem upper : ∃ c : (Fin 3 → ℝ) → Fin 4,
    ∀ x ∈ X, ∀ y ∈ X, c x = c y → sqd x y ≤ (7932058815501876 : ℝ) / 10 ^ 16 := by
  have := PtSet.upper (S := S) ok (by norm_num)
  norm_num [S, dc2] at this ⊢
  exact this

end LeanGale.P3
