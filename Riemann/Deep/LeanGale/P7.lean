import Riemann.Deep.LeanGale.Pts
import Riemann.Deep.LeanGale.Data.P7

/-!
# LeanGale: the 84-point set of points7.txt in ℝ⁷ (claim (1), n = 7)

Kernel checks (`decide +kernel`): the point data (`ok`), the spread masks of G(≥ dc²) (`sp_ok`) and the
refutation certificate (`run_ok`, 18887 tokens). Hence: 84 points, diam = √(D²), every partition into 8
parts has a part with two points at squared distance ≥ dc² (`lower`), and some partition into 8 parts has all
same-part squared distances ≤ dc² (`upper`); D² = 10000000213244525/10¹⁶, dc² = 8900856812318928/10¹⁶.
-/

namespace LeanGale.P7

open LeanGale Cert Gale7

/-- points7.txt -/
def S : PtSet := ⟨84, 7, pts, D2, dc2, iD, jD, col⟩

theorem ok : S.ok 8 = true := by decide +kernel

theorem sp_ok : SpreadOK sp 84 S.E = true := by decide +kernel

theorem run_ok : run sp 84 8 toks [initS 8 84] = true := by decide +kernel

/-- the 84-point set -/
noncomputable def X7 : Finset (Fin 7 → ℝ) := S.X

theorem card_X7 : X7.card = 84 := PtSet.card_X ok

theorem diam_X7 : diam X7 = Real.sqrt ((10000000213244525 : ℝ) / 10 ^ 16) := by
  have := PtSet.diam_X (S := S) ok
  norm_num [S, D2] at this ⊢
  exact this

theorem sqd_le (x y : Fin 7 → ℝ) (hx : x ∈ X7) (hy : y ∈ X7) :
    sqd x y ≤ (10000000213244525 : ℝ) / 10 ^ 16 := by
  have := PtSet.sqd_le ok hx hy
  norm_num [S, D2] at this ⊢
  exact this

theorem lower (c : (Fin 7 → ℝ) → Fin 8) :
    ∃ x ∈ X7, ∃ y ∈ X7, c x = c y ∧ (8900856812318928 : ℝ) / 10 ^ 16 ≤ sqd x y := by
  have := PtSet.lower_of_cert (S := S) ok (by norm_num) (by norm_num) sp_ok run_ok c
  norm_num [S, dc2] at this ⊢
  exact this

theorem upper : ∃ c : (Fin 7 → ℝ) → Fin 8,
    ∀ x ∈ X7, ∀ y ∈ X7, c x = c y → sqd x y ≤ (8900856812318928 : ℝ) / 10 ^ 16 := by
  have := PtSet.upper (S := S) ok (by norm_num)
  norm_num [S, dc2] at this ⊢
  exact this

end LeanGale.P7
