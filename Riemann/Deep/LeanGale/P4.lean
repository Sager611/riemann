import Riemann.Deep.LeanGale.Pts
import Riemann.Deep.LeanGale.Data.P4

/-!
# LeanGale: the 126-point set of Riemann/Deep/Gale45/points4.txt in ℝ^4

Kernel checks (`decide +kernel`): the point data (`ok`), the spread masks of G(≥ dc²) (`sp_ok`) and the
refutation certificate (`run_ok`). Hence: 126 points, every partition into 5 parts has a part with two points
at squared distance ≥ dc² (`lower`), some partition into 5 parts has all same-part squared distances ≤ dc²
(`upper`); D² = 10000000179265409/10¹⁶, dc² = 8399124363373578/10¹⁶. (Our R⁴ record set, discoveries/gale-r4-r5.md.)
-/

namespace LeanGale.P4

open LeanGale Cert Gale7

/-- Riemann/Deep/Gale45/points4.txt -/
def S : PtSet := ⟨126, 4, pts, D2, dc2, iD, jD, col⟩

theorem ok : S.ok 5 = true := by decide +kernel

theorem sp_ok : SpreadOK sp 126 S.E = true := by decide +kernel

theorem run_ok : run sp 126 5 toks [initS 5 126] = true := by decide +kernel

/-- the 126-point set -/
noncomputable def X : Finset (Fin 4 → ℝ) := S.X

theorem card_X : X.card = 126 := PtSet.card_X ok

theorem diam_X : diam X = Real.sqrt ((10000000179265409 : ℝ) / 10 ^ 16) := by
  have := PtSet.diam_X (S := S) ok
  norm_num [S, D2] at this ⊢
  exact this

theorem sqd_le (x y : Fin 4 → ℝ) (hx : x ∈ X) (hy : y ∈ X) :
    sqd x y ≤ (10000000179265409 : ℝ) / 10 ^ 16 := by
  have := PtSet.sqd_le ok hx hy
  norm_num [S, D2] at this ⊢
  exact this

theorem lower (c : (Fin 4 → ℝ) → Fin 5) :
    ∃ x ∈ X, ∃ y ∈ X, c x = c y ∧ (8399124363373578 : ℝ) / 10 ^ 16 ≤ sqd x y := by
  have := PtSet.lower_of_cert (S := S) ok (by norm_num) (by norm_num) sp_ok run_ok c
  norm_num [S, dc2] at this ⊢
  exact this

theorem upper : ∃ c : (Fin 4 → ℝ) → Fin 5,
    ∀ x ∈ X, ∀ y ∈ X, c x = c y → sqd x y ≤ (8399124363373578 : ℝ) / 10 ^ 16 := by
  have := PtSet.upper (S := S) ok (by norm_num)
  norm_num [S, dc2] at this ⊢
  exact this

end LeanGale.P4
