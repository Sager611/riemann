import Riemann.Deep.LeanGale.Pts
import Riemann.Deep.LeanGale.Data.P8

/-!
# LeanGale: the 85-point set of points8.txt in ℝ⁸ (claim (1), n = 8)

points8.txt is the family member V₈ on the 10⁻⁸ grid. Kernel checks (`decide +kernel`): the point data
(`ok`); the labels `idx` (P_i, Q_j, A_ij, B_ij, C_ij) send every edge of H(5, 5) to a pair at squared distance
≥ dc² (`H_ok`). Hence, by the colouring lemma (Gale7 `Vx.not_colourable`): every partition into 9 parts has a
part with two points at squared distance ≥ dc² (`lower`); some partition into 9 parts has all same-part
squared distances ≤ dc² (`upper`); D² = 10000000157507010/10¹⁶, dc² = 9023994457279939/10¹⁶.
-/

namespace LeanGale.P8

open LeanGale Cert Gale7 Vx

/-- points8.txt -/
def S : PtSet := ⟨85, 8, pts, D2, dc2, iD, jD, col⟩

theorem ok : S.ok 9 = true := by decide +kernel

/-- the H(5, 5) labels of the points (index of P_i, Q_j, A_ij, B_ij, C_ij in points8.txt) -/
def idx : Vx 5 5 → ℕ
  | P i => idxP.getD i 0
  | Q j => idxQ.getD j 0
  | A i j => idxA.getD (5 * i + j) 0
  | B i j => idxB.getD (5 * i + j) 0
  | C i j => idxC.getD (5 * i + j) 0

theorem idx_ok : (vxList 5 5).all (fun v => decide (idx v < 85)) = true := by decide +kernel

theorem H_ok : (vxList 5 5).all (fun v => (vxList 5 5).all fun w =>
    !adjD (orb v) (orb w) (decide (ix1 v = ix1 w)) (decide (ix2 v = ix2 w)) ||
      decide (dc2 ≤ S.sq (idx v) (idx w))) = true := by decide +kernel

/-- the 85-point set -/
noncomputable def X8 : Finset (Fin 8 → ℝ) := S.X

theorem card_X8 : X8.card = 85 := PtSet.card_X ok

theorem diam_X8 : diam X8 = Real.sqrt ((10000000157507010 : ℝ) / 10 ^ 16) := by
  have := PtSet.diam_X (S := S) ok
  norm_num [S, D2] at this ⊢
  exact this

theorem sqd_le (x y : Fin 8 → ℝ) (hx : x ∈ X8) (hy : y ∈ X8) :
    sqd x y ≤ (10000000157507010 : ℝ) / 10 ^ 16 := by
  have := PtSet.sqd_le ok hx hy
  norm_num [S, D2] at this ⊢
  exact this

theorem lower (c : (Fin 8 → ℝ) → Fin 9) :
    ∃ x ∈ X8, ∃ y ∈ X8, c x = c y ∧ (9023994457279939 : ℝ) / 10 ^ 16 ≤ sqd x y := by
  have := PtSet.lower_of_labels (S := S) (m₁ := 4) (m₂ := 4) ok (by norm_num) idx idx_ok H_ok c
  norm_num [S, dc2] at this ⊢
  exact this

theorem upper : ∃ c : (Fin 8 → ℝ) → Fin 9,
    ∀ x ∈ X8, ∀ y ∈ X8, c x = c y → sqd x y ≤ (9023994457279939 : ℝ) / 10 ^ 16 := by
  have := PtSet.upper (S := S) ok (by norm_num)
  norm_num [S, dc2] at this ⊢
  exact this

end LeanGale.P8
