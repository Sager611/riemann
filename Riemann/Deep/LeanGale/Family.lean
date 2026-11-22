import Riemann.Deep.LeanGale.Data.Rows00
import Riemann.Deep.LeanGale.Data.Rows01
import Riemann.Deep.LeanGale.Data.Rows02
import Riemann.Deep.LeanGale.Data.Rows03
import Riemann.Deep.LeanGale.Data.Rows04
import Riemann.Deep.LeanGale.Data.Rows05
import Riemann.Deep.LeanGale.Data.Rows06
import Riemann.Deep.LeanGale.Data.Rows07
import Riemann.Deep.LeanGale.Data.Rows08
import Riemann.Deep.LeanGale.Data.Rows09
import Riemann.Deep.LeanGale.Data.Rows10
import Riemann.Deep.LeanGale.Data.Rows11
import Riemann.Deep.LeanGale.Data.Rows12
import Riemann.Deep.LeanGale.Data.Rows13
import Riemann.Deep.LeanGale.Data.Rows14
import Riemann.Deep.LeanGale.Data.Rows15
import Riemann.Deep.LeanGale.Data.Rows16
import Riemann.Deep.LeanGale.Data.Rows17
import Riemann.Deep.LeanGale.Data.Rows18
import Riemann.Deep.LeanGale.Data.Rows19

/-!
# LeanGale: the family for every n = 5…1000 (claim (2))

`table`: the 996 rows of family_params.txt (n = 5…1000), each kernel-checked by `rowOK` in Data/RowsNN.lean.
`params n`: the row of n. `famSet n U W`: the family member Vₙ ⊂ ℝ^⌊n/2⌋ ⊕ ℝ^⌈n/2⌉ (NOTE §2) for unit
simplices U, W. `famSet_not_galeSplit`: for 5 ≤ n ≤ 1000, Vₙ has no partition into n + 1 parts of diameter
≤ σₙ · diam; `card_famSet`: it has (⌊n/2⌋+1) + (⌈n/2⌉+1) + 3(⌊n/2⌋+1)(⌈n/2⌉+1) points.
-/

namespace LeanGale

open Gale7 Vx

/-- H(p, q) has p + q + 3pq vertices (as in Gale7.Main) -/
theorem card_Vx' (p q : ℕ) : Fintype.card (Vx p q) = p + q + 3 * (p * q) := by
  have e : Vx p q ≃ Fin p ⊕ Fin q ⊕ (Fin p × Fin q) ⊕ (Fin p × Fin q) ⊕ (Fin p × Fin q) :=
    { toFun := fun v => match v with
        | P i => .inl i
        | Q j => .inr (.inl j)
        | A i j => .inr (.inr (.inl (i, j)))
        | B i j => .inr (.inr (.inr (.inl (i, j))))
        | C i j => .inr (.inr (.inr (.inr (i, j))))
      invFun := fun s => match s with
        | .inl i => P i
        | .inr (.inl j) => Q j
        | .inr (.inr (.inl (i, j))) => A i j
        | .inr (.inr (.inr (.inl (i, j)))) => B i j
        | .inr (.inr (.inr (.inr (i, j)))) => C i j
      left_inv := fun v => by cases v <;> rfl
      right_inv := fun s => by
        rcases s with i | j | ⟨i, j⟩ | ⟨i, j⟩ | ⟨i, j⟩ <;> rfl }
  rw [Fintype.card_congr e]
  simp only [Fintype.card_sum, Fintype.card_prod, Fintype.card_fin]
  ring

/-- the three metric facts behind `rowOK_sound`, for the family member of an accepted row -/
theorem rowOK_facts {m₁ m₂ : ℕ} {p : Params} (h : rowOK m₁ m₂ p = true)
    {U : Fin (m₁ + 1) → Fin m₁ → ℝ} {W : Fin (m₂ + 1) → Fin m₂ → ℝ}
    (hU : IsUnitSimplex U) (hW : IsUnitSimplex W) :
    (∀ x ∈ Finset.univ.image (famRow U W p), ∀ y ∈ Finset.univ.image (famRow U W p),
      sqd x y ≤ (p.DZ : ℝ) / ((m₁ : ℝ) * m₂ * 10 ^ 24)) ∧
    (∀ c : (Fin (m₁ + m₂) → ℝ) → Fin (m₁ + m₂ + 1), ∃ x ∈ Finset.univ.image (famRow U W p),
      ∃ y ∈ Finset.univ.image (famRow U W p), c x = c y ∧
        (p.tZ : ℝ) / ((m₁ : ℝ) * m₂ * 10 ^ 24) ≤ sqd x y) ∧
    sigma (m₁ + m₂) ^ 2 * ((p.DZ : ℝ) / ((m₁ : ℝ) * m₂ * 10 ^ 24)) <
      (p.tZ : ℝ) / ((m₁ : ℝ) * m₂ * 10 ^ 24) := by
  simp only [rowOK, Bool.and_eq_true, decide_eq_true_eq, List.all_eq_true] at h
  obtain ⟨⟨⟨hm₁, hm₂⟩, hσ⟩, hall⟩ := h
  have hM : (0 : ℝ) < (m₁ : ℝ) * m₂ * 10 ^ 24 := by
    have : (0 : ℝ) < m₁ := by exact_mod_cast hm₁
    have : (0 : ℝ) < m₂ := by exact_mod_cast hm₂
    positivity
  have hd : ∀ v w : Vx (m₁ + 1) (m₂ + 1), descOK m₁ m₂ p
      (orb v, orb w, decide (ix1 v = ix1 w), decide (ix2 v = ix2 w)) = true :=
    fun v w => hall _ (mem_descs _ _ _ _)
  refine ⟨?_, fun c => long_pair (by omega) (famRow U W p) (fun v w hvw => ?_) c, ?_⟩
  · intro x hx y hy
    obtain ⟨v, _, rfl⟩ := Finset.mem_image.mp hx
    obtain ⟨w, _, rfl⟩ := Finset.mem_image.mp hy
    rw [sqd_famRow (by omega) (by omega) hU hW, div_le_div_iff_of_pos_right hM]
    have := hd v w
    simp only [descOK, Bool.and_eq_true, decide_eq_true_eq] at this
    exact_mod_cast this.1.1
  · rw [sqd_famRow (by omega) (by omega) hU hW, div_le_div_iff_of_pos_right hM]
    have := hd v w
    simp only [descOK, Bool.and_eq_true, decide_eq_true_eq, Bool.or_eq_true, Bool.not_eq_true'] at this
    rcases this.1.2 with h1 | h1
    · rw [adj_adjD hvw] at h1; exact absurd h1 (by decide)
    · exact_mod_cast h1
  · rw [← mul_div_assoc, div_lt_div_iff_of_pos_right hM]
    exact sigma_sq_mul_lt (by omega) hσ

theorem chunkOK_append_of {n m : ℕ} {ps qs : List Params} (h1 : chunkOK n ps = true)
    (hm : n + ps.length = m) (h2 : chunkOK m qs = true) : chunkOK n (ps ++ qs) = true := by
  rw [chunkOK_append, h1, hm, h2]
  rfl

/-- family_params.txt, rows n = 5…1000 -/
def table : List Params :=
  Rows00.rows ++ (Rows01.rows ++ (Rows02.rows ++ (Rows03.rows ++ (Rows04.rows ++ (Rows05.rows ++ (Rows06.rows ++ (Rows07.rows ++ (Rows08.rows ++ (Rows09.rows ++ (Rows10.rows ++ (Rows11.rows ++ (Rows12.rows ++ (Rows13.rows ++ (Rows14.rows ++ (Rows15.rows ++ (Rows16.rows ++ (Rows17.rows ++ (Rows18.rows ++ (Rows19.rows)))))))))))))))))))

theorem table_ok : chunkOK 5 table = true :=
  chunkOK_append_of Rows00.rows_ok (by rfl : 5 + Rows00.rows.length = 55)
    (chunkOK_append_of Rows01.rows_ok (by rfl : 55 + Rows01.rows.length = 105)
    (chunkOK_append_of Rows02.rows_ok (by rfl : 105 + Rows02.rows.length = 155)
    (chunkOK_append_of Rows03.rows_ok (by rfl : 155 + Rows03.rows.length = 205)
    (chunkOK_append_of Rows04.rows_ok (by rfl : 205 + Rows04.rows.length = 255)
    (chunkOK_append_of Rows05.rows_ok (by rfl : 255 + Rows05.rows.length = 305)
    (chunkOK_append_of Rows06.rows_ok (by rfl : 305 + Rows06.rows.length = 355)
    (chunkOK_append_of Rows07.rows_ok (by rfl : 355 + Rows07.rows.length = 405)
    (chunkOK_append_of Rows08.rows_ok (by rfl : 405 + Rows08.rows.length = 455)
    (chunkOK_append_of Rows09.rows_ok (by rfl : 455 + Rows09.rows.length = 505)
    (chunkOK_append_of Rows10.rows_ok (by rfl : 505 + Rows10.rows.length = 555)
    (chunkOK_append_of Rows11.rows_ok (by rfl : 555 + Rows11.rows.length = 605)
    (chunkOK_append_of Rows12.rows_ok (by rfl : 605 + Rows12.rows.length = 655)
    (chunkOK_append_of Rows13.rows_ok (by rfl : 655 + Rows13.rows.length = 705)
    (chunkOK_append_of Rows14.rows_ok (by rfl : 705 + Rows14.rows.length = 755)
    (chunkOK_append_of Rows15.rows_ok (by rfl : 755 + Rows15.rows.length = 805)
    (chunkOK_append_of Rows16.rows_ok (by rfl : 805 + Rows16.rows.length = 855)
    (chunkOK_append_of Rows17.rows_ok (by rfl : 855 + Rows17.rows.length = 905)
    (chunkOK_append_of Rows18.rows_ok (by rfl : 905 + Rows18.rows.length = 955)
    (Rows19.rows_ok)))))))))))))))))))

theorem table_length : table.length = 996 := by decide +kernel

/-- the parameters of the family member in dimension n (family_params.txt row n) -/
def params (n : ℕ) : Params := table.getD (n - 5) default

theorem params_ok {n : ℕ} (h5 : 5 ≤ n) (h1000 : n ≤ 1000) : rowOK (n / 2) (n - n / 2) (params n) = true := by
  have := chunkOK_getD table_ok (i := n - 5) (by rw [table_length]; omega)
  rwa [show 5 + (n - 5) = n by omega] at this

/-- the family member Vₙ ⊂ ℝ^⌊n/2⌋ ⊕ ℝ^⌈n/2⌉ (NOTE §2) with the parameters of row n -/
noncomputable def famSet (n : ℕ) (U : Fin (n / 2 + 1) → Fin (n / 2) → ℝ)
    (W : Fin (n - n / 2 + 1) → Fin (n - n / 2) → ℝ) : Finset (Fin (n / 2 + (n - n / 2)) → ℝ) :=
  Finset.univ.image (famRow U W (params n))

/-- **The family is a counterexample for every n = 5…1000** (NOTE §2–§3): no partition of Vₙ into n + 1
parts of diameter ≤ σₙ · diam Vₙ (n = ⌊n/2⌋ + ⌈n/2⌉). -/
theorem famSet_not_galeSplit {n : ℕ} (h5 : 5 ≤ n) (h1000 : n ≤ 1000)
    {U : Fin (n / 2 + 1) → Fin (n / 2) → ℝ} {W : Fin (n - n / 2 + 1) → Fin (n - n / 2) → ℝ}
    (hU : IsUnitSimplex U) (hW : IsUnitSimplex W) : ¬ GaleSplit (n / 2 + (n - n / 2)) (famSet n U W) :=
  rowOK_sound (params_ok h5 h1000) hU hW

/-- Vₙ has (⌊n/2⌋+1) + (⌈n/2⌉+1) + 3(⌊n/2⌋+1)(⌈n/2⌉+1) points -/
theorem card_famSet {n : ℕ} (h5 : 5 ≤ n) (h1000 : n ≤ 1000)
    {U : Fin (n / 2 + 1) → Fin (n / 2) → ℝ} {W : Fin (n - n / 2 + 1) → Fin (n - n / 2) → ℝ}
    (hU : IsUnitSimplex U) (hW : IsUnitSimplex W) :
    (famSet n U W).card = (n / 2 + 1) + (n - n / 2 + 1) + 3 * ((n / 2 + 1) * (n - n / 2 + 1)) := by
  rw [famSet, Finset.card_image_of_injective _ (famRow_injective (params_ok h5 h1000) hU hW),
    Finset.card_univ]
  exact card_Vx' _ _

end LeanGale
