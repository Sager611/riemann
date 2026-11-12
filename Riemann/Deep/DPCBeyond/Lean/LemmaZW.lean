import Riemann.Deep.DPCBeyond.Lean.LemmaZ2

/-!
# dpcbeyond: Lemma Z's proved output wired into M3's first inequality and 6e′ (a)

m3_first_of_lemmaZ: M3Inputs' first inequality ∫|R_T|^q ≤ K^q(∫(1+W)^{qr})^{1/r}(∫((1+On)^{k−1})^{qr′})^{1−1/r} as a
THEOREM from Lemma Z's named inputs (ZPoint1 at y = 0, every t ∈ [T, 2T]) plus three NAMED integral facts (the integral
over [T, 2T] is monotone and positively homogeneous on |R_T|^q and (1 + W)^q((1 + On)^{k−1})^q; Hölder with r, r′ = r/(r − 1)).
hairpin_pointwise: 6e′ (a), |R_T R*|(x − iy) ≤ K²(X^y + Off)(1 + Off)((1 + On)^k)² as a THEOREM from Lemma Z's named
inputs (ZPoint1, ZPoint2 at the point) and W_split.
-/

namespace DPCBLean

open Finset DPCLean

variable {ex lg : ℝ → ℝ}

/-- **M3's first inequality from Lemma Z** (y = 0): Lemma Z's named inputs (ZPoint1 at y = 0) at every t ∈ [T, 2T]
and two NAMED integral facts — (hcmp) the integral over [T, 2T] is monotone and positively homogeneous on the two
integrable functions |R_T|^q and (1 + W)^q((1 + On)^{k−1})^q (Rq = ∫_T^{2T}|R_T|^q, IA = ∫_T^{2T}(1 + W)^q((1 + On)^{k−1})^q),
(hHol) Hölder IA ≤ B (B = (∫(1+W)^{qr})^{1/r}(∫((1+On)^{k−1})^{qr′})^{1−1/r} in M3InputsW) — give Rq ≤ K^q B,
K = NT K_pk C_I k + 1. -/
theorem m3_first_of_lemmaZ (h : ExpLog ex lg) {ι : Type*} [DecidableEq ι] {T X q CI Kpk Rq IA B : ℝ}
    {k NT : ℕ} (hX : 0 < X) (hq : 0 ≤ q) (hCI : 0 ≤ CI) (hKpk : 0 ≤ Kpk) (hk : 1 ≤ k)
    (S : ℝ → Finset ι) (d : ι → ℝ) (κ : ℝ → ι → ℝ) (hκ : ∀ t x, 0 ≤ κ t x) (Rabs : ℝ → ℝ)
    (hz : ∀ t, T ≤ t → t ≤ 2 * T → ZPoint1 ex lg X 0 CI Kpk k NT (S t) d (κ t) (Rabs t))
    (hcmp : ∀ c : ℝ, 0 ≤ c → (∀ t, T ≤ t → t ≤ 2 * T → rp ex lg (Rabs t) q ≤
      c * (rp ex lg (1 + Wsum ex lg X 0 (S t) d (κ t)) q *
        rp ex lg ((1 + Onsum (S t) (κ t)) ^ (k - 1)) q)) → Rq ≤ c * IA)
    (hHol : IA ≤ B) : Rq ≤ rp ex lg (NT * (Kpk * CI * k) + 1) q * B := by
  have hK : 0 < (NT : ℝ) * (Kpk * CI * k) + 1 := by
    have h0 : 0 ≤ (NT : ℝ) * (Kpk * CI * k) :=
      mul_nonneg (Nat.cast_nonneg _) (mul_nonneg (mul_nonneg hKpk hCI) (Nat.cast_nonneg _))
    linarith only [h0]
  have hpt : ∀ t, T ≤ t → t ≤ 2 * T → rp ex lg (Rabs t) q ≤ rp ex lg (NT * (Kpk * CI * k) + 1) q *
      (rp ex lg (1 + Wsum ex lg X 0 (S t) d (κ t)) q *
        rp ex lg ((1 + Onsum (S t) (κ t)) ^ (k - 1)) q) := by
    intro t h1 h2
    have hW : 0 ≤ Wsum ex lg X 0 (S t) d (κ t) :=
      sum_nonneg (fun x _ => mul_nonneg (rp_nonneg h _ _) (hκ t x))
    have hOn : 0 ≤ Onsum (S t) (κ t) := sum_nonneg (fun x _ => hκ t x)
    have hA : 0 < 1 + Wsum ex lg X 0 (S t) d (κ t) := by linarith only [hW]
    have hB : 0 < (1 + Onsum (S t) (κ t)) ^ (k - 1) := pow_pos (by linarith only [hOn]) _
    have hZ := lemmaZ_point1 h hX le_rfl hCI hKpk hk (hκ t) (hz t h1 h2)
    calc rp ex lg (Rabs t) q ≤ rp ex lg ((NT * (Kpk * CI * k) + 1) *
          ((1 + Wsum ex lg X 0 (S t) d (κ t)) * (1 + Onsum (S t) (κ t)) ^ (k - 1))) q :=
          rp_le_rp h hZ hq
      _ = rp ex lg (NT * (Kpk * CI * k) + 1) q * (rp ex lg (1 + Wsum ex lg X 0 (S t) d (κ t)) q *
          rp ex lg ((1 + Onsum (S t) (κ t)) ^ (k - 1)) q) := by
          rw [rp_mul h hK (mul_pos hA hB), rp_mul h hA hB]
  exact le_trans (hcmp _ (rp_nonneg h _ _) hpt) (mul_le_mul_of_nonneg_left hHol (rp_nonneg h _ _))

/-- **6e′ (a) from Lemma Z**: at a cut point x and 0 ≤ y, Lemma Z's named inputs for R_T(x − iy) (ZPoint1) and
R*(x − iy) (ZPoint2) give |R_T R*|(x − iy) ≤ K²(X^y + Off)(1 + Off)((1 + On)^k)², K = NT K_pk C_I k + 1. -/
theorem hairpin_pointwise (h : ExpLog ex lg) {ι : Type*} [DecidableEq ι] {X y CI Kpk : ℝ} {k NT : ℕ}
    {S : Finset ι} {d κ : ι → ℝ} {R1 R2 : ℝ} (hX : 1 ≤ X) (hy : 0 ≤ y) (hCI : 0 ≤ CI)
    (hKpk : 0 ≤ Kpk) (hk : 1 ≤ k) (hκ : ∀ x, 0 ≤ κ x) (hR1 : 0 ≤ R1)
    (h1 : ZPoint1 ex lg X y CI Kpk k NT S d κ R1) (h2 : ZPoint2 ex lg X y CI Kpk k NT S d κ R2) :
    R1 * R2 ≤ (NT * (Kpk * CI * k) + 1) ^ 2 * ((rp ex lg X y + Offsum ex lg X S d κ) *
      (1 + Offsum ex lg X S d κ) * ((1 + Onsum S κ) ^ k) ^ 2) := by
  have hX0 : 0 < X := by linarith only [hX]
  have hb1 := lemmaZ_point1 h hX0 hy hCI hKpk hk hκ h1
  have hb2 := lemmaZ_point2 h hX hy hCI hKpk hk hκ h2
  have hw1 := W_split h S d κ hκ hX hy
  have hw0 := W_split h S d κ hκ hX le_rfl
  rw [neg_zero, rp_zero h hX0, one_mul] at hw0
  obtain ⟨K, hK⟩ : ∃ K : ℝ, K = NT * (Kpk * CI * k) + 1 := ⟨_, rfl⟩
  obtain ⟨On, hOn⟩ : ∃ x : ℝ, x = Onsum S κ := ⟨_, rfl⟩
  obtain ⟨Off, hOff⟩ : ∃ x : ℝ, x = Offsum ex lg X S d κ := ⟨_, rfl⟩
  obtain ⟨Xy, hXy⟩ : ∃ x : ℝ, x = rp ex lg X y := ⟨_, rfl⟩
  obtain ⟨Xm, hXm⟩ : ∃ x : ℝ, x = rp ex lg X (-y) := ⟨_, rfl⟩
  rw [← hK] at hb1 hb2 ⊢
  rw [← hOn] at hb1 hb2 hw1 hw0 ⊢
  rw [← hOff] at hw1 hw0 ⊢
  rw [← hXy] at hb2 ⊢
  rw [← hXm] at hw1
  have hK0 : 0 ≤ K := by
    rw [hK]
    have h0 : 0 ≤ (NT : ℝ) * (Kpk * CI * k) :=
      mul_nonneg (Nat.cast_nonneg _) (mul_nonneg (mul_nonneg hKpk hCI) (Nat.cast_nonneg _))
    linarith only [h0]
  have hOn0 : 0 ≤ On := by rw [hOn]; exact sum_nonneg (fun x _ => hκ x)
  have hOff0 : 0 ≤ Off := by
    rw [hOff]; exact sum_nonneg (fun x _ => mul_nonneg (rp_nonneg h _ _) (hκ x))
  have hXy0 : 0 ≤ Xy := by rw [hXy]; exact rp_nonneg h _ _
  have hXyXm : Xy * Xm = 1 := by rw [hXy, hXm, mul_comm]; exact rp_neg h hX0 y
  have hW1 : 0 ≤ Wsum ex lg X y S d κ := sum_nonneg (fun x _ => mul_nonneg (rp_nonneg h _ _) (hκ x))
  have hW0 : 0 ≤ Wsum ex lg X 0 S d κ := sum_nonneg (fun x _ => mul_nonneg (rp_nonneg h _ _) (hκ x))
  obtain ⟨Q, hQ⟩ : ∃ x : ℝ, x = (1 + On) ^ (k - 1) := ⟨_, rfl⟩
  rw [← hQ] at hb1 hb2
  have hQ0 : 0 ≤ Q := by rw [hQ]; exact pow_nonneg (by linarith only [hOn0]) _
  have hQk : Q * (1 + On) = (1 + On) ^ k := by
    rw [hQ, ← pow_succ, Nat.sub_add_cancel hk]
  -- R1 R2 ≤ (K(1+W₁)Q)(K Xy (1+W₀) Q)
  have b2n : 0 ≤ K * (Xy * ((1 + Wsum ex lg X 0 S d κ) * Q)) := by
    have : 0 ≤ 1 + Wsum ex lg X 0 S d κ := by linarith only [hW0]
    positivity
  have s1 : R1 * R2 ≤ (K * ((1 + Wsum ex lg X y S d κ) * Q)) * (K * (Xy * ((1 + Wsum ex lg X 0 S d κ) * Q))) :=
    le_trans (mul_le_mul_of_nonneg_left hb2 hR1) (mul_le_mul_of_nonneg_right hb1 b2n)
  have A1 : 1 + Wsum ex lg X y S d κ ≤ 1 + On + Xm * Off := by linarith only [hw1]
  have A0 : 1 + Wsum ex lg X 0 S d κ ≤ 1 + On + Off := by linarith only [hw0]
  have hKQ : 0 ≤ K * Q := mul_nonneg hK0 hQ0
  have s2 : (K * ((1 + Wsum ex lg X y S d κ) * Q)) * (K * (Xy * ((1 + Wsum ex lg X 0 S d κ) * Q))) ≤
      (K * Q) ^ 2 * ((Xy * (1 + On + Xm * Off)) * (1 + On + Off)) := by
    have e1 : (K * ((1 + Wsum ex lg X y S d κ) * Q)) * (K * (Xy * ((1 + Wsum ex lg X 0 S d κ) * Q))) =
        (K * Q) ^ 2 * ((Xy * (1 + Wsum ex lg X y S d κ)) * (1 + Wsum ex lg X 0 S d κ)) := by ring
    rw [e1]
    apply mul_le_mul_of_nonneg_left _ (sq_nonneg _)
    apply mul_le_mul (mul_le_mul_of_nonneg_left A1 hXy0) A0 (by linarith only [hW0])
    have : 0 ≤ 1 + On + Xm * Off := by linarith only [hw1, hW1]
    exact mul_nonneg hXy0 this
  have B1 : Xy * (1 + On + Xm * Off) = Xy * (1 + On) + Off := by
    rw [mul_add, ← mul_assoc, hXyXm, one_mul]
  have C1 : Xy * (1 + On) + Off ≤ (Xy + Off) * (1 + On) := by
    have e : (Xy + Off) * (1 + On) = Xy * (1 + On) + Off + Off * On := by ring
    have h0 := mul_nonneg hOff0 hOn0
    linarith only [e, h0]
  have C2 : 1 + On + Off ≤ (1 + Off) * (1 + On) := by
    have e : (1 + Off) * (1 + On) = 1 + On + Off + Off * On := by ring
    have h0 := mul_nonneg hOff0 hOn0
    linarith only [e, h0]
  have s3 : (Xy * (1 + On + Xm * Off)) * (1 + On + Off) ≤ ((Xy + Off) * (1 + On)) * ((1 + Off) * (1 + On)) := by
    rw [B1]
    exact mul_le_mul C1 C2 (by linarith only [hOn0, hOff0]) (mul_nonneg (by linarith only [hXy0, hOff0])
      (by linarith only [hOn0]))
  have e4 : (K * Q) ^ 2 * (((Xy + Off) * (1 + On)) * ((1 + Off) * (1 + On))) =
      K ^ 2 * ((Xy + Off) * (1 + Off) * (Q * (1 + On)) ^ 2) := by ring
  rw [hQk] at e4
  have s4 := mul_le_mul_of_nonneg_left s3 (sq_nonneg (K * Q))
  linarith only [s1, s2, s4, e4]

end DPCBLean

