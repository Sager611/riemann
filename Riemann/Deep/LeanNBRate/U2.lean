import Riemann.Deep.LeanNBRate.U2asym

/-!
# U2: the dual bound for zero configurations (Theorem R′ and its Corollary)
-/

open MeasureTheory Filter Topology Set Complex ComplexConjugate

noncomputable section

namespace Riemann.Deep.LeanNBRate

theorem ZeroConfig.mono {P : ExpSum} {T T' : ℝ} (hP : ZeroConfig P T) (h : T ≤ T') : ZeroConfig P T' :=
  ⟨hP.inj, hP.zero, fun j => (hP.height j).trans h, hP.mult_pos, hP.mult_le⟩

theorem card_groups_le (P : ExpSum) (T : Set ℝ) (s₀ : ℝ) :
    ((groups P T s₀).card : ℝ) ≤ (nodeSet P T).card := by
  have h := sum_groups_sum P T s₀ (fun _ => (1 : ℝ))
  simp only [Finset.sum_const, nsmul_eq_mul, mul_one] at h
  rw [← h, Finset.card_eq_sum_ones, Nat.cast_sum]
  refine Finset.sum_le_sum fun g hg => ?_
  have := (groups_nonempty hg).card_pos
  simp only [Nat.cast_one]
  exact_mod_cast this

theorem part_singleton_M (P : ExpSum) (γ : ℝ) : (P.part {γ} γ).M = nw P γ := by
  classical
  rw [ExpSum.part_M]; unfold nw
  refine Finset.sum_congr rfl fun j _ => ?_
  by_cases h : P.node j = γ <;> simp [h]

theorem inv_sq_sub_le {γ c Hs d : ℝ} (hd : 0 < d) (hdH : 2 * d ≤ Hs + d) (hγ : |γ| ≤ Hs - d)
    (hc : Hs + d ≤ |c|) : 1 / (γ - c) ^ 2 ≤ ((Hs + d) / (2 * d)) ^ 2 / c ^ 2 := by
  have hHd : 0 < Hs + d := by linarith
  have hc0 : 0 < |c| := by linarith
  have h1 : |c| - |γ| ≤ |γ - c| := by
    have := abs_sub_abs_le_abs_sub c γ
    rw [abs_sub_comm] at this
    linarith
  have h2 : 2 * d / (Hs + d) * |c| ≤ |c| - |γ| := by
    have e : |c| - (Hs - d) - 2 * d / (Hs + d) * |c| = (|c| - (Hs + d)) * (1 - 2 * d / (Hs + d)) := by
      field_simp; ring
    have h3 : 0 ≤ (|c| - (Hs + d)) * (1 - 2 * d / (Hs + d)) := by
      apply mul_nonneg (by linarith)
      rw [sub_nonneg, div_le_one hHd]; exact hdH
    linarith
  have h4 : 0 < 2 * d / (Hs + d) * |c| := by positivity
  have h5 : (2 * d / (Hs + d) * |c|) ^ 2 ≤ (γ - c) ^ 2 := by
    rw [← sq_abs (γ - c)]
    exact pow_le_pow_left₀ h4.le (h2.trans h1) 2
  have h6 : 0 < (γ - c) ^ 2 := lt_of_lt_of_le (by positivity) h5
  rw [div_le_div_iff₀ h6 (by rw [← sq_abs]; positivity)]
  have e2 : ((Hs + d) / (2 * d)) ^ 2 * (γ - c) ^ 2 ≥ ((Hs + d) / (2 * d)) ^ 2 * (2 * d / (Hs + d) * |c|) ^ 2 :=
    mul_le_mul_of_nonneg_left h5 (by positivity)
  have e3 : ((Hs + d) / (2 * d)) ^ 2 * (2 * d / (Hs + d) * |c|) ^ 2 = c ^ 2 := by
    rw [← sq_abs c]; field_simp
  linarith

theorem F1_ge_inv_sq {T : ℝ} (hT : 1 ≤ T) : 1 / T ^ 2 ≤ F1 T := by
  rw [F1_of_ge_one hT]
  have h1 : 1 ≤ Real.sqrt 2 := Real.one_le_sqrt.mpr (by norm_num)
  have h2 : 1 ≤ Real.sqrt (T + 1) := Real.one_le_sqrt.mpr (by linarith)
  apply div_le_div_of_nonneg_right _ (by positivity)
  nlinarith


/-- The θ-sum of Lemma E for a zero configuration. -/
theorem u2_theta {P : ExpSum} {T : ℝ} (hP : ZeroConfig P T) {Tlo Thi : Set ℝ} {s₀ q Y κe ε₁ Hs d CF1 mk : ℝ}
    {Ku J₀ : ℕ} (hq1 : 1 ≤ q) (hκe : 1 ≤ κe) (hε₁ : 0 ≤ ε₁) (hsm4 : 4 * ε₁ * Real.log (2 * κe) ≤ 1 / 2)
    (hd : 0 < d) (hdH : 2 * d ≤ Hs + d) (hHs : 0 ≤ Hs) (hmk : 0 ≤ mk)
    (hgapl : ∀ x ∈ nodeSet P Tlo, |x| ≤ Hs - d) (hcen : ∀ g ∈ groups P Thi s₀, Hs + d ≤ |gC g|)
    (hThi2 : ∀ γ ∈ nodeSet P Thi, 2 ≤ |γ|)
    (hgM : ∀ g ∈ groups P Thi s₀, (P.part ↑g (gC g)).M ≤ Ku ∧
      ((P.part ↑g (gC g)).M : ℝ) ≤ ε₁ * Real.log (|gC g| + 2))
    (hgW : ∀ g ∈ groups P Thi s₀, 0 ≤ gW g ∧ gW g ≤ 1 / q ^ 2 ∧ gW g ≤ 1)
    (hpow : ∀ M : ℕ, M ≤ Ku → (2 * κe) ^ (6 * M) ≤ (2 * κe) ^ 6 * q) (hKuY : (Ku : ℝ) + 1 ≤ Y)
    (hcard : (nodeSet P Tlo).card ≤ J₀)
    (hCF1 : ∑ γ ∈ nodeSet P Thi, (nw P γ : ℝ) * F1 |γ| ≤ CF1) :
    ∑ γ ∈ nodeSet P Tlo, ∑ g ∈ groups P Thi s₀,
        phiE κe Hs mk (q ^ 4) (P.part ↑g (gC g)).M (gW g) (1 / q ^ 2) ^ 2 / (γ - gC g) ^ 2 ≤
      (6 * (Hs + mk + 1) * (2 + mk + 4 * CB * (2 * κe) ^ 6 * q ^ 3 * Y ^ 3) / q ^ 2) ^ 2 *
        (J₀ * ((Hs + d) / (2 * d)) ^ 2 * CF1) := by
  set Φ := (6 * (Hs + mk + 1) * (2 + mk + 4 * CB * (2 * κe) ^ 6 * q ^ 3 * Y ^ 3) / q ^ 2) ^ 2
  set R := ((Hs + d) / (2 * d)) ^ 2
  have hΦ : 0 ≤ Φ := sq_nonneg _
  have hR : 0 ≤ R := sq_nonneg _
  have hW := groups_weight_le hP (Th := Thi) (s₀ := s₀) hκe hε₁ hsm4 (fun g hg => (hgM g hg).2)
    (fun g hg => (hgW g hg).2.2) hThi2
  have hterm : ∀ γ ∈ nodeSet P Tlo, ∀ g ∈ groups P Thi s₀,
      phiE κe Hs mk (q ^ 4) (P.part ↑g (gC g)).M (gW g) (1 / q ^ 2) ^ 2 / (γ - gC g) ^ 2 ≤
        Φ * R * ((2 * κe) ^ (4 * ((P.part ↑g (gC g)).M - 1)) / (gC g) ^ 2) := by
    intro γ hγ g hg
    obtain ⟨hW0, hWη, _⟩ := hgW g hg
    have hMK := (hgM g hg).1
    have hMY : ((P.part ↑g (gC g)).M : ℝ) + 1 ≤ Y := by
      have : ((P.part ↑g (gC g)).M : ℝ) ≤ Ku := by exact_mod_cast hMK
      linarith
    have h1 := phiE_sq_le (H₀ := Hs) (mk := mk) hκe hq1 hW0 hWη (hpow _ hMK) hMY hHs hmk
    have h2 := inv_sq_sub_le hd hdH (hgapl γ hγ) (hcen g hg)
    have h3 : 0 ≤ phiE κe Hs mk (q ^ 4) (P.part ↑g (gC g)).M (gW g) (1 / q ^ 2) ^ 2 := sq_nonneg _
    have h4 : 0 ≤ (2 * κe) ^ (4 * ((P.part ↑g (gC g)).M - 1)) := by positivity
    calc phiE κe Hs mk (q ^ 4) (P.part ↑g (gC g)).M (gW g) (1 / q ^ 2) ^ 2 / (γ - gC g) ^ 2
        = phiE κe Hs mk (q ^ 4) (P.part ↑g (gC g)).M (gW g) (1 / q ^ 2) ^ 2 * (1 / (γ - gC g) ^ 2) := by
          ring
      _ ≤ (Φ * (2 * κe) ^ (4 * ((P.part ↑g (gC g)).M - 1))) * (R / (gC g) ^ 2) := by
          apply mul_le_mul h1 h2 (by positivity) (by positivity)
      _ = Φ * R * ((2 * κe) ^ (4 * ((P.part ↑g (gC g)).M - 1)) / (gC g) ^ 2) := by ring
  calc ∑ γ ∈ nodeSet P Tlo, ∑ g ∈ groups P Thi s₀,
        phiE κe Hs mk (q ^ 4) (P.part ↑g (gC g)).M (gW g) (1 / q ^ 2) ^ 2 / (γ - gC g) ^ 2
      ≤ ∑ γ ∈ nodeSet P Tlo, ∑ g ∈ groups P Thi s₀,
          Φ * R * ((2 * κe) ^ (4 * ((P.part ↑g (gC g)).M - 1)) / (gC g) ^ 2) :=
        Finset.sum_le_sum fun γ hγ => Finset.sum_le_sum fun g hg => hterm γ hγ g hg
    _ = ∑ γ ∈ nodeSet P Tlo, Φ * R * ∑ g ∈ groups P Thi s₀,
          (2 * κe) ^ (4 * ((P.part ↑g (gC g)).M - 1)) / (gC g) ^ 2 := by
        refine Finset.sum_congr rfl fun γ _ => ?_
        rw [Finset.mul_sum]
    _ ≤ ∑ γ ∈ nodeSet P Tlo, Φ * R * CF1 :=
        Finset.sum_le_sum fun γ _ => mul_le_mul_of_nonneg_left (hW.trans hCF1) (by positivity)
    _ = (nodeSet P Tlo).card * (Φ * R * CF1) := by rw [Finset.sum_const, nsmul_eq_mul]
    _ ≤ J₀ * (Φ * R * CF1) := by
        have hCF0 : 0 ≤ CF1 := le_trans (Finset.sum_nonneg fun γ _ =>
          mul_nonneg (Nat.cast_nonneg _) (F1_nonneg _)) hCF1
        have : ((nodeSet P Tlo).card : ℝ) ≤ J₀ := by exact_mod_cast hcard
        exact mul_le_mul_of_nonneg_right this (by positivity)
    _ = Φ * (J₀ * R * CF1) := by ring

/-- The number of high groups is at most CF1·q⁵⁶. -/
theorem u2_card {P : ExpSum} {T q CF1 s₀ : ℝ} (hP : ZeroConfig P T) (hT1 : 1 ≤ T) (hTq : T ≤ q ^ 28)
    (Thi : Set ℝ) (hCF1 : ∑ γ ∈ nodeSet P Thi, (nw P γ : ℝ) * F1 |γ| ≤ CF1) :
    ((groups P Thi s₀).card : ℝ) ≤ CF1 * q ^ 56 := by
  have h1 := card_groups_le P Thi s₀
  have h2 := card_nodes_le hP hT1 Thi
  have hF := F1_ge_inv_sq hT1
  have hFp : 0 < F1 T := lt_of_lt_of_le (by positivity) hF
  have hS0 : 0 ≤ ∑ γ ∈ nodeSet P Thi, (nw P γ : ℝ) * F1 |γ| :=
    Finset.sum_nonneg fun γ _ => mul_nonneg (Nat.cast_nonneg _) (F1_nonneg _)
  have hCF0 : 0 ≤ CF1 := hS0.trans hCF1
  have h3 : (∑ γ ∈ nodeSet P Thi, (nw P γ : ℝ) * F1 |γ|) / F1 T ≤ CF1 * T ^ 2 := by
    rw [div_le_iff₀ hFp]
    have : CF1 * T ^ 2 * (1 / T ^ 2) = CF1 := by field_simp
    calc ∑ γ ∈ nodeSet P Thi, (nw P γ : ℝ) * F1 |γ| ≤ CF1 := hCF1
      _ = CF1 * T ^ 2 * (1 / T ^ 2) := this.symm
      _ ≤ CF1 * T ^ 2 * F1 T := mul_le_mul_of_nonneg_left hF (by positivity)
  have h4 : T ^ 2 ≤ q ^ 56 := by
    have : T ^ 2 ≤ (q ^ 28) ^ 2 := pow_le_pow_left₀ (by linarith) hTq 2
    calc T ^ 2 ≤ (q ^ 28) ^ 2 := this
      _ = q ^ 56 := by ring
  calc ((groups P Thi s₀).card : ℝ) ≤ (nodeSet P Thi).card := h1
    _ ≤ _ := h2
    _ ≤ CF1 * T ^ 2 := h3
    _ ≤ CF1 * q ^ 56 := mul_le_mul_of_nonneg_left h4 hCF0

/-- The Lemma D defect of the high part. -/
theorem u2_cdelta {κ q Y CF1 G w : ℝ} {M : ℕ} (hκ : 1 ≤ κ) (hq : 1 ≤ q) (hw0 : 0 ≤ w)
    (hw : w + 1 / q ^ 3 / 4 ≤ Y / q ^ 3) (hMY : (M : ℝ) + 1 + 8 * A₀ ≤ Y) (hA : 1 + 16 * A₀ ≤ 2 * Y)
    (hpow : (2 * κ) ^ (6 * M) ≤ (2 * κ) ^ 6 * q) (hG0 : 0 ≤ G) (hG : G ≤ CF1 * q ^ 56) :
    ((1 + (1 / q ^ 4)⁻¹) * G + 1) * (2 * Real.pi) * tailC κ M w (q ^ 4) (1 / q ^ 3) 40 ≤
      (2 * CF1 + 1) * (2 * Real.pi) * KT * (2 * κ) ^ 6 * Y ^ 84 / q ^ 16 := by
  have hq0 : 0 < q := by linarith
  have ht := tailC_bound hκ hq hw0 hw hMY hA
  have hKT := KT_nonneg
  have hY0 : 0 ≤ Y := by have := A₀_pos; linarith
  have hCF0 : 0 ≤ CF1 := by
    have : 0 ≤ CF1 * q ^ 56 := hG0.trans hG
    exact nonneg_of_mul_nonneg_left this (by positivity)
  have ht0 : 0 ≤ tailC κ M w (q ^ 4) (1 / q ^ 3) 40 := by
    unfold tailC; have := B₀_pos; have := A₀_pos; positivity
  have h1 : (1 + (1 / q ^ 4)⁻¹) * G + 1 ≤ (2 * CF1 + 1) * q ^ 60 := by
    rw [one_div, inv_inv]
    have hq4 : 1 ≤ q ^ 4 := one_le_pow₀ hq
    have hq60 : 1 ≤ q ^ 60 := one_le_pow₀ hq
    have e : (1 + q ^ 4) * (CF1 * q ^ 56) ≤ 2 * CF1 * q ^ 60 := by
      have : (1 + q ^ 4) * (CF1 * q ^ 56) = CF1 * q ^ 56 + CF1 * q ^ 60 := by ring
      rw [this]
      have : CF1 * q ^ 56 ≤ CF1 * q ^ 60 := mul_le_mul_of_nonneg_left
        (pow_le_pow_right₀ hq (by norm_num)) hCF0
      linarith
    have : (1 + q ^ 4) * G ≤ (1 + q ^ 4) * (CF1 * q ^ 56) := mul_le_mul_of_nonneg_left hG (by positivity)
    nlinarith
  have h2 : tailC κ M w (q ^ 4) (1 / q ^ 3) 40 ≤ KT * ((2 * κ) ^ 6 * q) * Y ^ 84 / q ^ 77 := by
    refine ht.trans ?_
    gcongr
  have hL0 : 0 ≤ (1 + (1 / q ^ 4)⁻¹) * G + 1 := by positivity
  calc ((1 + (1 / q ^ 4)⁻¹) * G + 1) * (2 * Real.pi) * tailC κ M w (q ^ 4) (1 / q ^ 3) 40
      ≤ ((2 * CF1 + 1) * q ^ 60) * (2 * Real.pi) * (KT * ((2 * κ) ^ 6 * q) * Y ^ 84 / q ^ 77) := by
        gcongr
    _ = (2 * CF1 + 1) * (2 * Real.pi) * KT * (2 * κ) ^ 6 * Y ^ 84 / q ^ 16 := by
        field_simp


end Riemann.Deep.LeanNBRate
