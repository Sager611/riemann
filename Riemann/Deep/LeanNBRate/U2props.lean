import Riemann.Deep.LeanNBRate.U2

/-!
# U2: the three per-configuration estimates as predicates, and the final arithmetic
-/

open MeasureTheory Filter Topology Set Complex ComplexConjugate

noncomputable section

namespace Riemann.Deep.LeanNBRate

/-- The conclusion of highR_parts for a given constant. -/
def HRprop (κ : ℝ) : Prop :=
    ∀ (P : ExpSum) (T : Set ℝ) (s₀ ℓ H η η' wmax : ℝ) (N Mmax : ℕ),
      0 < s₀ → 0 < ℓ → s₀ < 2 * H → 0 < η → 0 < η' → 1 ≤ N →
      (∀ j, P.node j ∈ T → H ≤ |P.node j|) → (∀ j, P.node j ∈ T → 1 ≤ P.mult j) →
      (∀ g ∈ groups P T s₀, (P.part ↑g (gC g)).M ≤ Mmax ∧ gW g ≤ wmax) →
      4 * A₀ / η * ((Mmax : ℝ) + 1) ^ 2 ≤ ℓ →
      (Mmax : ℝ) * (4 * A₀ / (s₀ / 4)) * ((Mmax : ℝ) + 1) ^ 2 ≤ ℓ →
      smallB κ Mmax wmax η ≤ H ^ 2 / 2 →
      ‖(P.part T 0).eval 0‖ ^ 2 ≤
        ((∑ g ∈ groups P T s₀, 16 * (2 * κ) ^ (4 * ((P.part ↑g (gC g)).M - 1)) / (gC g) ^ 2) / ℓ) *
          (∑ g ∈ groups P T s₀, ∫ x : ℝ, ‖grpQ (P.part ↑g (gC g)) (gC g) ℓ x‖ ^ 2) ∧
      (1 - ((1 + η'⁻¹) * ((groups P T s₀).card : ℝ) + 1) * (2 * Real.pi) *
          tailC κ Mmax wmax ℓ s₀ N) *
          (∑ g ∈ groups P T s₀, ∫ x : ℝ, ‖grpQ (P.part ↑g (gC g)) (gC g) ℓ x‖ ^ 2) ≤
        (1 + η') * ∫ u : ℝ, ‖((P.part T 0).dmod 0).eval u * psiC ℓ u‖ ^ 2

/-- The conclusion of lemmaE for given constants. -/
def LEprop (κ c : ℝ) : Prop :=
    ∀ (P : ExpSum) (Tlo Thi : Set ℝ) (s₀ ℓ η H₀ : ℝ) (M₀ : ℕ), 1 ≤ ℓ → 0 < η → 0 ≤ H₀ →
      (∀ γ ∈ nodeSet P Tlo, |γ| ≤ H₀) →
      (∀ γ ∈ nodeSet P Tlo, 1 ≤ (P.part {γ} γ).M ∧ (P.part {γ} γ).M ≤ M₀) →
      (∀ g ∈ groups P Thi s₀, 1 ≤ (P.part ↑g (gC g)).M ∧ 1 ≤ |gC g| ∧ 0 ≤ gW g ∧
        (∀ j, (P.part ↑g (gC g)).node j ∈ Icc (-gW g) (gW g)) ∧
        4 * A₀ / η * (((P.part ↑g (gC g)).M : ℝ) + 1) ^ 2 ≤ ℓ ∧
        smallB κ (P.part ↑g (gC g)).M (gW g) η ≤ (gC g) ^ 2 / 2) →
      (∀ γ ∈ nodeSet P Tlo, ∀ g ∈ groups P Thi s₀, γ ≠ gC g) →
      ‖oubil ℓ (P.part Tlo 0).eval ((P.part Tlo 0).dmod 0).eval (P.part Thi 0).eval
          ((P.part Thi 0).dmod 0).eval‖ ≤
        2 * Real.sqrt (2 * (2 * c) ^ (2 * (M₀ - 1)) / ℓ) *
          Real.sqrt (∑ γ ∈ nodeSet P Tlo, ∑ g ∈ groups P Thi s₀,
            phiE κ H₀ (markovC c M₀) ℓ (P.part ↑g (gC g)).M (gW g) η ^ 2 / (γ - gC g) ^ 2) *
          Real.sqrt (∑ γ ∈ nodeSet P Tlo, (γ ^ 2 + 1 / 4) * ∫ u in (0 : ℝ)..ℓ, ‖(P.part {γ} γ).eval u‖ ^ 2) *
          Real.sqrt (∑ g ∈ groups P Thi s₀, ∫ x : ℝ, ‖grpQ (P.part ↑g (gC g)) (gC g) ℓ x‖ ^ 2)

/-- The conclusion of claimF_core for a given constant. -/
def CFprop (c : ℝ) : Prop :=
    ∀ (P : ExpSum) (T : Set ℝ) (ℓ g₀ H₀ : ℝ) (M₀ : ℕ), 1 ≤ ℓ → 0 < g₀ → 0 ≤ H₀ →
      (∀ γ ∈ nodeSet P T, ∀ γ' ∈ nodeSet P T, γ ≠ γ' → g₀ ≤ |γ - γ'|) →
      (∀ γ ∈ nodeSet P T, |γ| ≤ H₀) →
      (∀ γ ∈ nodeSet P T, 1 ≤ (P.part {γ} γ).M ∧ (P.part {γ} γ).M ≤ M₀) →
      (1 - claimFC c M₀ H₀ g₀ (nodeSet P T).card / ℓ) *
          (∑ γ ∈ nodeSet P T, (γ ^ 2 + 1 / 4) * ∫ u in (0 : ℝ)..ℓ, ‖(P.part {γ} γ).eval u‖ ^ 2) ≤
        (oubil ℓ (P.part T 0).eval ((P.part T 0).dmod 0).eval (P.part T 0).eval
          ((P.part T 0).dmod 0).eval).re ∧
      ‖(P.part T 0).eval 0‖ ^ 2 ≤
        (∑ γ ∈ nodeSet P T, ((P.part {γ} γ).M : ℝ) ^ 2 / (γ ^ 2 + 1 / 4)) *
          ((∑ γ ∈ nodeSet P T, (γ ^ 2 + 1 / 4) * ∫ u in (0 : ℝ)..ℓ, ‖(P.part {γ} γ).eval u‖ ^ 2) / ℓ)

theorem highR_parts' (hT : Nazarov1993_Turan) : ∃ κ : ℝ, 1 ≤ κ ∧ HRprop κ := highR_parts hT

theorem lemmaE' (hT : Nazarov1993_Turan) : ∃ κ c : ℝ, 1 ≤ κ ∧ 1 ≤ c ∧ LEprop κ c := lemmaE hT

theorem claimF_core' (hT : Nazarov1993_Turan) : ∃ c : ℝ, 1 ≤ c ∧ CFprop c := claimF_core hT

theorem theta_final {σ K H m C κ q Y δ₀ SS : ℝ} (hσ : 0 ≤ σ) (hK : 0 ≤ K) (hH : 0 ≤ H) (hm : 0 ≤ m)
    (hC : 0 ≤ C) (hκ : 0 ≤ κ) (hq1 : 1 ≤ q) (hY1 : 1 ≤ Y) (h3 : Y ^ 3 ≤ δ₀ * q)
    (hSS : SS ≤ (6 * (H + m + 1) * (2 + m + 4 * C * (2 * κ) ^ 6 * q ^ 3 * Y ^ 3) / q ^ 2) ^ 2 * K) :
    2 * Real.sqrt (σ / q ^ 4) * Real.sqrt SS ≤
      (12 * Real.sqrt (σ * K) * (H + m + 1) * (2 + m + 4 * C * (2 * κ) ^ 6)) * δ₀ := by
  have hqp : 0 < q := by linarith
  set Φ₁ := 6 * (H + m + 1) * (2 + m + 4 * C * (2 * κ) ^ 6 * q ^ 3 * Y ^ 3) / q ^ 2 with hΦdef
  have hΦ₁ : 0 ≤ Φ₁ := by positivity
  have hs1 : Real.sqrt SS ≤ Φ₁ * Real.sqrt K := by
    calc Real.sqrt SS ≤ Real.sqrt (Φ₁ ^ 2 * K) := Real.sqrt_le_sqrt hSS
      _ = Φ₁ * Real.sqrt K := by rw [Real.sqrt_mul (sq_nonneg _), Real.sqrt_sq hΦ₁]
  have hs2 : Real.sqrt (σ / q ^ 4) = Real.sqrt σ / q ^ 2 := by
    rw [Real.sqrt_div' _ (by positivity), show q ^ 4 = (q ^ 2) ^ 2 by ring, Real.sqrt_sq (by positivity)]
  have hq3Y : 1 ≤ q ^ 3 * Y ^ 3 := one_le_mul_of_one_le_of_one_le (one_le_pow₀ hq1) (one_le_pow₀ hY1)
  have hsK : Real.sqrt σ * Real.sqrt K = Real.sqrt (σ * K) := (Real.sqrt_mul hσ _).symm
  have hK0 : 0 ≤ Real.sqrt (σ * K) := Real.sqrt_nonneg _
  have hA : 0 ≤ 4 * C * (2 * κ) ^ 6 := by positivity
  have hB : 2 + m + 4 * C * (2 * κ) ^ 6 * q ^ 3 * Y ^ 3 ≤ (2 + m + 4 * C * (2 * κ) ^ 6) * (q ^ 3 * Y ^ 3) := by
    nlinarith
  have hKθ : 0 ≤ 12 * Real.sqrt (σ * K) * (H + m + 1) * (2 + m + 4 * C * (2 * κ) ^ 6) := by positivity
  calc 2 * Real.sqrt (σ / q ^ 4) * Real.sqrt SS = 2 * (Real.sqrt σ / q ^ 2) * Real.sqrt SS := by rw [hs2]
    _ ≤ 2 * (Real.sqrt σ / q ^ 2) * (Φ₁ * Real.sqrt K) := by gcongr
    _ = 12 * Real.sqrt (σ * K) * (H + m + 1) * (2 + m + 4 * C * (2 * κ) ^ 6 * q ^ 3 * Y ^ 3) / q ^ 4 := by
        rw [← hsK, hΦdef]; field_simp; ring
    _ ≤ 12 * Real.sqrt (σ * K) * (H + m + 1) * ((2 + m + 4 * C * (2 * κ) ^ 6) * (q ^ 3 * Y ^ 3)) / q ^ 4 := by
        gcongr
    _ = (12 * Real.sqrt (σ * K) * (H + m + 1) * (2 + m + 4 * C * (2 * κ) ^ 6)) * (Y ^ 3 / q) := by
        first | (field_simp; done) | (field_simp; ring)
    _ ≤ (12 * Real.sqrt (σ * K) * (H + m + 1) * (2 + m + 4 * C * (2 * κ) ^ 6)) * δ₀ := by
        apply mul_le_mul_of_nonneg_left _ hKθ
        rw [div_le_iff₀ hqp]; exact h3

theorem b_final {Cδ η' δ₀ Kb : ℝ} (hC0 : 0 ≤ Cδ) (hη' : 0 < η') (hηδ : η' ≤ δ₀) (hC : Cδ ≤ Kb * δ₀) :
    0 ≤ 1 - (1 - Cδ) / (1 + η') ∧ 1 - (1 - Cδ) / (1 + η') ≤ δ₀ + Kb * δ₀ := by
  have hden : 0 < 1 + η' := by linarith
  have e : 1 - (1 - Cδ) / (1 + η') = (η' + Cδ) / (1 + η') := by
    first | (field_simp; done) | (field_simp; ring)
  rw [e]
  refine ⟨by positivity, ?_⟩
  rw [div_le_iff₀ hden]
  have h1 : η' + Cδ ≤ δ₀ + Kb * δ₀ := by linarith
  have h2 : 0 ≤ (δ₀ + Kb * δ₀) * η' := by
    have : 0 ≤ δ₀ + Kb * δ₀ := le_trans (by positivity) h1
    positivity
  nlinarith

theorem cdelta_final {Kb Y q δ₀ : ℝ} (hKb : 0 ≤ Kb) (hq1 : 1 ≤ q) (h84 : Y ^ 84 ≤ δ₀ * q) (hδ0 : 0 ≤ δ₀) :
    Kb * Y ^ 84 / q ^ 16 ≤ Kb * δ₀ := by
  have hqp : 0 < q := by linarith
  have hq16 : q ≤ q ^ 16 := by
    calc q = q ^ 1 := (pow_one q).symm
      _ ≤ q ^ 16 := pow_le_pow_right₀ hq1 (by norm_num)
  rw [mul_div_assoc]
  apply mul_le_mul_of_nonneg_left _ hKb
  rw [div_le_iff₀ (by positivity)]
  calc Y ^ 84 ≤ δ₀ * q := h84
    _ ≤ δ₀ * q ^ 16 := mul_le_mul_of_nonneg_left hq16 hδ0

theorem smallB_half {κ' κ q Y δ₀ : ℝ} (h1 : 1 ≤ κ') (h2 : κ' ≤ κ) (hq1 : 1 ≤ q) (h6 : Y ^ 6 ≤ δ₀ * q)
    (hδ0 : 0 ≤ δ₀) (hδB : δ₀ ≤ 1 / (128 * CB ^ 2 * (2 * κ) ^ 6 + 1)) :
    64 * CB ^ 2 * (2 * κ') ^ 6 * Y ^ 6 / q ^ 3 ≤ 1 / 2 := by
  have hqp : 0 < q := by linarith
  have hCB := CB_nonneg
  have hk6 : (2 * κ') ^ 6 ≤ (2 * κ) ^ 6 := pow_le_pow_left₀ (by linarith) (by linarith) 6
  have hq3 : q ≤ q ^ 3 := by
    calc q = q ^ 1 := (pow_one q).symm
      _ ≤ q ^ 3 := pow_le_pow_right₀ hq1 (by norm_num)
  have hD : 0 < 128 * CB ^ 2 * (2 * κ) ^ 6 + 1 := by positivity
  have hX : 0 ≤ 64 * CB ^ 2 * (2 * κ) ^ 6 := by positivity
  calc 64 * CB ^ 2 * (2 * κ') ^ 6 * Y ^ 6 / q ^ 3 ≤ 64 * CB ^ 2 * (2 * κ) ^ 6 * (δ₀ * q) / q ^ 3 := by
        gcongr
    _ ≤ 64 * CB ^ 2 * (2 * κ) ^ 6 * δ₀ := by
        rw [div_le_iff₀ (by positivity)]
        calc 64 * CB ^ 2 * (2 * κ) ^ 6 * (δ₀ * q) = (64 * CB ^ 2 * (2 * κ) ^ 6 * δ₀) * q := by ring
          _ ≤ (64 * CB ^ 2 * (2 * κ) ^ 6 * δ₀) * q ^ 3 :=
              mul_le_mul_of_nonneg_left hq3 (by positivity)
    _ ≤ 64 * CB ^ 2 * (2 * κ) ^ 6 * (1 / (128 * CB ^ 2 * (2 * κ) ^ 6 + 1)) :=
        mul_le_mul_of_nonneg_left hδB hX
    _ ≤ 1 / 2 := by
        rw [mul_one_div, div_le_iff₀ hD]; nlinarith

theorem final_arith {p2 Slo Shi Z ε Q ℓ x x₀ : ℝ} (hε : 0 < ε) (hZ : 0 ≤ Z) (hQ : 0 ≤ Q) (hℓ : 0 < ℓ)
    (hx0 : 0 ≤ x) (hxx : x ≤ x₀) (hx₀ : x₀ = ε / (2 * (Z + ε))) (hSlo : Slo ≤ Z) (hShi : Shi ≤ ε / 3)
    (h : p2 ≤ (Slo + Shi) / ℓ * Q / (1 - x)) : p2 ≤ (Z + ε) / ℓ * Q := by
  have hZε : 0 < Z + ε := by linarith
  have hx₀1 : x₀ ≤ 1 / 2 := by
    rw [hx₀, div_le_iff₀ (by positivity)]; linarith
  have hx1 : 0 < 1 - x := by linarith
  have hxZ : x * (Z + ε) ≤ ε / 2 := by
    calc x * (Z + ε) ≤ x₀ * (Z + ε) := mul_le_mul_of_nonneg_right hxx hZε.le
      _ = ε / 2 := by rw [hx₀]; field_simp
  have hQq : 0 ≤ Q / ℓ := div_nonneg hQ hℓ.le
  refine h.trans ?_
  rw [div_le_iff₀ hx1]
  have e2 : (Slo + Shi) / ℓ * Q = (Slo + Shi) * (Q / ℓ) := by ring
  have e3 : (Z + ε) / ℓ * Q * (1 - x) = ((Z + ε) * (1 - x)) * (Q / ℓ) := by ring
  rw [e2, e3]
  apply mul_le_mul_of_nonneg_right _ hQq
  nlinarith

end Riemann.Deep.LeanNBRate
