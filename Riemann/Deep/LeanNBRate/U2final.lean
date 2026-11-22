import Riemann.Deep.LeanNBRate.U2main

/-!
# U2: the dual bound for all large ℓ
-/

open MeasureTheory Filter Topology Set Complex ComplexConjugate

noncomputable section

namespace Riemann.Deep.LeanNBRate

/-- **U2** (Theorem R′ with its Corollary, [N] §§4–5): under RH, for every ε > 0 and all large ℓ,
every zero configuration of height at most ℓ⁷ satisfies |p(0)|² ≤ (S + ε)/ℓ · Q_ℓ(p). -/
theorem u2 (hT : Nazarov1993_Turan) (hGG : GoldstonGonek2007_Thm1) (hc : Titchmarsh_Thm9_2)
    (hRH : RiemannHypothesis) : ∀ ε > 0, ∃ ℓ₀ : ℝ, 1 ≤ ℓ₀ ∧ ∀ ℓ ≥ ℓ₀, ∀ (T : ℝ) (P : ExpSum),
      T ≤ ℓ ^ 7 → ZeroConfig P T →
      ‖P.eval 0‖ ^ 2 ≤ (zeroSum + ε) / ℓ * (oubil ℓ P.eval (P.dmod 0).eval P.eval (P.dmod 0).eval).re := by
  intro ε hε
  classical
  obtain ⟨κh, hκh, HR⟩ := highR_parts' hT
  obtain ⟨κe, ce, hκe, hce, LE⟩ := lemmaE' hT
  obtain ⟨cf, hcf, CF⟩ := claimF_core' hT
  obtain ⟨κ, hκdef⟩ : ∃ κ : ℝ, κ = max κh κe := ⟨_, rfl⟩
  have hκhκ : κh ≤ κ := by rw [hκdef]; exact le_max_left _ _
  have hκeκ : κe ≤ κ := by rw [hκdef]; exact le_max_right _ _
  have hκ1 : 1 ≤ κ := hκh.trans hκhκ
  have hlk : 0 ≤ Real.log (2 * κ) := Real.log_nonneg (by linarith)
  obtain ⟨ε₁, hε₁def⟩ : ∃ e : ℝ, e = 1 / (180 * (Real.log (2 * κ) + 1)) := ⟨_, rfl⟩
  have hε₁ : 0 < ε₁ := by rw [hε₁def]; positivity
  have hε₁s : 180 * ε₁ * Real.log (2 * κ) ≤ 1 := by
    have e : 180 * ε₁ * Real.log (2 * κ) = Real.log (2 * κ) / (Real.log (2 * κ) + 1) := by
      rw [hε₁def]; field_simp
    rw [e, div_le_one (by linarith)]; linarith
  have hε₁30 : 30 * ε₁ ≤ 1 := by
    have e : 30 * ε₁ = 1 / (6 * (Real.log (2 * κ) + 1)) := by rw [hε₁def]; field_simp; ring
    rw [e, div_le_one (by positivity)]; nlinarith
  have hε₁1 : ε₁ ≤ 1 := by linarith
  have hsm : ∀ κ' : ℝ, 1 ≤ κ' → κ' ≤ κ → 180 * ε₁ * Real.log (2 * κ') ≤ 1 := by
    intro κ' h1 h2
    have : Real.log (2 * κ') ≤ Real.log (2 * κ) := Real.log_le_log (by linarith) (by linarith)
    have : 180 * ε₁ * Real.log (2 * κ') ≤ 180 * ε₁ * Real.log (2 * κ) := by gcongr
    linarith
  have hsm4 : ∀ κ' : ℝ, 1 ≤ κ' → κ' ≤ κ → 4 * ε₁ * Real.log (2 * κ') ≤ 1 / 2 := by
    intro κ' h1 h2
    have h3 := hsm κ' h1 h2
    have h4 : 0 ≤ Real.log (2 * κ') := Real.log_nonneg (by linarith)
    have : 0 ≤ ε₁ * Real.log (2 * κ') := mul_nonneg hε₁.le h4
    nlinarith
  obtain ⟨t₁, ht₁, hcl⟩ := cluster_bound hGG hRH ε₁ hε₁
  obtain ⟨H₁, hH₁⟩ := node_sum_tail F1 F1_nonneg (summable_wF1 hc) (ε / 48) (by positivity)
  obtain ⟨Hs, d, hHs1, hHs2, hd, hd2, hgap⟩ := gap_exists (max H₁ (max (t₁ + 2) 2))
  have hHsH₁ : H₁ ≤ Hs := (le_max_left _ _).trans hHs1
  have hHst : t₁ + 2 ≤ Hs := ((le_max_left _ _).trans (le_max_right _ _)).trans hHs1
  have hHs2' : 2 ≤ Hs := ((le_max_right _ _).trans (le_max_right _ _)).trans hHs1
  obtain ⟨g₀, hg₀, hgap₀⟩ := finset_gap (lowOrd Hs)
  obtain ⟨M₀, hM₀def⟩ : ∃ M : ℕ, M = (lowOrd Hs).sup (fun γ => mult ((1 / 2 : ℂ) + (γ : ℂ) * I)) :=
    ⟨_, rfl⟩
  have hM₀ : ∀ γ ∈ lowOrd Hs, mult ((1 / 2 : ℂ) + (γ : ℂ) * I) ≤ M₀ := fun γ hγ => by
    rw [hM₀def]; exact Finset.le_sup (f := fun γ : ℝ => mult ((1 / 2 : ℂ) + (γ : ℂ) * I)) hγ
  obtain ⟨CF1, hCF1def⟩ : ∃ C : ℝ, C = ∑' ρ : NTZ, (mult (ρ : ℂ) : ℝ) * F1 |(ρ : ℂ).im| := ⟨_, rfl⟩
  have hCF1 : 0 ≤ CF1 := by
    rw [hCF1def]; exact tsum_nonneg fun ρ => mul_nonneg (Nat.cast_nonneg _) (F1_nonneg _)
  have hCF1s : ∀ (P : ExpSum) (T : ℝ), ZeroConfig P T → ∀ T' : Set ℝ,
      ∑ γ ∈ nodeSet P T', (nw P γ : ℝ) * F1 |γ| ≤ CF1 := fun P T hP T' => by
    rw [hCF1def]; exact node_sum_le_tsum hP T' F1 F1_nonneg (summable_wF1 hc) _ (Finset.Subset.refl _)
  obtain ⟨Kall, hKall⟩ : ∃ K : ℝ, K = 3 * (claimFC cf M₀ Hs g₀ (lowOrd Hs).card +
      12 * Real.sqrt (2 * (2 * ce) ^ (2 * (M₀ - 1)) * ((lowOrd Hs).card * ((Hs + d) / (2 * d)) ^ 2 * CF1)) *
        (Hs + markovC ce M₀ + 1) * (2 + markovC ce M₀ + 4 * CB * (2 * κe) ^ 6) +
      (2 * CF1 + 1) * (2 * Real.pi) * KT * (2 * κh) ^ 6 + 1) := ⟨_, rfl⟩
  have hCBp := CB_nonneg
  have hKT := KT_nonneg
  have hKall0 : 0 < Kall := by
    rw [hKall]
    have h1 := claimFC_nonneg (H₀ := Hs) hcf M₀ (by linarith) hg₀ (lowOrd Hs).card
    have h2 := markovC_nonneg (by linarith : (0 : ℝ) < ce) M₀
    have h3 : 0 ≤ 12 * Real.sqrt (2 * (2 * ce) ^ (2 * (M₀ - 1)) * ((lowOrd Hs).card * ((Hs + d) / (2 * d)) ^ 2 * CF1)) *
        (Hs + markovC ce M₀ + 1) * (2 + markovC ce M₀ + 4 * CB * (2 * κe) ^ 6) := by
      have : 0 ≤ Hs + markovC ce M₀ + 1 := by linarith
      positivity
    have h4 : 0 ≤ (2 * CF1 + 1) * (2 * Real.pi) * KT * (2 * κh) ^ 6 := by positivity
    linarith
  set x₀ := ε / (2 * (zeroSum + ε)) with hx₀
  have hx₀p : 0 < x₀ := by have := zeroSum_nonneg; positivity
  obtain ⟨δ₀, hδ₀def⟩ : ∃ δ : ℝ, δ = min (min (1 / 16) (ε₁ / 4))
      (min (1 / (128 * CB ^ 2 * (2 * κ) ^ 6 + 1)) (x₀ / Kall)) := ⟨_, rfl⟩
  have hδ₀ : 0 < δ₀ := by
    rw [hδ₀def]
    exact lt_min (lt_min (by norm_num) (by positivity)) (lt_min (by positivity) (by positivity))
  have hδ1 : δ₀ ≤ 1 / 16 := by rw [hδ₀def]; exact (min_le_left _ _).trans (min_le_left _ _)
  have hδε : δ₀ ≤ ε₁ / 4 := by rw [hδ₀def]; exact (min_le_left _ _).trans (min_le_right _ _)
  have hδB : δ₀ ≤ 1 / (128 * CB ^ 2 * (2 * κ) ^ 6 + 1) := by
    rw [hδ₀def]; exact (min_le_right _ _).trans (min_le_left _ _)
  have hδx : δ₀ * Kall ≤ x₀ := by
    have : δ₀ ≤ x₀ / Kall := by rw [hδ₀def]; exact (min_le_right _ _).trans (min_le_right _ _)
    rwa [le_div_iff₀ hKall0] at this
  rw [hKall] at hδx
  have hev : ∀ᶠ q in atTop, (Real.log q + (8 * A₀ + 10)) ^ 90 / q ∈ Iic δ₀ ∧ 2 ≤ q :=
    ((tendsto_logpow_div 90 (8 * A₀ + 10)).eventually (Iic_mem_nhds hδ₀)).and (eventually_ge_atTop 2)
  obtain ⟨q₁, hq₁⟩ := eventually_atTop.mp hev
  refine ⟨max (q₁ ^ 4) 16, le_max_of_le_right (by norm_num), fun ℓ hℓ T P hTℓ hP => ?_⟩
  have hℓ16 : 16 ≤ ℓ := le_trans (le_max_right _ _) hℓ
  have hℓ0 : 0 < ℓ := by linarith
  obtain ⟨q, hqdef⟩ : ∃ q : ℝ, q = ℓ ^ (1 / 4 : ℝ) := ⟨_, rfl⟩
  have hq0 : 0 ≤ q := by rw [hqdef]; exact Real.rpow_nonneg hℓ0.le _
  have hq4 : q ^ 4 = ℓ := by
    rw [hqdef, ← Real.rpow_natCast, ← Real.rpow_mul hℓ0.le]; norm_num
  have hqq₁ : q₁ ≤ q := by
    by_contra h
    push Not at h
    have : q ^ 4 < q₁ ^ 4 := pow_lt_pow_left₀ h hq0 (by norm_num)
    have : q₁ ^ 4 ≤ ℓ := le_trans (le_max_left _ _) hℓ
    linarith
  have hq2 : 2 ≤ q := by
    by_contra h
    push Not at h
    have : q ^ 4 < 2 ^ 4 := pow_lt_pow_left₀ h hq0 (by norm_num)
    linarith
  obtain ⟨hYq, -⟩ := hq₁ q hqq₁
  rw [← hq4] at hTℓ ⊢
  have hTq : T ≤ q ^ 28 := by rw [show q ^ 28 = (q ^ 4) ^ 7 by ring]; exact hTℓ
  exact u2_config hc hRH HR LE CF hκh hκe hcf hce hκhκ hκeκ hε₁ hε₁30 hε₁1 hsm hsm4 ht₁ hcl hH₁ hHsH₁ hHst
    hHs2' hd hd2 hgap hg₀ hgap₀ hM₀ le_rfl hCF1 hCF1s hε hδ₀ hδ1 hδε hδB hδx hq2 hYq hTq hP

end Riemann.Deep.LeanNBRate

