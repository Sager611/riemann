import Riemann.Deep.LeanNBRate.U2props

/-!
# U2 main theorem
-/

open MeasureTheory Filter Topology Set Complex ComplexConjugate

noncomputable section

namespace Riemann.Deep.LeanNBRate

set_option maxHeartbeats 2000000 in
/-- U2 for one configuration at ℓ = q⁴, all constants fixed. -/
theorem u2_config (hc : Titchmarsh_Thm9_2) (hRH : RiemannHypothesis)
    {κh κe ce cf κ ε₁ t₁ H₁ Hs d g₀ CF1 ε δ₀ q T : ℝ} {M₀ J₀ : ℕ} {P : ExpSum}
    (HR : HRprop κh) (LE : LEprop κe ce) (CF : CFprop cf)
    (hκh : 1 ≤ κh) (hκe : 1 ≤ κe) (hcf : 1 ≤ cf) (hce : 1 ≤ ce) (hκhκ : κh ≤ κ) (hκeκ : κe ≤ κ)
    (hε₁ : 0 < ε₁) (hε₁30 : 30 * ε₁ ≤ 1) (hε₁1 : ε₁ ≤ 1)
    (hsm : ∀ κ' : ℝ, 1 ≤ κ' → κ' ≤ κ → 180 * ε₁ * Real.log (2 * κ') ≤ 1)
    (hsm4 : ∀ κ' : ℝ, 1 ≤ κ' → κ' ≤ κ → 4 * ε₁ * Real.log (2 * κ') ≤ 1 / 2)
    (ht₁ : 3 ≤ t₁) (hcl : ∀ t ≥ t₁, ∀ h, 0 < h → h ≤ ε₁ → (Ncount t (t + h) : ℝ) ≤ ε₁ * Real.log t)
    (hH₁ : ∀ H ≥ H₁, ∀ (P : ExpSum) (T : ℝ), ZeroConfig P T → ∀ (T' : Set ℝ) (A : Finset ℝ),
      A ⊆ nodeSet P T' → (∀ γ ∈ A, H < |γ|) → ∑ γ ∈ A, (nw P γ : ℝ) * F1 |γ| ≤ ε / 48)
    (hHsH₁ : H₁ ≤ Hs) (hHst : t₁ + 2 ≤ Hs) (hHs2 : 2 ≤ Hs) (hd : 0 < d) (hd2 : d ≤ 1 / 2)
    (hgap : ∀ ρ ∈ NTZ, |ρ.im| ≤ Hs - d ∨ Hs + d ≤ |ρ.im|)
    (hg₀ : 0 < g₀) (hgap₀ : ∀ γ ∈ lowOrd Hs, ∀ γ' ∈ lowOrd Hs, γ ≠ γ' → g₀ ≤ |γ - γ'|)
    (hM₀ : ∀ γ ∈ lowOrd Hs, mult ((1 / 2 : ℂ) + (γ : ℂ) * I) ≤ M₀) (hJ₀ : (lowOrd Hs).card ≤ J₀)
    (hCF1 : 0 ≤ CF1) (hCF1s : ∀ (P : ExpSum) (T : ℝ), ZeroConfig P T → ∀ T' : Set ℝ,
      ∑ γ ∈ nodeSet P T', (nw P γ : ℝ) * F1 |γ| ≤ CF1)
    (hε : 0 < ε) (hδ₀ : 0 < δ₀) (hδ1 : δ₀ ≤ 1 / 16) (hδε : δ₀ ≤ ε₁ / 4)
    (hδB : δ₀ ≤ 1 / (128 * CB ^ 2 * (2 * κ) ^ 6 + 1))
    (hδx : δ₀ * (3 * (claimFC cf M₀ Hs g₀ J₀ +
        12 * Real.sqrt (2 * (2 * ce) ^ (2 * (M₀ - 1)) * (J₀ * ((Hs + d) / (2 * d)) ^ 2 * CF1)) *
          (Hs + markovC ce M₀ + 1) * (2 + markovC ce M₀ + 4 * CB * (2 * κe) ^ 6) +
        (2 * CF1 + 1) * (2 * Real.pi) * KT * (2 * κh) ^ 6 + 1)) ≤ ε / (2 * (zeroSum + ε)))
    (hq2 : 2 ≤ q) (hYq : (Real.log q + (8 * A₀ + 10)) ^ 90 / q ≤ δ₀)
    (hTq : T ≤ q ^ 28) (hP : ZeroConfig P T) :
    ‖P.eval 0‖ ^ 2 ≤ (zeroSum + ε) / q ^ 4 *
      (oubil (q ^ 4) P.eval (P.dmod 0).eval P.eval (P.dmod 0).eval).re := by
  classical
  have hq1 : 1 ≤ q := by linarith
  have hqp : 0 < q := by linarith
  have hA0 := A₀_pos
  have hCBp := CB_nonneg
  -- height
  obtain ⟨T', hT'def⟩ : ∃ T', T' = max T 1 := ⟨_, rfl⟩
  have hP' : ZeroConfig P T' := hP.mono (by rw [hT'def]; exact le_max_left _ _)
  have hT'1 : 1 ≤ T' := by rw [hT'def]; exact le_max_right _ _
  have hT'q : T' ≤ q ^ 28 := by rw [hT'def]; exact max_le hTq (one_le_pow₀ hq1)
  -- Y and the multiplicity cutoff
  obtain ⟨Y, hYdef⟩ : ∃ Y, Y = Real.log q + (8 * A₀ + 10) := ⟨_, rfl⟩
  rw [← hYdef] at hYq
  have hlq : 0 ≤ Real.log q := Real.log_nonneg hq1
  have hY1 : 1 ≤ Y := by rw [hYdef]; linarith
  have hYk : ∀ k ≤ 90, Y ^ k ≤ δ₀ * q := fun k hk => Ypow_le hq1 hY1 hYq hk
  obtain ⟨Ku, hKudef⟩ : ∃ Ku, Ku = KuOf ε₁ q := ⟨_, rfl⟩
  have hKu1 : (Ku : ℝ) ≤ Real.log q + 1 := by rw [hKudef]; exact KuOf_le hε₁.le hε₁30 hq1
  have hKuY : (Ku : ℝ) + 1 + 8 * A₀ ≤ Y := by rw [hYdef]; linarith
  have hKuY' : (Ku : ℝ) + 1 ≤ Y := by linarith
  have hAY : A₀ ≤ Y := by rw [hYdef]; linarith
  have hA16 : 1 + 16 * A₀ ≤ 2 * Y := by rw [hYdef]; linarith
  have hKuge : ε₁ * Real.log (T' + 2) ≤ Ku := by rw [hKudef]; exact KuOf_ge hε₁.le hq2 (by linarith) hT'q
  have hpowK : ∀ κ' : ℝ, 1 ≤ κ' → κ' ≤ κ → ∀ M : ℕ, M ≤ Ku → (2 * κ') ^ (6 * M) ≤ (2 * κ') ^ 6 * q := by
    intro κ' h1 h2 M hM
    have h3 : (2 * κ') ^ (6 * M) ≤ (2 * κ') ^ (6 * Ku) := pow_le_pow_right₀ (by linarith) (by omega)
    refine h3.trans ?_
    rw [hKudef]; exact pow_KuOf_le h1 hε₁.le (hsm κ' h1 h2) hq1
  have hKuq : (Ku : ℝ) ≤ q := condD hq1 hY1 hYq (by linarith) hKuY'
  -- scales
  have hs₀ : (0 : ℝ) < 1 / q ^ 3 := by positivity
  have hη : (0 : ℝ) < 1 / q ^ 2 := by positivity
  have hη' : (0 : ℝ) < 1 / q ^ 4 := by positivity
  have hwmax0 : (0 : ℝ) ≤ Ku * (1 / q ^ 3) := by positivity
  have hwmaxη : (Ku : ℝ) * (1 / q ^ 3) ≤ 1 / q ^ 2 := by
    rw [mul_one_div, div_le_div_iff₀ (by positivity) (by positivity)]
    calc (Ku : ℝ) * q ^ 2 ≤ q * q ^ 2 := mul_le_mul_of_nonneg_right hKuq (by positivity)
      _ = 1 * q ^ 3 := by ring
  have hs₀1 : (1 : ℝ) / q ^ 3 ≤ 1 := by rw [div_le_one (by positivity)]; exact one_le_pow₀ hq1
  have hη1 : (1 : ℝ) / q ^ 2 ≤ 1 := by rw [div_le_one (by positivity)]; exact one_le_pow₀ hq1
  have hq4ge : (1 : ℝ) ≤ q ^ 4 := one_le_pow₀ hq1
  have hq4p : (0 : ℝ) < q ^ 4 := by positivity
  have hKu2 : ((Ku : ℝ) + 2) * (1 / q ^ 3) ≤ ε₁ / 2 := condE hq1 hY1 hYq hδε hKuY'
  -- high groups
  have hG := hi_group_facts hP' (H := Hs) hs₀ hε₁ hε₁1 hHst (by linarith) hcl hKu2 hKuge
  have hgapd : ∀ x ∈ nodeSet P {γ | Hs < |γ|}, Hs + d ≤ |x| := by
    intro x hx
    have h1 := node_mem_NTZ hP' hx
    have h2 : Hs < |x| := (mem_nodeSet.mp hx).2
    rcases hgap _ h1 with h | h
    · rw [im_half_add] at h; linarith
    · rwa [im_half_add] at h
  have hgapl : ∀ x ∈ nodeSet P {γ | |γ| ≤ Hs}, |x| ≤ Hs - d := by
    intro x hx
    have h1 := node_mem_NTZ hP' hx
    have h2 : |x| ≤ Hs := (mem_nodeSet.mp hx).2
    rcases hgap _ h1 with h | h
    · rwa [im_half_add] at h
    · rw [im_half_add] at h; linarith
  have hcen : ∀ g ∈ groups P {γ | Hs < |γ|} (1 / q ^ 3), Hs + d ≤ |gC g| := fun g hg =>
    groups_center P _ hs₀ (by linarith) hgapd hg
  have hgW : ∀ g ∈ groups P {γ | Hs < |γ|} (1 / q ^ 3), gW g ≤ 1 / q ^ 2 := by
    intro g hg
    obtain ⟨_, _, hMK, hW, _⟩ := hG g hg
    have : ((P.part ↑g (gC g)).M : ℝ) * (1 / q ^ 3) ≤ Ku * (1 / q ^ 3) := by gcongr
    linarith
  have hgW1 : ∀ g ∈ groups P {γ | Hs < |γ|} (1 / q ^ 3), gW g ≤ 1 := fun g hg => (hgW g hg).trans hη1
  have hsmallB : ∀ κ' : ℝ, 1 ≤ κ' → κ' ≤ κ → ∀ M : ℕ, M ≤ Ku → ∀ w, 0 ≤ w → w ≤ 1 / q ^ 2 →
      smallB κ' M w (1 / q ^ 2) ≤ 1 / 2 := by
    intro κ' h1 h2 M hM w hw0 hw
    have hMY : (M : ℝ) + 1 ≤ Y := by
      have : (M : ℝ) ≤ Ku := by exact_mod_cast hM
      linarith
    exact (smallB_le h1 hq1 hw0 hw (hpowK κ' h1 h2 M hM) hMY).trans
      (smallB_half h1 h2 hq1 (hYk 6 (by norm_num)) hδ₀.le hδB)
  -- low facts
  have hlo_sub : nodeSet P {γ | |γ| ≤ Hs} ⊆ lowOrd Hs := nodeSet_low_subset hP' Hs
  have hlo_gap : ∀ γ ∈ nodeSet P {γ | |γ| ≤ Hs}, ∀ γ' ∈ nodeSet P {γ | |γ| ≤ Hs}, γ ≠ γ' →
      g₀ ≤ |γ - γ'| := fun γ hγ γ' hγ' hne => hgap₀ γ (hlo_sub hγ) γ' (hlo_sub hγ') hne
  have hlo_H : ∀ γ ∈ nodeSet P {γ | |γ| ≤ Hs}, |γ| ≤ Hs := fun γ hγ => (mem_nodeSet.mp hγ).2
  have hlo_M : ∀ γ ∈ nodeSet P {γ | |γ| ≤ Hs}, 1 ≤ (P.part {γ} γ).M ∧ (P.part {γ} γ).M ≤ M₀ := by
    intro γ hγ
    rw [part_singleton_M]
    exact ⟨nw_pos hP' hγ, (nw_le_mult hP' hγ).trans (hM₀ γ (hlo_sub hγ))⟩
  have hlo_card : (nodeSet P {γ | |γ| ≤ Hs}).card ≤ J₀ := (Finset.card_le_card hlo_sub).trans hJ₀
  -- the three estimates
  obtain ⟨hCFa, hCFb⟩ := CF P {γ | |γ| ≤ Hs} (q ^ 4) g₀ Hs M₀ hq4ge hg₀ (by linarith) hlo_gap hlo_H hlo_M
  have hgrp1 : ∀ g ∈ groups P {γ | Hs < |γ|} (1 / q ^ 3),
      (P.part ↑g (gC g)).M ≤ Ku ∧ gW g ≤ Ku * (1 / q ^ 3) := by
    intro g hg
    obtain ⟨_, _, hMK, hW, _⟩ := hG g hg
    refine ⟨hMK, hW.trans ?_⟩
    gcongr
  have hcA : 4 * A₀ / (1 / q ^ 2) * ((Ku : ℝ) + 1) ^ 2 ≤ q ^ 4 := condA hq1 hY1 hYq hδ1 hKuY' hAY
  have hcB : (Ku : ℝ) * (4 * A₀ / (1 / q ^ 3 / 4)) * ((Ku : ℝ) + 1) ^ 2 ≤ q ^ 4 :=
    condB hq1 hY1 hYq hδ1 hKuY' hAY
  have hsBh : smallB κh Ku (Ku * (1 / q ^ 3)) (1 / q ^ 2) ≤ Hs ^ 2 / 2 := by
    have := hsmallB κh hκh hκhκ Ku le_rfl _ hwmax0 hwmaxη
    nlinarith
  obtain ⟨hHRa, hHRb⟩ := HR P {γ | Hs < |γ|} (1 / q ^ 3) (q ^ 4) Hs (1 / q ^ 2) (1 / q ^ 4)
    (Ku * (1 / q ^ 3)) 40 Ku hs₀ hq4p (by linarith) hη hη' (by norm_num)
    (fun j hj => le_of_lt hj) (fun j _ => hP.mult_pos j) hgrp1 hcA hcB hsBh
  have hLEg : ∀ g ∈ groups P {γ | Hs < |γ|} (1 / q ^ 3), 1 ≤ (P.part ↑g (gC g)).M ∧ 1 ≤ |gC g| ∧
      0 ≤ gW g ∧ (∀ j, (P.part ↑g (gC g)).node j ∈ Icc (-gW g) (gW g)) ∧
      4 * A₀ / (1 / q ^ 2) * (((P.part ↑g (gC g)).M : ℝ) + 1) ^ 2 ≤ q ^ 4 ∧
      smallB κe (P.part ↑g (gC g)).M (gW g) (1 / q ^ 2) ≤ (gC g) ^ 2 / 2 := by
    intro g hg
    obtain ⟨hM1, _, hMK, _, _⟩ := hG g hg
    have hc1 : 1 ≤ |gC g| := by have := hcen g hg; linarith
    have hW0 : 0 ≤ gW g := gW_nonneg (groups_nonempty hg)
    have hMY : ((P.part ↑g (gC g)).M : ℝ) + 1 ≤ Y := by
      have : ((P.part ↑g (gC g)).M : ℝ) ≤ Ku := by exact_mod_cast hMK
      linarith
    refine ⟨hM1, hc1, hW0, P.part_node_mem (↑g) hW0 fun j hj => abs_sub_gC_le (by exact_mod_cast hj),
      condA hq1 hY1 hYq hδ1 hMY hAY, ?_⟩
    have := hsmallB κe hκe hκeκ _ hMK (gW g) hW0 (hgW g hg)
    have hc2 : 1 ≤ (gC g) ^ 2 := by rw [← sq_abs]; nlinarith
    linarith
  have hne : ∀ γ ∈ nodeSet P {γ | |γ| ≤ Hs}, ∀ g ∈ groups P {γ | Hs < |γ|} (1 / q ^ 3), γ ≠ gC g := by
    intro γ hγ g hg h
    have h1 := hlo_H γ hγ
    have h2 := hcen g hg
    rw [h] at h1
    linarith
  have hLE := LE P {γ | |γ| ≤ Hs} {γ | Hs < |γ|} (1 / q ^ 3) (q ^ 4) (1 / q ^ 2) Hs M₀ hq4ge hη
    (by linarith) hlo_H hlo_M hLEg hne
  -- the split p = p_lo + p_hi
  have hdisj : Disjoint {γ : ℝ | |γ| ≤ Hs} {γ | Hs < |γ|} := Set.disjoint_left.mpr fun x h1 h2 => by
    simp only [Set.mem_ofPred_eq] at h1 h2; linarith
  have hcov : ∀ j, P.node j ∈ {γ : ℝ | |γ| ≤ Hs} ∪ {γ | Hs < |γ|} := fun j => by
    by_cases h : |P.node j| ≤ Hs
    · exact Or.inl h
    · exact Or.inr (lt_of_not_ge h)
  have hpe : P.eval = fun u => (P.part {γ | |γ| ≤ Hs} 0).eval u + (P.part {γ | Hs < |γ|} 0).eval u :=
    funext (part_eval_union P _ _ hdisj hcov)
  have hpe' : (P.dmod 0).eval = fun u => ((P.part {γ | |γ| ≤ Hs} 0).dmod 0).eval u +
      ((P.part {γ | Hs < |γ|} 0).dmod 0).eval u := funext (part_dmod_eval_union P _ _ hdisj hcov)
  have hQ := oubil_add_re (q ^ 4) _ _ _ _ (P.part {γ | |γ| ≤ Hs} 0).continuous_eval
    ((P.part {γ | |γ| ≤ Hs} 0).dmod 0).continuous_eval (P.part {γ | Hs < |γ|} 0).continuous_eval
    ((P.part {γ | Hs < |γ|} 0).dmod 0).continuous_eval
  rw [← hpe, ← hpe'] at hQ
  have hQ0 : 0 ≤ (oubil (q ^ 4) P.eval (P.dmod 0).eval P.eval (P.dmod 0).eval).re := by
    rw [oubil_self_re _ hq4p.le _ _ P.continuous_eval (P.dmod 0).continuous_eval]
    exact intervalIntegral.integral_nonneg hq4p.le fun u _ => by positivity
  have hpsi := psi_mass_le_oubil hq4p _ _ (P.part {γ | Hs < |γ|} 0).continuous_eval
    ((P.part {γ | Hs < |γ|} 0).dmod 0).continuous_eval
  have hp : ‖P.eval 0‖ ≤ ‖(P.part {γ | |γ| ≤ Hs} 0).eval 0‖ + ‖(P.part {γ | Hs < |γ|} 0).eval 0‖ := by
    rw [hpe]; exact norm_add_le _ _
  -- the sums
  have hThi2 : ∀ γ ∈ nodeSet P {γ | Hs < |γ|}, 2 ≤ |γ| := fun γ hγ => by have := hgapd γ hγ; linarith
  have hCF1' := hCF1s P T' hP' {γ | Hs < |γ|}
  have hSloZ := low_sum_le_zeroSum hc hRH hP' {γ | |γ| ≤ Hs}
  have hShiε : ∑ g ∈ groups P {γ | Hs < |γ|} (1 / q ^ 3),
      16 * (2 * κh) ^ (4 * ((P.part ↑g (gC g)).M - 1)) / (gC g) ^ 2 ≤ ε / 3 := by
    have h1 := groups_weight_le hP' (Th := {γ | Hs < |γ|}) (s₀ := 1 / q ^ 3) hκh hε₁.le
      (hsm4 κh hκh hκhκ) (fun g hg => (hG g hg).2.1) hgW1 hThi2
    have h2 := hH₁ Hs hHsH₁ P T' hP' {γ | Hs < |γ|} (nodeSet P {γ | Hs < |γ|}) (Finset.Subset.refl _)
      (fun γ hγ => (mem_nodeSet.mp hγ).2)
    have e : ∑ g ∈ groups P {γ | Hs < |γ|} (1 / q ^ 3),
        16 * (2 * κh) ^ (4 * ((P.part ↑g (gC g)).M - 1)) / (gC g) ^ 2 =
        16 * ∑ g ∈ groups P {γ | Hs < |γ|} (1 / q ^ 3),
          (2 * κh) ^ (4 * ((P.part ↑g (gC g)).M - 1)) / (gC g) ^ 2 := by
      rw [Finset.mul_sum]; refine Finset.sum_congr rfl fun g _ => ?_; ring
    rw [e]
    linarith
  have hθ := theta_final (σ := 2 * (2 * ce) ^ (2 * (M₀ - 1)))
    (K := J₀ * ((Hs + d) / (2 * d)) ^ 2 * CF1) (H := Hs) (m := markovC ce M₀) (C := CB) (κ := κe)
    (by positivity) (by positivity) (by linarith) (markovC_nonneg (by linarith) M₀) hCBp (by linarith)
    hq1 hY1 (hYk 3 (by norm_num))
    (u2_theta hP' (Tlo := {γ | |γ| ≤ Hs}) (Thi := {γ | Hs < |γ|}) (s₀ := 1 / q ^ 3) (Y := Y) (ε₁ := ε₁)
      (mk := markovC ce M₀) (CF1 := CF1) (Ku := Ku) (J₀ := J₀) hq1 hκe hε₁.le (hsm4 κe hκe hκeκ) hd
      (by linarith) (by linarith) (markovC_nonneg (by linarith) M₀) hgapl hcen hThi2
      (fun g hg => ⟨(hG g hg).2.2.1, (hG g hg).2.1⟩)
      (fun g hg => ⟨gW_nonneg (groups_nonempty hg), hgW g hg, hgW1 g hg⟩)
      (hpowK κe hκe hκeκ) hKuY' hlo_card hCF1')
  -- the Lemma D defect
  have hw' : (Ku : ℝ) * (1 / q ^ 3) + 1 / q ^ 3 / 4 ≤ Y / q ^ 3 := by
    rw [show (Ku : ℝ) * (1 / q ^ 3) + 1 / q ^ 3 / 4 = ((Ku : ℝ) + 1 / 4) / q ^ 3 by ring]
    apply div_le_div_of_nonneg_right _ (by positivity)
    linarith
  have hCδle := (u2_cdelta hκh hq1 hwmax0 hw' hKuY hA16 (hpowK κh hκh hκhκ Ku le_rfl)
    (Nat.cast_nonneg _) (u2_card hP' hT'1 hT'q {γ | Hs < |γ|} (s₀ := 1 / q ^ 3) hCF1')).trans
    (cdelta_final (by have := KT_nonneg; positivity) hq1 (hYk 84 (by norm_num)) hδ₀.le)
  have hCδ0 : 0 ≤ ((1 + (1 / q ^ 4)⁻¹) * ((groups P {γ | Hs < |γ|} (1 / q ^ 3)).card : ℝ) + 1) *
      (2 * Real.pi) * tailC κh Ku (Ku * (1 / q ^ 3)) (q ^ 4) (1 / q ^ 3) 40 := by
    have : 0 ≤ tailC κh Ku (Ku * (1 / q ^ 3)) (q ^ 4) (1 / q ^ 3) 40 := by
      unfold tailC; have := B₀_pos; positivity
    positivity
  have hinvq : 1 / q ≤ δ₀ := by
    have h1 := hYk 0 (by norm_num)
    rw [pow_zero] at h1
    rw [div_le_iff₀ hqp]; linarith
  have hq4q : 1 / q ^ 4 ≤ 1 / q := by
    apply div_le_div_of_nonneg_left (by norm_num) hqp
    calc q = q ^ 1 := (pow_one q).symm
      _ ≤ q ^ 4 := pow_le_pow_right₀ hq1 (by norm_num)
  obtain ⟨hb0, hble⟩ := b_final hCδ0 hη' (hq4q.trans hinvq) hCδle
  -- name the quantities
  generalize hCδ : ((1 + (1 / q ^ 4)⁻¹) * ((groups P {γ | Hs < |γ|} (1 / q ^ 3)).card : ℝ) + 1) *
      (2 * Real.pi) * tailC κh Ku (Ku * (1 / q ^ 3)) (q ^ 4) (1 / q ^ 3) 40 = Cδ at hHRb hb0 hble
  have hhi : (1 - (1 - (1 - Cδ) / (1 + 1 / q ^ 4))) *
      (∑ g ∈ groups P {γ | Hs < |γ|} (1 / q ^ 3), ∫ x : ℝ, ‖grpQ (P.part ↑g (gC g)) (gC g) (q ^ 4) x‖ ^ 2) ≤
      (oubil (q ^ 4) (P.part {γ | Hs < |γ|} 0).eval ((P.part {γ | Hs < |γ|} 0).dmod 0).eval
        (P.part {γ | Hs < |γ|} 0).eval ((P.part {γ | Hs < |γ|} 0).dmod 0).eval).re := by
    rw [sub_sub_cancel, div_mul_eq_mul_div, div_le_iff₀ (by positivity)]
    calc (1 - Cδ) * _ ≤ _ := hHRb
      _ ≤ (1 + 1 / q ^ 4) * (oubil (q ^ 4) (P.part {γ | Hs < |γ|} 0).eval
          ((P.part {γ | Hs < |γ|} 0).dmod 0).eval (P.part {γ | Hs < |γ|} 0).eval
          ((P.part {γ | Hs < |γ|} 0).dmod 0).eval).re := mul_le_mul_of_nonneg_left hpsi (by positivity)
      _ = _ := by ring
  have hBle := (Complex.abs_re_le_norm _).trans hLE
  have hDlo : 0 ≤ ∑ γ ∈ nodeSet P {γ | |γ| ≤ Hs}, (γ ^ 2 + 1 / 4) *
      ∫ u in (0 : ℝ)..q ^ 4, ‖(P.part {γ} γ).eval u‖ ^ 2 := Finset.sum_nonneg fun γ _ =>
    mul_nonneg (by positivity) (intervalIntegral.integral_nonneg hq4p.le fun u _ => by positivity)
  have hYg : 0 ≤ ∑ g ∈ groups P {γ | Hs < |γ|} (1 / q ^ 3),
      ∫ x : ℝ, ‖grpQ (P.part ↑g (gC g)) (gC g) (q ^ 4) x‖ ^ 2 :=
    Finset.sum_nonneg fun g _ => integral_nonneg fun x => by positivity
  have hSlo0 : 0 ≤ ∑ γ ∈ nodeSet P {γ | |γ| ≤ Hs}, ((P.part {γ} γ).M : ℝ) ^ 2 / (γ ^ 2 + 1 / 4) :=
    Finset.sum_nonneg fun γ _ => by positivity
  have hShi0 : 0 ≤ ∑ g ∈ groups P {γ | Hs < |γ|} (1 / q ^ 3),
      16 * (2 * κh) ^ (4 * ((P.part ↑g (gC g)).M - 1)) / (gC g) ^ 2 :=
    Finset.sum_nonneg fun g _ => by positivity
  have ha0 : 0 ≤ claimFC cf M₀ Hs g₀ (nodeSet P {γ | |γ| ≤ Hs}).card / q ^ 4 :=
    div_nonneg (claimFC_nonneg (H₀ := Hs) hcf M₀ (by linarith) hg₀ _) hq4p.le
  have ha_le : claimFC cf M₀ Hs g₀ (nodeSet P {γ | |γ| ≤ Hs}).card / q ^ 4 ≤ claimFC cf M₀ Hs g₀ J₀ * δ₀ := by
    have h1 := claimFC_mono (H₀ := Hs) hcf M₀ (by linarith) hg₀ hlo_card
    calc claimFC cf M₀ Hs g₀ (nodeSet P {γ | |γ| ≤ Hs}).card / q ^ 4 ≤ claimFC cf M₀ Hs g₀ J₀ / q ^ 4 := by
          gcongr
      _ = claimFC cf M₀ Hs g₀ J₀ * (1 / q ^ 4) := by ring
      _ ≤ claimFC cf M₀ Hs g₀ J₀ * δ₀ :=
          mul_le_mul_of_nonneg_left (hq4q.trans hinvq) (claimFC_nonneg (H₀ := Hs) hcf M₀ (by linarith) hg₀ _)
  have hθ0 : 0 ≤ 2 * Real.sqrt (2 * (2 * ce) ^ (2 * (M₀ - 1)) / q ^ 4) *
      Real.sqrt (∑ γ ∈ nodeSet P {γ | |γ| ≤ Hs}, ∑ g ∈ groups P {γ | Hs < |γ|} (1 / q ^ 3),
        phiE κe Hs (markovC ce M₀) (q ^ 4) (P.part ↑g (gC g)).M (gW g) (1 / q ^ 2) ^ 2 / (γ - gC g) ^ 2) := by
    positivity
  generalize hDdef : (∑ γ ∈ nodeSet P {γ | |γ| ≤ Hs}, (γ ^ 2 + 1 / 4) *
      ∫ u in (0 : ℝ)..q ^ 4, ‖(P.part {γ} γ).eval u‖ ^ 2) = Dlo at hCFa hCFb hBle hDlo
  generalize hYdef' : (∑ g ∈ groups P {γ | Hs < |γ|} (1 / q ^ 3),
      ∫ x : ℝ, ‖grpQ (P.part ↑g (gC g)) (gC g) (q ^ 4) x‖ ^ 2) = Yg at hHRa hhi hBle hYg
  generalize hθdef : 2 * Real.sqrt (2 * (2 * ce) ^ (2 * (M₀ - 1)) / q ^ 4) *
      Real.sqrt (∑ γ ∈ nodeSet P {γ | |γ| ≤ Hs}, ∑ g ∈ groups P {γ | Hs < |γ|} (1 / q ^ 3),
        phiE κe Hs (markovC ce M₀) (q ^ 4) (P.part ↑g (gC g)).M (gW g) (1 / q ^ 2) ^ 2 / (γ - gC g) ^ 2) = θ
    at hBle hθ hθ0
  generalize haDef : claimFC cf M₀ Hs g₀ (nodeSet P {γ | |γ| ≤ Hs}).card / q ^ 4 = a at hCFa ha0 ha_le
  generalize hbDef : 1 - (1 - Cδ) / (1 + 1 / q ^ 4) = b at hhi hb0 hble
  have hsum : a + b + θ ≤ ε / (2 * (zeroSum + ε)) / 3 := by
    have e : δ₀ * (3 * (claimFC cf M₀ Hs g₀ J₀ +
        12 * Real.sqrt (2 * (2 * ce) ^ (2 * (M₀ - 1)) * (J₀ * ((Hs + d) / (2 * d)) ^ 2 * CF1)) *
          (Hs + markovC ce M₀ + 1) * (2 + markovC ce M₀ + 4 * CB * (2 * κe) ^ 6) +
        (2 * CF1 + 1) * (2 * Real.pi) * KT * (2 * κh) ^ 6 + 1)) =
        3 * (claimFC cf M₀ Hs g₀ J₀ * δ₀ +
        12 * Real.sqrt (2 * (2 * ce) ^ (2 * (M₀ - 1)) * (J₀ * ((Hs + d) / (2 * d)) ^ 2 * CF1)) *
          (Hs + markovC ce M₀ + 1) * (2 + markovC ce M₀ + 4 * CB * (2 * κe) ^ 6) * δ₀ +
        (2 * CF1 + 1) * (2 * Real.pi) * KT * (2 * κh) ^ 6 * δ₀ + δ₀) := by ring
    rw [e] at hδx
    linarith
  have hx₀1 : ε / (2 * (zeroSum + ε)) ≤ 1 / 2 := by
    have hZ := zeroSum_nonneg
    rw [div_le_iff₀ (by positivity)]; linarith
  have hab : a + b + θ < 1 := by linarith
  have hmain := combine_bound hq4p hDlo hYg hSlo0 hShi0 hab hθ0 ha0 hb0 hQ hCFa hhi hBle hp
    (norm_nonneg _) (norm_nonneg _) (norm_nonneg _) hCFb hHRa
  rw [show 1 - a - b - θ = 1 - (a + b + θ) by ring] at hmain
  exact final_arith hε zeroSum_nonneg hQ0 hq4p (by linarith) (by linarith : a + b + θ ≤ ε / (2 * (zeroSum + ε)))
    rfl hSloZ hShiε hmain

end Riemann.Deep.LeanNBRate

