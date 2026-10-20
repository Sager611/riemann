import Riemann.Fayed.Asymptotics

/-!
# Fayed, arXiv:2209.01890v55 — the whole argument in Lean

s = σ + it, 0 < σ < 1; p ∈ (0,1); η ∈ (0, min(σ,1-σ)); γ ∈ (0,1/4); T = u⁸.

* Literature facts the paper cites (all true, standard) → fields of Inputs (Riemann.Fayed.Defs).
* Paper-internal steps → lemmas eqNN, numbered as in the paper. All proved from Inputs
  (asymptotic machinery in Riemann.Fayed.Asymptotics).
* fayed_RH : Inputs → Theorem1C → RiemannHypothesis — the paper's logic, checked.
* theorem1C_false : ¬ Theorem1C — the argument breaks exactly at (55).
-/

open Complex Filter Topology Set MeasureTheory intervalIntegral
open scoped Real ComplexConjugate

namespace Fayed

noncomputable section

/-! ## Asymptotic steps (proved from Inputs) -/

/-- (39): F/G = u^p F(u) → 0 along D_γ. From (34),(35) and the identity
u^{s+p}J(u,s) - u^{1-s̄+p}J(u,1-s̄) = u^p F(u): u^p F = (1-2σ₀)B̄₂ u^{p-1}/2 + O(u^{p-2}). -/
lemma eq39 (hin : Inputs) (X : Setting) (M₀ : ℕ) :
    Tendsto (fun u => Ffun X.η X.σ₀ X.t₀ u / Gfun X.p u) (along (Dγ X.γ M₀)) (𝓝 0) := by
  obtain ⟨K, U, hU, hK, hFG⟩ := FG_bound hin X X.hγ.1 X.hγ.2
  have hlim : Tendsto (fun u : ℝ => K * u ^ (-(1 - X.p))) atTop (𝓝 0) := by
    simpa using (tendsto_rpow_neg_atTop (by linarith [X.hp.2] : 0 < 1 - X.p)).const_mul K
  refine squeeze_zero_norm' ?_ (hlim.mono_left (along_le_atTop _))
  filter_upwards [(eventually_ge_atTop U).filter_mono (along_le_atTop _),
    eventually_mem_along _] with u hUu huD
  exact hFG M₀ u huD hUu

/-- (42): F → 0 along D_γ. -/
lemma eq42 (hin : Inputs) (X : Setting) (M₀ : ℕ) :
    Tendsto (Ffun X.η X.σ₀ X.t₀) (along (Dγ X.γ M₀)) (𝓝 0) := by
  have h := (eq39 hin X M₀).mul ((eq43 X.hp.1).mono_left (along_le_atTop _))
  rw [zero_mul] at h
  refine h.congr' ?_
  filter_upwards [(eventually_gt_atTop 0).filter_mono (along_le_atTop _)] with u hu
  have : Gfun X.p u ≠ 0 := by
    rw [Gfun]; exact_mod_cast (Real.rpow_pos_of_pos hu _).ne'
  exact div_mul_cancel₀ _ this

/-- (50): σ₀ ≠ 1/2 ⇒ |F'/G'| → ∞ along D_γ; indeed F'/G' = (1-2σ₀)u^p B̄₁(u)/(-p) + o(1)
and B̄₁ ≥ γ on D_γ. From (44)-(49). -/
lemma eq50 (hin : Inputs) (X : Setting) (M₀ : ℕ) (hσ : X.σ₀ ≠ 1 / 2) :
    Tendsto (fun u => deriv (Ffun X.η X.σ₀ X.t₀) u / deriv (Gfun X.p) u) (along (Dγ X.γ M₀))
      (Bornology.cobounded ℂ) := by
  obtain ⟨K1, U1, hU1, hK1, hX⟩ := X_bound hin X X.hγ.1 X.hγ.2
  obtain ⟨K2, U2, hU2, hK2, hB⟩ := bdry_bound hin X.t₀ X.hη.1 X.hη.2.1 X.hη.2.2
  have hp := X.hp.1
  have hm : 0 < |1 - 2 * X.σ₀| * X.γ :=
    mul_pos (abs_pos.2 (fun h => hσ (by linarith))) X.hγ.1
  have hc1 : 1 / (2 * π) ≤ 1 := by
    rw [div_le_one (by positivity)]; linarith [Real.pi_gt_three]
  have hev : ∀ᶠ u in along (Dγ X.γ M₀), |1 - 2 * X.σ₀| * X.γ / (2 * X.p) * u ^ X.p ≤
      ‖deriv (Ffun X.η X.σ₀ X.t₀) u / deriv (Gfun X.p) u‖ := by
    filter_upwards [(eventually_ge_atTop (max (max U1 U2)
      (2 * (K1 + K2) / (|1 - 2 * X.σ₀| * X.γ)))).filter_mono (along_le_atTop _),
      eventually_mem_along _] with u hUu huD
    have hUu1 : U1 ≤ u := le_trans (le_max_left _ _) (le_trans (le_max_left _ _) hUu)
    have hUu2 : U2 ≤ u := le_trans (le_max_right _ _) (le_trans (le_max_left _ _) hUu)
    have hUm : 2 * (K1 + K2) / (|1 - 2 * X.σ₀| * X.γ) ≤ u := le_trans (le_max_right _ _) hUu
    have hu0 : 0 < u := by linarith
    have hw : 0 < u ^ (-X.p - 1) := Real.rpow_pos_of_pos hu0 _
    have hGn : ‖((-X.p * u ^ (-X.p - 1) : ℝ) : ℂ)‖ = X.p * u ^ (-X.p - 1) := by
      rw [Complex.norm_real, Real.norm_eq_abs, abs_of_neg (mul_neg_of_neg_of_pos (by linarith) hw)]
      ring
    have hT : |1 - 2 * X.σ₀| * X.γ ≤ ‖((1 - 2 * X.σ₀ : ℝ) : ℂ) * (Bbar 1 u : ℂ)‖ := by
      have hb := Bbar_one_gt X.hγ.1 huD
      rw [norm_mul, Complex.norm_real, Complex.norm_real, Real.norm_eq_abs, Real.norm_eq_abs,
        abs_of_pos (by linarith [X.hγ.1] : 0 < Bbar 1 u)]
      exact mul_le_mul_of_nonneg_left hb.le (abs_nonneg _)
    have hlow := lower_aux (c := 1 / (2 * π)) hu0 (by positivity) hc1 hT (hX M₀ u huD hUu1)
      (hB u hUu2) hUm hm
    rw [(hin.leibniz44 _ _ _ X.hη.1 X.hη.2.1 X.hη.2.2 u hu0).deriv,
      (hasDerivAt_Gfun X.p hu0).deriv, Fder_eq hu0 X.hη.1 X.hη.2.1 X.hη.2.2, norm_div, hGn]
    calc |1 - 2 * X.σ₀| * X.γ / (2 * X.p) * u ^ X.p
        = |1 - 2 * X.σ₀| * X.γ / (2 * X.p) * (u⁻¹ / u ^ (-X.p - 1)) := by
          rw [inv_div_rpow hu0, show (-1 - (-X.p - 1)) = X.p by ring]
      _ = u⁻¹ * (|1 - 2 * X.σ₀| * X.γ / 2) / (X.p * u ^ (-X.p - 1)) := by ring
      _ ≤ _ := div_le_div_of_nonneg_right hlow (mul_pos hp hw).le
  rw [← tendsto_norm_atTop_iff_cobounded]
  refine tendsto_atTop_mono' _ hev ?_
  exact ((tendsto_rpow_atTop hp).const_mul_atTop (div_pos hm (by linarith))).mono_left
    (along_le_atTop _)

/-- (51)–(54): oscillation of F/G on each D_M is ≤ δ_M → 0. -/
lemma eq54 (hin : Inputs) (X : Setting) (M₀ : ℕ) :
    ∃ δ : ℕ → ℝ, Tendsto δ atTop (𝓝 0) ∧ ∀ M, ∀ u ∈ Ioo (aγ X.γ M₀ M) (bγ X.γ M₀ M),
      ‖Ffun X.η X.σ₀ X.t₀ u / Gfun X.p u -
        Ffun X.η X.σ₀ X.t₀ (bγ X.γ M₀ M) / Gfun X.p (bγ X.γ M₀ M)‖ ≤ δ M := by
  have hγ0 := X.hγ.1
  have hγ1 := X.hγ.2
  obtain ⟨K, U, hU, hK, hFG⟩ := FG_bound hin X (γ := X.γ / 2) (by linarith) (by linarith)
  have hbd : ∀ M, ∃ B, ∀ u ∈ Icc (aγ X.γ M₀ M) (bγ X.γ M₀ M),
      ‖Ffun X.η X.σ₀ X.t₀ u / Gfun X.p u‖ ≤ B := by
    intro M
    apply isCompact_Icc.exists_bound_of_continuousOn
    intro u hu
    have hu0 : 0 < u := (aγ_pos hγ0 M₀ M).trans_le hu.1
    have hG : Gfun X.p u ≠ 0 := by
      rw [Gfun]; exact_mod_cast (Real.rpow_pos_of_pos hu0 _).ne'
    exact ((hin.leibniz44 _ _ _ X.hη.1 X.hη.2.1 X.hη.2.2 u hu0).continuousAt.div
      (hasDerivAt_Gfun X.p hu0).continuousAt hG).continuousWithinAt
  choose B hB using hbd
  have hp1 : 0 < 1 - X.p := by linarith [X.hp.2]
  classical
  refine ⟨fun M => if U ≤ aγ X.γ M₀ M then 2 * (K * aγ X.γ M₀ M ^ (-(1 - X.p))) else 2 * B M,
    ?_, ?_⟩
  · have h1 : Tendsto (fun M => 2 * (K * aγ X.γ M₀ M ^ (-(1 - X.p)))) atTop (𝓝 0) := by
      have := (((tendsto_rpow_neg_atTop hp1).comp (tendsto_aγ X.γ M₀)).const_mul K).const_mul 2
      simpa using this
    refine h1.congr' ?_
    filter_upwards [(tendsto_aγ X.γ M₀).eventually_ge_atTop U] with M hM
    simp [hM]
  · intro M u hu
    have hab : aγ X.γ M₀ M ≤ bγ X.γ M₀ M := (hu.1.trans hu.2).le
    have hnorm := norm_sub_le (Ffun X.η X.σ₀ X.t₀ u / Gfun X.p u)
      (Ffun X.η X.σ₀ X.t₀ (bγ X.γ M₀ M) / Gfun X.p (bγ X.γ M₀ M))
    beta_reduce
    split_ifs with hM
    · have ha0 := aγ_pos hγ0 M₀ M
      have huD : u ∈ Dγ (X.γ / 2) M₀ := mem_iUnion.2 ⟨M, by
        have h1 := hu.1
        have h2 := hu.2
        simp only [aγ, bγ, mem_Ioo] at h1 h2 ⊢
        constructor <;> linarith⟩
      have hbD : bγ X.γ M₀ M ∈ Dγ (X.γ / 2) M₀ := mem_iUnion.2 ⟨M, by
        simp only [aγ, bγ, mem_Ioo]
        constructor <;> linarith⟩
      have e1 := hFG M₀ u huD (hM.trans hu.1.le)
      have e2 := hFG M₀ _ hbD (hM.trans hab)
      have m1 : u ^ (-(1 - X.p)) ≤ aγ X.γ M₀ M ^ (-(1 - X.p)) :=
        Real.rpow_le_rpow_of_nonpos ha0 hu.1.le (by linarith)
      have m2 : bγ X.γ M₀ M ^ (-(1 - X.p)) ≤ aγ X.γ M₀ M ^ (-(1 - X.p)) :=
        Real.rpow_le_rpow_of_nonpos ha0 hab (by linarith)
      have k1 := mul_le_mul_of_nonneg_left m1 hK
      have k2 := mul_le_mul_of_nonneg_left m2 hK
      linarith
    · have e1 := hB M u (Ioo_subset_Icc_self hu)
      have e2 := hB M _ (right_mem_Icc.2 hab)
      linarith

/-! ## The paper's conclusion, (55)–(56) -/

/-- The paper's argument, checked: given the literature inputs AND its Theorem 1, RH follows. -/
theorem fayed_RH (hin : Inputs) (thm1 : Theorem1C) : RiemannHypothesis := by
  intro s hs hnt hs1
  obtain ⟨hσ0, hσ1⟩ := hin.strip s hs hnt hs1
  by_contra hne
  have hm : 0 < min s.re (1 - s.re) := lt_min hσ0 (by linarith)
  have hmle1 : min s.re (1 - s.re) ≤ s.re := min_le_left _ _
  have hmle2 : min s.re (1 - s.re) ≤ 1 - s.re := min_le_right _ _
  let X : Setting :=
    { σ₀ := s.re, t₀ := s.im, p := 1 / 2, η := min s.re (1 - s.re) / 2, γ := 1 / 8,
      hσ₀ := ⟨hσ0, hσ1⟩, hzero := by rw [Complex.re_add_im]; exact hs,
      hp := by norm_num, hη := ⟨by linarith, by linarith, by linarith⟩, hγ := by norm_num }
  obtain ⟨-, -, ha0, hab, hba, hb, -⟩ := counterexample (γ := 1 / 8) (by norm_num) (by norm_num) 0
  obtain ⟨δ, hδ, hosc⟩ := eq54 hin X 0
  have hneBot : (along (Dγ X.γ 0)).NeBot := along_Dγ_neBot (by norm_num) 0
  -- (55): Theorem 1 turns (42),(43),(50),(54) into F/G → ∞
  have h55 := thm1 (Ffun X.η X.σ₀ X.t₀) (Gfun X.p) (aγ X.γ 0) (bγ X.γ 0) δ
    (Bornology.cobounded ℂ) (Or.inl rfl)
    (fun u hu => (hin.leibniz44 _ _ _ X.hη.1 X.hη.2.1 X.hη.2.2 u hu).differentiableAt)
    (fun u hu => (hasDerivAt_Gfun X.p hu).differentiableAt)
    ha0 hab hba hb (eq42 hin X 0) ((eq43 X.hp.1).mono_left (along_le_atTop _))
    (eq50 hin X 0 hne) hosc hδ
  -- (56): contradiction with (39)
  rw [← tendsto_norm_atTop_iff_cobounded] at h55
  have h39 := (eq39 hin X 0).norm
  rw [norm_zero] at h39
  exact not_tendsto_nhds_of_tendsto_atTop h55 0 h39

/-- **Where the paper breaks.** Its argument is valid given its Theorem 1 (fayed_RH),
and Theorem 1 is false (theorem1C_false). The failing step is (55). -/
theorem fayed_gap : (Inputs → Theorem1C → RiemannHypothesis) ∧ ¬ Theorem1C :=
  ⟨fayed_RH, theorem1C_false⟩

end

end Fayed
