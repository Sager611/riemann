import Riemann.Deep.LeanNBRate.LemmaD
import Riemann.Deep.LeanNBRate.Config

/-!
# [N] Theorem R″, high part: endpoint value of the high groups by their derivative mass

|p_hi(0)|² ≤ (Σ_g B(M_g)/γ_g²)/ℓ · (1+η')/(1 - C δ) · ∫|p_hi′ ψ|²  (Lemma B per group, Lemma D).
-/

open MeasureTheory Filter Topology Set Complex

noncomputable section

namespace Riemann.Deep.LeanNBRate

/-- The smallness quantity of Lemma B. -/
def smallB (κ : ℝ) (M : ℕ) (w η : ℝ) : ℝ :=
  4 * (2 * κ) ^ (2 * (M - 1)) * ((2 * κ) ^ (M - 1) * CB * ((M : ℝ) + 1) ^ 3 *
    ((w + η) ^ 2 / η)) ^ 2 * ((2 * κ) ^ (M - 1)) ^ 2

theorem CB_nonneg : 0 ≤ CB := by unfold CB; have := B₀_pos; have := A₀_pos; positivity

theorem pow_mono_two {κ₁ κ : ℝ} (h1 : 0 < κ₁) (hκ : κ₁ ≤ κ) (hκ1 : 1 ≤ κ) {m n : ℕ} (hmn : m ≤ n) :
    (2 * κ₁) ^ m ≤ (2 * κ) ^ n :=
  (pow_le_pow_left₀ (by linarith) (by linarith) m).trans (pow_le_pow_right₀ (by linarith) hmn)

theorem smallB_mono {κ₁ κ : ℝ} (h1 : 0 < κ₁) (hκ : κ₁ ≤ κ) (hκ1 : 1 ≤ κ) {M M' : ℕ}
    (hM : M ≤ M') {w w' η : ℝ} (hw : 0 ≤ w) (hww : w ≤ w') (hη : 0 < η) :
    smallB κ₁ M w η ≤ smallB κ M' w' η := by
  unfold smallB
  have hCB := CB_nonneg
  have a1 := pow_mono_two h1 hκ hκ1 (Nat.mul_le_mul_left 2 (Nat.sub_le_sub_right hM 1))
  have a2 := pow_mono_two h1 hκ hκ1 (Nat.sub_le_sub_right hM 1)
  have a3 : ((M : ℝ) + 1) ^ 3 ≤ ((M' : ℝ) + 1) ^ 3 := by gcongr
  have a4 : (w + η) ^ 2 / η ≤ (w' + η) ^ 2 / η := by gcongr
  have b1 : 0 ≤ (2 * κ₁) ^ (M - 1) := by positivity
  gcongr

theorem tailC_mono {κ₁ κ : ℝ} (h1 : 0 < κ₁) (hκ : κ₁ ≤ κ) (hκ1 : 1 ≤ κ) {M M' : ℕ}
    (hM : M ≤ M') {w w' ℓ s₀ : ℝ} (hw : 0 ≤ w) (hww : w ≤ w') (hℓ : 0 < ℓ) (hs₀ : 0 < s₀)
    (N : ℕ) : tailC κ₁ M w ℓ s₀ N ≤ tailC κ M' w' ℓ s₀ N := by
  unfold tailC
  have hA0 := A₀_pos
  have hB0 := B₀_pos
  have a1 := pow_mono_two h1 hκ hκ1 (Nat.sub_le_sub_right hM 1)
  have a2 := pow_mono_two h1 hκ hκ1 (Nat.mul_le_mul_left 4 (Nat.sub_le_sub_right hM 1))
  have a3 : ((M : ℝ) + 1 + 8 * A₀) ≤ ((M' : ℝ) + 1 + 8 * A₀) := by
    have : (M : ℝ) ≤ M' := by exact_mod_cast hM
    linarith
  have hw' : 0 ≤ w' := hw.trans hww
  gcongr

theorem norm_sum_sq_le_of_le {ι : Type*} (s : Finset ι) (z : ι → ℂ) (β X : ι → ℝ)
    (hβ : ∀ i ∈ s, 0 ≤ β i) (hX : ∀ i ∈ s, 0 ≤ X i) (h : ∀ i ∈ s, ‖z i‖ ^ 2 ≤ β i * X i) :
    ‖∑ i ∈ s, z i‖ ^ 2 ≤ (∑ i ∈ s, β i) * ∑ i ∈ s, X i := by
  have h1 : ‖∑ i ∈ s, z i‖ ≤ ∑ i ∈ s, Real.sqrt (β i) * Real.sqrt (X i) := by
    refine (norm_sum_le _ _).trans (Finset.sum_le_sum fun i hi => ?_)
    rw [← Real.sqrt_mul (hβ i hi)]
    exact Real.le_sqrt_of_sq_le (h i hi)
  have h2 := Finset.sum_mul_sq_le_sq_mul_sq s (fun i => Real.sqrt (β i)) (fun i => Real.sqrt (X i))
  have e1 : ∑ i ∈ s, Real.sqrt (β i) ^ 2 = ∑ i ∈ s, β i :=
    Finset.sum_congr rfl fun i hi => Real.sq_sqrt (hβ i hi)
  have e2 : ∑ i ∈ s, Real.sqrt (X i) ^ 2 = ∑ i ∈ s, X i :=
    Finset.sum_congr rfl fun i hi => Real.sq_sqrt (hX i hi)
  rw [e1, e2] at h2
  calc ‖∑ i ∈ s, z i‖ ^ 2 ≤ (∑ i ∈ s, Real.sqrt (β i) * Real.sqrt (X i)) ^ 2 :=
        pow_le_pow_left₀ (norm_nonneg _) h1 2
    _ ≤ _ := h2

/-- A group whose nodes all have |γ| ≥ H > s₀/2 lies on one side: |centre| ≥ H. -/
theorem groups_center (P : ExpSum) (T : Set ℝ) {s₀ H : ℝ} (hs₀ : 0 < s₀) (hH : s₀ < 2 * H)
    (hTH : ∀ x ∈ nodeSet P T, H ≤ |x|) {g : Finset ℝ} (hg : g ∈ groups P T s₀) :
    H ≤ |gC g| := by
  obtain ⟨x, hx, rfl⟩ := mem_groups.mp hg
  set g := grp (nodeSet P T) s₀ x
  have hne : g.Nonempty := ⟨x, self_mem_grp hx⟩
  have hmin := gMin_mem hne
  have hmax := gMax_mem hne
  have hHa := hTH _ (grp_subset x hmin)
  have hHb := hTH _ (grp_subset x hmax)
  have hab := gMin_le hmax
  by_cases hpos : 0 < gMin g
  · rw [abs_of_pos hpos] at hHa
    unfold gC; rw [abs_of_pos (by linarith)]; linarith
  · push Not at hpos
    by_cases hneg : gMax g < 0
    · rw [abs_of_neg hneg] at hHb
      unfold gC; rw [abs_of_neg (by linarith)]; linarith
    · push Not at hneg
      exfalso
      have h0 : (0 : ℝ) ∈ comp (nodeSet P T) s₀ x := Icc_subset_comp hs₀ hmin hmax ⟨hpos, hneg⟩
      have := comp_subset_balls hs₀ h0
      simp only [mem_iUnion] at this
      obtain ⟨γ, hγ, hγ0⟩ := this
      have hγH := hTH γ (grp_subset x hγ)
      have h1 := hγ0.1
      have h2 := hγ0.2
      have : |γ| < s₀ / 2 := abs_lt.mpr ⟨by linarith, by linarith⟩
      linarith

open Classical in
theorem part_M_pos (P : ExpSum) {g : Finset ℝ} (hg : g.Nonempty) (hsub : ∀ x ∈ g, ∃ j, P.node j = x)
    (hm : ∀ j, P.node j ∈ (g : Set ℝ) → 1 ≤ P.mult j) (c : ℝ) : 1 ≤ (P.part ↑g c).M := by
  obtain ⟨x, hx⟩ := hg
  obtain ⟨j₀, hj₀⟩ := hsub x hx
  have hj : P.node j₀ ∈ (g : Set ℝ) := by rw [hj₀]; exact_mod_cast hx
  have key : P.mult j₀ ≤ ∑ j, (if P.node j ∈ (g : Set ℝ) then P.mult j else 0) := by
    have := Finset.single_le_sum (f := fun j => if P.node j ∈ (g : Set ℝ) then P.mult j else 0)
      (fun j _ => Nat.zero_le _) (Finset.mem_univ j₀)
    simpa only [hj, ite_true] using this
  rw [ExpSum.part_M]
  refine (hm j₀ hj).trans (le_of_le_of_eq key (Finset.sum_congr rfl fun j _ => ?_))
  split_ifs <;> rfl

theorem integral_mid_le_grpQ (r : ExpSum) (c : ℝ) {ℓ : ℝ} (hℓ : 0 < ℓ) :
    ∫ x in Icc (ℓ / 4) (3 * ℓ / 4), ‖(r.dmod c).eval x‖ ^ 2 ≤ ∫ u : ℝ, ‖grpQ r c ℓ u‖ ^ 2 := by
  have e : ∫ x in Icc (ℓ / 4) (3 * ℓ / 4), ‖(r.dmod c).eval x‖ ^ 2 =
      ∫ x in Icc (ℓ / 4) (3 * ℓ / 4), ‖grpQ r c ℓ x‖ ^ 2 := by
    refine setIntegral_congr_fun measurableSet_Icc fun x hx => ?_
    rw [norm_grpQ, twin_eq_one hℓ hx, mul_one]
  rw [e]
  exact setIntegral_le_integral (grpQ_isTest r c hℓ).integrable_norm_sq
    (Eventually.of_forall fun x => by positivity)

/-- The derivative of the T-part splits into the group derivatives. -/
theorem part_dmod_eval_groups (P : ExpSum) (T : Set ℝ) (s₀ u : ℝ) :
    ∑ g ∈ groups P T s₀, cexp (I * gC g * u) * ((P.part ↑g (gC g)).dmod (gC g)).eval u =
      ((P.part T 0).dmod 0).eval u := by
  have hF : (fun v : ℝ => ∑ g ∈ groups P T s₀, cexp (I * gC g * v) * (P.part ↑g (gC g)).eval v) =
      (P.part T 0).eval := by
    funext v
    rw [← part_eval_groups (P := P) (T := T) s₀ v]
    exact Finset.sum_congr rfl fun g _ => P.part_eval_shift _ _ v
  have h1 : HasDerivAt (fun v : ℝ => ∑ g ∈ groups P T s₀, cexp (I * gC g * v) * (P.part ↑g (gC g)).eval v)
      (∑ g ∈ groups P T s₀, cexp (I * gC g * u) * ((P.part ↑g (gC g)).dmod (gC g)).eval u) u :=
    HasDerivAt.fun_sum fun g _ => (P.part ↑g (gC g)).hasDerivAt_mod_eval (gC g) u
  rw [hF] at h1
  exact h1.unique ((P.part T 0).hasDerivAt_eval u)

/-- **Theorem R″, high part.** -/
theorem highR (hT : Nazarov1993_Turan) : ∃ κ : ℝ, 1 ≤ κ ∧
    ∀ (P : ExpSum) (T : Set ℝ) (s₀ ℓ H η η' wmax : ℝ) (N Mmax : ℕ),
      0 < s₀ → 0 < ℓ → s₀ < 2 * H → 0 < η → 0 < η' → 1 ≤ N →
      (∀ j, P.node j ∈ T → H ≤ |P.node j|) → (∀ j, P.node j ∈ T → 1 ≤ P.mult j) →
      (∀ g ∈ groups P T s₀, (P.part ↑g (gC g)).M ≤ Mmax ∧ gW g ≤ wmax) →
      4 * A₀ / η * ((Mmax : ℝ) + 1) ^ 2 ≤ ℓ →
      (Mmax : ℝ) * (4 * A₀ / (s₀ / 4)) * ((Mmax : ℝ) + 1) ^ 2 ≤ ℓ →
      smallB κ Mmax wmax η ≤ H ^ 2 / 2 →
      (1 - ((1 + η'⁻¹) * ((groups P T s₀).card : ℝ) + 1) * (2 * Real.pi) *
          tailC κ Mmax wmax ℓ s₀ N) * ‖(P.part T 0).eval 0‖ ^ 2 ≤
        (∑ g ∈ groups P T s₀, 16 * (2 * κ) ^ (4 * ((P.part ↑g (gC g)).M - 1)) / (gC g) ^ 2) / ℓ *
          ((1 + η') * ∫ u : ℝ, ‖((P.part T 0).dmod 0).eval u * psiC ℓ u‖ ^ 2) := by
  obtain ⟨κB, hκB, hB⟩ := lemmaB hT
  obtain ⟨κT, hκT, hTl⟩ := tail_bound hT
  set κ := max (max κB κT) 1 with hκdef
  have hκ1 : 1 ≤ κ := le_max_right _ _
  have hκB' : κB ≤ κ := le_trans (le_max_left _ _) (le_max_left _ _)
  have hκT' : κT ≤ κ := le_trans (le_max_right _ _) (le_max_left _ _)
  refine ⟨κ, hκ1, fun P T s₀ ℓ H η η' wmax N Mmax hs₀ hℓ hsH hη hη' hN hTH hTm hgrp hbern htail
    hsmall => ?_⟩
  classical
  set 𝒢 := groups P T s₀
  have hH : 0 < H := by linarith
  have hA0 := A₀_pos
  have hTH' : ∀ x ∈ nodeSet P T, H ≤ |x| := by
    intro x hx
    obtain ⟨⟨j, hjx⟩, hj⟩ := mem_nodeSet.mp hx
    rw [← hjx] at hj ⊢
    exact hTH j hj
  have hg_ne : ∀ g ∈ 𝒢, g.Nonempty := fun g hg => groups_nonempty hg
  have hg_w : ∀ g ∈ 𝒢, 0 ≤ gW g := fun g hg => gW_nonneg (hg_ne g hg)
  have hg_nodes : ∀ g ∈ 𝒢, ∀ j, (P.part ↑g (gC g)).node j ∈ Icc (-gW g) (gW g) := by
    intro g hg
    exact P.part_node_mem (↑g) (hg_w g hg) fun j hj => abs_sub_gC_le (by exact_mod_cast hj)
  have hg_M : ∀ g ∈ 𝒢, 1 ≤ (P.part ↑g (gC g)).M := by
    intro g hg
    refine part_M_pos P (hg_ne g hg) (fun x hx => (mem_nodeSet.mp (groups_subset hg hx)).1)
      (fun j hj => ?_) _
    exact hTm j (mem_nodeSet.mp (groups_subset hg (by exact_mod_cast hj))).2
  have hg_C : ∀ g ∈ 𝒢, H ≤ |gC g| := fun g hg => groups_center P T hs₀ hsH hTH' hg
  have hg_C0 : ∀ g ∈ 𝒢, gC g ≠ 0 := by
    intro g hg h
    have := hg_C g hg
    rw [h, abs_zero] at this; linarith
  set Q : Finset ℝ → ℝ → ℂ := fun g => grpQ (P.part ↑g (gC g)) (gC g) ℓ
  have hQt : ∀ g, IsTest (Q g) := fun g => grpQ_isTest _ _ hℓ
  set X : Finset ℝ → ℝ := fun g => ∫ u : ℝ, ‖Q g u‖ ^ 2
  have hX0 : ∀ g, 0 ≤ X g := fun g => integral_nonneg fun u => by positivity
  set β : Finset ℝ → ℝ := fun g =>
    16 * (2 * κ) ^ (4 * ((P.part ↑g (gC g)).M - 1)) / (gC g) ^ 2 / ℓ
  have hβ0 : ∀ g, 0 ≤ β g := fun g => by simp only [β]; positivity
  -- Lemma B per group
  have hval : ∀ g ∈ 𝒢, ‖(P.part ↑g (gC g)).eval 0‖ ^ 2 ≤ β g * X g := by
    intro g hg
    obtain ⟨hMg, hwg⟩ := hgrp g hg
    have hsm : smallB κB (P.part ↑g (gC g)).M (gW g) η ≤ (gC g) ^ 2 / 2 := by
      have h1 := smallB_mono hκB hκB' hκ1 hMg (hg_w g hg) hwg hη
      have h2 : H ^ 2 ≤ (gC g) ^ 2 := by
        have h3 := hg_C g hg
        rw [← sq_abs (gC g)]
        exact pow_le_pow_left₀ hH.le h3 2
      linarith
    have hbern' : 4 * A₀ / η * (((P.part ↑g (gC g)).M : ℝ) + 1) ^ 2 ≤ ℓ := by
      have hc : ((P.part ↑g (gC g)).M : ℝ) ≤ Mmax := by exact_mod_cast hMg
      calc 4 * A₀ / η * (((P.part ↑g (gC g)).M : ℝ) + 1) ^ 2 ≤ 4 * A₀ / η * ((Mmax : ℝ) + 1) ^ 2 := by
            gcongr
        _ ≤ ℓ := hbern
    have h := hB (P.part ↑g (gC g)) (gC g) (gW g) η ℓ (hg_w g hg) hη hℓ (hg_M g hg) (hg_C0 g hg)
      (hg_nodes g hg) hbern' hsm 0 ⟨le_rfl, hℓ.le⟩
    have hmid := integral_mid_le_grpQ (P.part ↑g (gC g)) (gC g) hℓ
    set M := (P.part ↑g (gC g)).M
    have hc2 : 0 < (gC g) ^ 2 := by have := hg_C0 g hg; positivity
    have hp1 : ((2 * κB) ^ (M - 1)) ^ 2 * (2 * κB) ^ (2 * (M - 1)) ≤ (2 * κ) ^ (4 * (M - 1)) := by
      rw [← pow_mul, ← pow_add]
      exact pow_mono_two hκB hκB' hκ1 (by omega)
    have hcoef : 4 * ((2 * κB) ^ (M - 1)) ^ 2 * (4 * (2 * κB) ^ (2 * (M - 1))) / (ℓ * gC g ^ 2) ≤
        β g := by
      have e : β g = 16 * (2 * κ) ^ (4 * (M - 1)) / (ℓ * gC g ^ 2) := by
        simp only [β, M]; rw [div_div, mul_comm (gC g ^ 2) ℓ]
      rw [e]
      have e2 : 4 * ((2 * κB) ^ (M - 1)) ^ 2 * (4 * (2 * κB) ^ (2 * (M - 1))) =
          16 * (((2 * κB) ^ (M - 1)) ^ 2 * (2 * κB) ^ (2 * (M - 1))) := by ring
      rw [e2]
      gcongr
    calc ‖(P.part ↑g (gC g)).eval 0‖ ^ 2
        ≤ 4 * ((2 * κB) ^ (M - 1)) ^ 2 * (4 * (2 * κB) ^ (2 * (M - 1))) / (ℓ * gC g ^ 2) *
          ∫ x in Icc (ℓ / 4) (3 * ℓ / 4), ‖((P.part ↑g (gC g)).dmod (gC g)).eval x‖ ^ 2 := h
      _ ≤ β g * X g := by
          have hI0 : 0 ≤ ∫ x in Icc (ℓ / 4) (3 * ℓ / 4), ‖((P.part ↑g (gC g)).dmod (gC g)).eval x‖ ^ 2 :=
            setIntegral_nonneg measurableSet_Icc fun x _ => by positivity
          exact mul_le_mul hcoef hmid hI0 (hβ0 g)
  -- frequency tails per group
  have htailg : ∀ g ∈ 𝒢, ∫ ξ : ℝ, (farSet (gC g) (gW g) s₀).indicator
      (fun ξ => ‖kern (Q g) ξ‖ ^ 2) ξ ≤ tailC κ Mmax wmax ℓ s₀ N * X g := by
    intro g hg
    obtain ⟨hMg, hwg⟩ := hgrp g hg
    have hc : ((P.part ↑g (gC g)).M : ℝ) ≤ Mmax := by exact_mod_cast hMg
    have htail' : ((P.part ↑g (gC g)).M : ℝ) * (4 * A₀ / (s₀ / 4)) *
        (((P.part ↑g (gC g)).M : ℝ) + 1) ^ 2 ≤ ℓ := by
      calc ((P.part ↑g (gC g)).M : ℝ) * (4 * A₀ / (s₀ / 4)) * (((P.part ↑g (gC g)).M : ℝ) + 1) ^ 2
          ≤ (Mmax : ℝ) * (4 * A₀ / (s₀ / 4)) * ((Mmax : ℝ) + 1) ^ 2 := by gcongr
        _ ≤ ℓ := htail
    have h := hTl (P.part ↑g (gC g)) (gC g) (gW g) ℓ s₀ N (hg_w g hg) hs₀ hℓ (hg_M g hg) hN
      (hg_nodes g hg) htail'
    exact h.trans (mul_le_mul_of_nonneg_right
      (tailC_mono hκT hκT' hκ1 hMg (hg_w g hg) hwg hℓ hs₀ N) (hX0 g))
  -- Lemma D
  have hD := lemmaD (ι := 𝒢) (fun g => Q g.1) (fun g => hQt g.1) (fun g => gC g.1)
    (fun g => gW g.1) hs₀ (fun g g' hne => groups_sep hs₀ g.2 g'.2 (fun h => hne (Subtype.ext h)))
    (δ := tailC κ Mmax wmax ℓ s₀ N) (fun g => htailg g.1 g.2) hη'
  rw [Fintype.card_coe] at hD
  rw [Finset.sum_coe_sort 𝒢 (fun g => ∫ u : ℝ, ‖Q g u‖ ^ 2)] at hD
  have hsumQ : ∀ u, ∑ g : 𝒢, Q g.1 u = ((P.part T 0).dmod 0).eval u * psiC ℓ u := by
    intro u
    rw [Finset.sum_coe_sort 𝒢 (fun g => Q g u), ← part_dmod_eval_groups P T s₀ u, Finset.sum_mul]
    rfl
  simp_rw [hsumQ] at hD
  -- value at 0
  have hp0 : (P.part T 0).eval 0 = ∑ g ∈ 𝒢, (P.part ↑g (gC g)).eval 0 := by
    rw [← part_eval_groups (P := P) (T := T) s₀ 0]
    exact Finset.sum_congr rfl fun g _ => (P.part_eval_zero_eq _ _).symm
  have hCS := norm_sum_sq_le_of_le 𝒢 (fun g => (P.part ↑g (gC g)).eval 0) β X
    (fun g _ => hβ0 g) (fun g _ => hX0 g) hval
  rw [← hp0] at hCS
  set A := 1 - ((1 + η'⁻¹) * (𝒢.card : ℝ) + 1) * (2 * Real.pi) * tailC κ Mmax wmax ℓ s₀ N
  have hsumβ : ∑ g ∈ 𝒢, β g =
      (∑ g ∈ 𝒢, 16 * (2 * κ) ^ (4 * ((P.part ↑g (gC g)).M - 1)) / (gC g) ^ 2) / ℓ := by
    simp only [β]; rw [Finset.sum_div]
  rw [← hsumβ]
  have hRHS : 0 ≤ (∑ g ∈ 𝒢, β g) * ((1 + η') * ∫ u : ℝ, ‖((P.part T 0).dmod 0).eval u * psiC ℓ u‖ ^ 2) := by
    have : 0 ≤ ∫ u : ℝ, ‖((P.part T 0).dmod 0).eval u * psiC ℓ u‖ ^ 2 :=
      integral_nonneg fun u => by positivity
    exact mul_nonneg (Finset.sum_nonneg fun g _ => hβ0 g) (by positivity)
  rcases le_or_gt 0 A with hA | hA
  · calc A * ‖(P.part T 0).eval 0‖ ^ 2 ≤ A * ((∑ g ∈ 𝒢, β g) * ∑ g ∈ 𝒢, X g) :=
          mul_le_mul_of_nonneg_left hCS hA
      _ = (∑ g ∈ 𝒢, β g) * (A * ∑ g ∈ 𝒢, X g) := by ring
      _ ≤ (∑ g ∈ 𝒢, β g) * ((1 + η') * ∫ u : ℝ, ‖((P.part T 0).dmod 0).eval u * psiC ℓ u‖ ^ 2) :=
          mul_le_mul_of_nonneg_left hD (Finset.sum_nonneg fun g _ => hβ0 g)
  · have : A * ‖(P.part T 0).eval 0‖ ^ 2 ≤ 0 :=
      mul_nonpos_of_nonpos_of_nonneg hA.le (by positivity)
    linarith

end Riemann.Deep.LeanNBRate


