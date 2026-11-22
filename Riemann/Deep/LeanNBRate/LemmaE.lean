import Riemann.Deep.LeanNBRate.LowPart

/-!
# [N] Lemma E (without endpoint terms): low/high decoupling of the OU form

|B(p_lo, p_hi)| ≤ θ √(Σ_γ (γ²+¼)∫|r_γ|²) √(Σ_g ∫|Q_g|²), θ = 2√(σ/ℓ) (Σ_γ Σ_g Φ_g²/(γ - c_g)²)^{1/2}.
-/

open MeasureTheory Filter Topology Set Complex

noncomputable section

namespace Riemann.Deep.LeanNBRate

/-- Bernstein's factor. -/
def bernC (κ : ℝ) (M : ℕ) (w η : ℝ) : ℝ := (2 * κ) ^ (M - 1) * CB * ((M : ℝ) + 1) ^ 3 * ((w + η) ^ 2 / η)

theorem bernC_nonneg {κ : ℝ} (hκ : 0 < κ) (M : ℕ) {w η : ℝ} (hη : 0 < η) : 0 ≤ bernC κ M w η := by
  unfold bernC; have := CB_nonneg; positivity

/-- Sup bounds for one high group in terms of X = ∫|Q|². -/
theorem group_bounds (hT : Nazarov1993_Turan) : ∃ κ : ℝ, 1 ≤ κ ∧
    ∀ (r : ExpSum) (c w η ℓ : ℝ), 0 ≤ w → 0 < η → 0 < ℓ → 1 ≤ r.M → c ≠ 0 →
      (∀ j, r.node j ∈ Icc (-w) w) → 4 * A₀ / η * ((r.M : ℝ) + 1) ^ 2 ≤ ℓ →
      smallB κ r.M w η ≤ c ^ 2 / 2 →
      (∀ u ∈ Icc 0 ℓ, ‖r.eval u‖ ^ 2 ≤
        16 * (2 * κ) ^ (4 * (r.M - 1)) / (ℓ * c ^ 2) * ∫ x : ℝ, ‖grpQ r c ℓ x‖ ^ 2) ∧
      (∀ u ∈ Icc 0 ℓ, ‖(r.dmod c).eval u‖ ^ 2 ≤
        4 * (2 * κ) ^ (4 * (r.M - 1)) / ℓ * ∫ x : ℝ, ‖grpQ r c ℓ x‖ ^ 2) ∧
      (∀ K, (∀ u ∈ Icc 0 ℓ, ‖r.eval u‖ ≤ K) →
        ∀ u ∈ Icc 0 ℓ, ‖(r.dmod 0).eval u‖ ≤ bernC κ r.M w η * K) ∧
      (∀ K, (∀ u ∈ Icc 0 ℓ, ‖(r.dmod c).eval u‖ ≤ K) →
        ∀ u ∈ Icc 0 ℓ, ‖((r.dmod c).dmod 0).eval u‖ ≤ bernC κ r.M w η * K) := by
  obtain ⟨κB, hκB, hB⟩ := lemmaB hT
  obtain ⟨κS, hκS, hS⟩ := sup_sq_le hT
  obtain ⟨κE, hκE, hE⟩ := bernstein hT
  set κ := max (max κB κS) (max κE 1)
  have hκ1 : 1 ≤ κ := le_trans (le_max_right κE 1) (le_max_right _ _)
  have hκB' : κB ≤ κ := le_trans (le_max_left _ _) (le_max_left _ _)
  have hκS' : κS ≤ κ := le_trans (le_max_right _ _) (le_max_left _ _)
  have hκE' : κE ≤ κ := le_trans (le_max_left κE 1) (le_max_right _ _)
  refine ⟨κ, hκ1, fun r c w η ℓ hw hη hℓ hM hc hnodes hbern hsmall => ?_⟩
  set X := ∫ x : ℝ, ‖grpQ r c ℓ x‖ ^ 2
  have hmid := integral_mid_le_grpQ r c hℓ
  have hc2 : 0 < c ^ 2 := by positivity
  have hCB := CB_nonneg
  refine ⟨?_, ?_, ?_, ?_⟩
  · intro u hu
    have hsm : smallB κB r.M w η ≤ c ^ 2 / 2 :=
      (smallB_mono hκB hκB' hκ1 le_rfl hw le_rfl hη).trans hsmall
    have h := hB r c w η ℓ hw hη hℓ hM hc hnodes hbern hsm u hu
    have hp1 : ((2 * κB) ^ (r.M - 1)) ^ 2 * (2 * κB) ^ (2 * (r.M - 1)) ≤ (2 * κ) ^ (4 * (r.M - 1)) := by
      rw [← pow_mul, ← pow_add]
      exact pow_mono_two hκB hκB' hκ1 (by omega)
    have hI0 : 0 ≤ ∫ x in Icc (ℓ / 4) (3 * ℓ / 4), ‖(r.dmod c).eval x‖ ^ 2 :=
      setIntegral_nonneg measurableSet_Icc fun x _ => by positivity
    have hcoef : 4 * ((2 * κB) ^ (r.M - 1)) ^ 2 * (4 * (2 * κB) ^ (2 * (r.M - 1))) / (ℓ * c ^ 2) ≤
        16 * (2 * κ) ^ (4 * (r.M - 1)) / (ℓ * c ^ 2) := by
      have e2 : 4 * ((2 * κB) ^ (r.M - 1)) ^ 2 * (4 * (2 * κB) ^ (2 * (r.M - 1))) =
          16 * (((2 * κB) ^ (r.M - 1)) ^ 2 * (2 * κB) ^ (2 * (r.M - 1))) := by ring
      rw [e2]; gcongr
    calc ‖r.eval u‖ ^ 2 ≤ 4 * ((2 * κB) ^ (r.M - 1)) ^ 2 * (4 * (2 * κB) ^ (2 * (r.M - 1))) / (ℓ * c ^ 2) *
          ∫ x in Icc (ℓ / 4) (3 * ℓ / 4), ‖(r.dmod c).eval x‖ ^ 2 := h
      _ ≤ 16 * (2 * κ) ^ (4 * (r.M - 1)) / (ℓ * c ^ 2) * X :=
          mul_le_mul hcoef hmid hI0 (by positivity)
  · intro u hu
    have h := hS (r.dmod c) ℓ hℓ (by rw [r.dmod_M]; exact hM) u hu
    rw [r.dmod_M] at h
    have hp : (2 * κS) ^ (4 * (r.M - 1)) ≤ (2 * κ) ^ (4 * (r.M - 1)) :=
      pow_le_pow_left₀ (by positivity) (by linarith) _
    have hI0 : 0 ≤ ∫ x in Icc (ℓ / 4) (3 * ℓ / 4), ‖(r.dmod c).eval x‖ ^ 2 :=
      setIntegral_nonneg measurableSet_Icc fun x _ => by positivity
    calc ‖(r.dmod c).eval u‖ ^ 2 ≤ 4 * (2 * κS) ^ (4 * (r.M - 1)) / ℓ *
          ∫ x in Icc (ℓ / 4) (3 * ℓ / 4), ‖(r.dmod c).eval x‖ ^ 2 := h
      _ ≤ 4 * (2 * κ) ^ (4 * (r.M - 1)) / ℓ * X := by
          apply mul_le_mul _ hmid hI0 (by positivity)
          gcongr
  · intro K hK u hu
    have h := hE r w η ℓ K hw hη hℓ hM hnodes hbern hK u hu
    have hK0 : 0 ≤ K := (norm_nonneg _).trans (hK 0 ⟨le_rfl, hℓ.le⟩)
    refine h.trans ?_
    unfold bernC
    have hp : (2 * κE) ^ (r.M - 1) ≤ (2 * κ) ^ (r.M - 1) := pow_le_pow_left₀ (by positivity) (by linarith) _
    gcongr
  · intro K hK u hu
    have hnq : ∀ j, (r.dmod c).node j ∈ Icc (-w) w := fun j => by rw [r.dmod_node]; exact hnodes j
    have h := hE (r.dmod c) w η ℓ K hw hη hℓ (by rw [r.dmod_M]; exact hM) hnq
      (by rw [r.dmod_M]; exact hbern) hK u hu
    rw [r.dmod_M] at h
    have hK0 : 0 ≤ K := (norm_nonneg _).trans (hK 0 ⟨le_rfl, hℓ.le⟩)
    refine h.trans ?_
    unfold bernC
    have hp : (2 * κE) ^ (r.M - 1) ≤ (2 * κ) ^ (r.M - 1) := pow_le_pow_left₀ (by positivity) (by linarith) _
    gcongr

open ComplexConjugate in
theorem oubil_sum2 {ι κ' : Type*} (s : Finset ι) (t : Finset κ') (f f' : ι → ℝ → ℂ) (g g' : κ' → ℝ → ℂ)
    (hf : ∀ i, Continuous (f i)) (hf' : ∀ i, Continuous (f' i))
    (hg : ∀ j, Continuous (g j)) (hg' : ∀ j, Continuous (g' j)) (ℓ : ℝ) :
    oubil ℓ (fun u => ∑ i ∈ s, f i u) (fun u => ∑ i ∈ s, f' i u) (fun u => ∑ j ∈ t, g j u)
        (fun u => ∑ j ∈ t, g' j u) =
      ∑ i ∈ s, ∑ j ∈ t, oubil ℓ (f i) (f' i) (g j) (g' j) := by
  unfold oubil
  have e : ∀ u : ℝ, (∑ i ∈ s, f' i u) * conj (∑ j ∈ t, g' j u) +
      (1 / 4 : ℂ) * ((∑ i ∈ s, f i u) * conj (∑ j ∈ t, g j u)) =
      ∑ i ∈ s, ∑ j ∈ t, (f' i u * conj (g' j u) + (1 / 4 : ℂ) * (f i u * conj (g j u))) := by
    intro u
    rw [map_sum, map_sum, Finset.sum_mul_sum, Finset.sum_mul_sum, Finset.mul_sum,
      ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [Finset.mul_sum, ← Finset.sum_add_distrib]
  simp_rw [e]
  rw [intervalIntegral.integral_finsetSum]
  · refine Finset.sum_congr rfl fun i _ => ?_
    rw [intervalIntegral.integral_finsetSum]
    intro j _
    exact (((hf' i).mul ((hg' j).star)).add (continuous_const.mul ((hf i).mul (hg j).star))).intervalIntegrable _ _
  · intro i _
    exact (continuous_finsetSum _ fun j _ => (((hf' i).mul ((hg' j).star)).add
      (continuous_const.mul ((hf i).mul (hg j).star)))).intervalIntegrable _ _

theorem sum_mul_le_sqrt {ι : Type*} (s : Finset ι) (a b : ι → ℝ) (ha : ∀ i ∈ s, 0 ≤ a i)
    (hb : ∀ i ∈ s, 0 ≤ b i) :
    ∑ i ∈ s, a i * b i ≤ Real.sqrt (∑ i ∈ s, a i ^ 2) * Real.sqrt (∑ i ∈ s, b i ^ 2) := by
  have h := Finset.sum_mul_sq_le_sq_mul_sq s a b
  have h0 : 0 ≤ ∑ i ∈ s, a i * b i := Finset.sum_nonneg fun i hi => mul_nonneg (ha i hi) (hb i hi)
  rw [← Real.sqrt_mul (Finset.sum_nonneg fun i _ => sq_nonneg _)]
  exact Real.le_sqrt_of_sq_le h

/-- The Lemma E weight of one high group. -/
def phiE (κ : ℝ) (H₀ mk ℓ : ℝ) (M : ℕ) (w η : ℝ) : ℝ :=
  3 * (H₀ + mk + 1) * (2 + mk + ℓ * bernC κ M w η) * Real.sqrt (4 * (2 * κ) ^ (4 * (M - 1)) / ℓ)

open ComplexConjugate in
/-- **[N] Lemma E** (no endpoint terms): the cross form of the low part and the high part. -/
theorem lemmaE (hT : Nazarov1993_Turan) : ∃ κ c : ℝ, 1 ≤ κ ∧ 1 ≤ c ∧
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
          Real.sqrt (∑ g ∈ groups P Thi s₀, ∫ x : ℝ, ‖grpQ (P.part ↑g (gC g)) (gC g) ℓ x‖ ^ 2) := by
  obtain ⟨κ, hκ1, hG⟩ := group_bounds hT
  obtain ⟨c₀, hc₀, hPB⟩ := poly_bounds hT
  set c := max c₀ 1
  have hc1 : 1 ≤ c := le_max_right _ _
  have hcc : c₀ ≤ c := le_max_left _ _
  refine ⟨κ, c, hκ1, hc1, fun P Tlo Thi s₀ ℓ η H₀ M₀ hℓ hη hH₀ hH hM hgrp hne => ?_⟩
  classical
  have hℓ0 : 0 < ℓ := by linarith
  set S := nodeSet P Tlo
  set 𝒢 := groups P Thi s₀
  set r : ℝ → ExpSum := fun γ => P.part {γ} γ
  set rg : Finset ℝ → ExpSum := fun g => P.part ↑g (gC g)
  set Ig : ℝ → ℝ := fun γ => ∫ u in (0 : ℝ)..ℓ, ‖(r γ).eval u‖ ^ 2
  set Sg : ℝ → ℝ := fun γ => supS c₀ ℓ (r γ)
  set X : Finset ℝ → ℝ := fun g => ∫ x : ℝ, ‖grpQ (rg g) (gC g) ℓ x‖ ^ 2
  set σ := 2 * (2 * c) ^ (2 * (M₀ - 1))
  set mk := markovC c M₀
  have hmk0 : 0 ≤ mk := markovC_nonneg (by linarith) M₀
  have hIg0 : ∀ γ, 0 ≤ Ig γ := fun γ => intervalIntegral.integral_nonneg hℓ0.le fun u _ => by positivity
  have hSg0 : ∀ γ, 0 ≤ Sg γ := fun γ => supS_nonneg _ _ _
  have hX0 : ∀ g, 0 ≤ X g := fun g => integral_nonneg fun u => by positivity
  have hSg2 : ∀ γ ∈ S, Sg γ ^ 2 ≤ σ / ℓ * Ig γ := by
    intro γ hγ
    simp only [Sg]
    rw [supS_sq hc₀ hℓ0]
    have hp := pow_mono_two hc₀ hcc hc1 (Nat.mul_le_mul_left 2 (Nat.sub_le_sub_right (hM γ hγ).2 1))
    have := hIg0 γ
    simp only [σ, Ig]
    gcongr
  have hmkg : ∀ γ ∈ S, markovC c₀ (r γ).M ≤ mk := fun γ hγ =>
    markovC_mono2 hc₀ hcc hc1 (hM γ hγ).2
  have hB : ∀ γ ∈ S, (∀ u ∈ Icc 0 ℓ, ‖(r γ).eval u‖ ≤ Sg γ) ∧
      (∀ u ∈ Icc 0 ℓ, ‖((r γ).dmod 0).eval u‖ ≤ mk / ℓ * Sg γ) ∧
      (∀ u ∈ Icc 0 ℓ, ‖((r γ).dmod γ).eval u‖ ≤ (H₀ + mk) * Sg γ) ∧
      (∀ u ∈ Icc 0 ℓ, ‖(((r γ).dmod γ).dmod 0).eval u‖ ≤ mk / ℓ * ((H₀ + mk) * Sg γ)) := by
    intro γ hγ
    obtain ⟨b1, b2, b3, b4⟩ := hPB (r γ) ℓ γ (part_singleton_isPoly P γ) hℓ0 (hM γ hγ).1
    have hmg := hmkg γ hγ
    have hq : markovC c₀ (r γ).M / ℓ ≤ mk / ℓ := by gcongr
    have hq2 : mk / ℓ ≤ mk := div_le_self hmk0 hℓ
    have hγH := hH γ hγ
    have hS0 := hSg0 γ
    have hq0 : 0 ≤ markovC c₀ (r γ).M / ℓ := div_nonneg (markovC_nonneg hc₀ _) hℓ0.le
    refine ⟨b1, fun u hu => (b2 u hu).trans (by gcongr), fun u hu => (b3 u hu).trans ?_,
      fun u hu => (b4 u hu).trans ?_⟩
    · gcongr; linarith
    · have : 0 ≤ (|γ| + markovC c₀ (r γ).M / ℓ) * Sg γ := by positivity
      calc markovC c₀ (r γ).M / ℓ * ((|γ| + markovC c₀ (r γ).M / ℓ) * Sg γ)
          ≤ mk / ℓ * ((|γ| + markovC c₀ (r γ).M / ℓ) * Sg γ) := by gcongr
        _ ≤ mk / ℓ * ((H₀ + mk) * Sg γ) := by gcongr; linarith
  -- the functions
  set f : ℝ → ℝ → ℂ := fun γ u => cexp (I * γ * u) * (r γ).eval u
  set f' : ℝ → ℝ → ℂ := fun γ u => cexp (I * γ * u) * ((r γ).dmod γ).eval u
  set h : Finset ℝ → ℝ → ℂ := fun g u => cexp (I * gC g * u) * (rg g).eval u
  set h' : Finset ℝ → ℝ → ℂ := fun g u => cexp (I * gC g * u) * ((rg g).dmod (gC g)).eval u
  have hfc : ∀ γ, Continuous (f γ) := fun γ => by
    have := (r γ).continuous_eval; simp only [f]; fun_prop
  have hf'c : ∀ γ, Continuous (f' γ) := fun γ => by
    have := ((r γ).dmod γ).continuous_eval; simp only [f']; fun_prop
  have hhc : ∀ g, Continuous (h g) := fun g => by
    have := (rg g).continuous_eval; simp only [h]; fun_prop
  have hh'c : ∀ g, Continuous (h' g) := fun g => by
    have := ((rg g).dmod (gC g)).continuous_eval; simp only [h']; fun_prop
  have hp : (P.part Tlo 0).eval = fun u => ∑ γ ∈ S, f γ u := by
    funext u; exact (part_eval_singletons P Tlo u).symm
  have hp' : ((P.part Tlo 0).dmod 0).eval = fun u => ∑ γ ∈ S, f' γ u := by
    funext u; exact (part_dmod_eval_singletons P Tlo u).symm
  have hq : (P.part Thi 0).eval = fun u => ∑ g ∈ 𝒢, h g u := by
    funext u
    rw [← part_eval_groups (P := P) (T := Thi) s₀ u]
    exact Finset.sum_congr rfl fun g _ => (P.part_eval_shift _ _ u).symm
  have hq' : ((P.part Thi 0).dmod 0).eval = fun u => ∑ g ∈ 𝒢, h' g u := by
    funext u; exact (part_dmod_eval_groups P Thi s₀ u).symm
  rw [hp, hp', hq, hq', oubil_sum2 S 𝒢 f f' h h' hfc hf'c hhc hh'c ℓ]
  -- per-pair bound
  set Φ : Finset ℝ → ℝ := fun g => phiE κ H₀ mk ℓ (rg g).M (gW g) η
  have hpair : ∀ γ ∈ S, ∀ g ∈ 𝒢, ‖oubil ℓ (f γ) (f' γ) (h g) (h' g)‖ ≤
      Sg γ * (Real.sqrt (X g) * (Φ g / |γ - gC g|)) := by
    intro γ hγ g hg
    obtain ⟨b1, b2, b3, b4⟩ := hB γ hγ
    obtain ⟨gM, gC1, gW0, gnodes, gbern, gsmall⟩ := hgrp g hg
    have hc0 : gC g ≠ 0 := by intro h0; rw [h0, abs_zero] at gC1; linarith
    obtain ⟨G1, G2, G3, G4⟩ := hG (rg g) (gC g) (gW g) η ℓ gW0 hη hℓ0 gM hc0 gnodes gbern gsmall
    set E := (2 * κ) ^ (4 * ((rg g).M - 1))
    have hE0 : 0 ≤ E := by positivity
    set ag := Real.sqrt (4 * E / ℓ * X g)
    set sg := Real.sqrt (16 * E / (ℓ * gC g ^ 2) * X g)
    have hag : ∀ u ∈ Icc 0 ℓ, ‖((rg g).dmod (gC g)).eval u‖ ≤ ag := fun u hu =>
      Real.le_sqrt_of_sq_le (G2 u hu)
    have hsg : ∀ u ∈ Icc 0 ℓ, ‖(rg g).eval u‖ ≤ sg := fun u hu => Real.le_sqrt_of_sq_le (G1 u hu)
    set be := bernC κ (rg g).M (gW g) η
    have hbe0 : 0 ≤ be := bernC_nonneg (by linarith) _ hη
    have hcb := cross_bound (r γ) (rg g) hℓ0.le (hne γ hγ g hg) b3 hag b4 (G4 ag hag) b1 hsg b2 (G3 sg hsg)
    have hag0 : 0 ≤ ag := Real.sqrt_nonneg _
    have hsg0 : 0 ≤ sg := Real.sqrt_nonneg _
    have hS0 := hSg0 γ
    -- ag + sg ≤ 3 √(4E/ℓ) √X
    have hsum : ag + sg ≤ 3 * Real.sqrt (4 * E / ℓ) * Real.sqrt (X g) := by
      have e1 : ag = Real.sqrt (4 * E / ℓ) * Real.sqrt (X g) := Real.sqrt_mul (by positivity) _
      have e2 : sg = 2 / |gC g| * (Real.sqrt (4 * E / ℓ) * Real.sqrt (X g)) := by
        simp only [sg]
        rw [show 16 * E / (ℓ * gC g ^ 2) * X g = (2 / |gC g|) ^ 2 * (4 * E / ℓ * X g) by
          rw [div_pow, sq_abs]; field_simp; ring]
        rw [Real.sqrt_mul (by positivity), Real.sqrt_sq (by positivity), Real.sqrt_mul (by positivity)]
      have h2 : 2 / |gC g| ≤ 2 := by
        rw [div_le_iff₀ (by linarith)]; nlinarith
      have h3 : 0 ≤ Real.sqrt (4 * E / ℓ) * Real.sqrt (X g) := by positivity
      rw [e1, e2]
      nlinarith
    have hnum : 2 * ((H₀ + mk) * Sg γ * ag + Sg γ * sg / 4) +
        ℓ * (mk / ℓ * ((H₀ + mk) * Sg γ) * ag + (H₀ + mk) * Sg γ * (be * ag) +
          (mk / ℓ * Sg γ * sg + Sg γ * (be * sg)) / 4) ≤
        Sg γ * ((H₀ + mk + 1) * (2 + mk + ℓ * be) * (ag + sg)) := by
      have e : 2 * ((H₀ + mk) * Sg γ * ag + Sg γ * sg / 4) +
          ℓ * (mk / ℓ * ((H₀ + mk) * Sg γ) * ag + (H₀ + mk) * Sg γ * (be * ag) +
            (mk / ℓ * Sg γ * sg + Sg γ * (be * sg)) / 4) =
          Sg γ * ((H₀ + mk) * (2 + mk + ℓ * be) * ag + (1 / 2 + mk / 4 + ℓ * be / 4) * sg) := by
        field_simp; ring
      rw [e]
      apply mul_le_mul_of_nonneg_left _ hS0
      have hlb : 0 ≤ ℓ * be := by positivity
      have k1 : (H₀ + mk) * (2 + mk + ℓ * be) * ag ≤ (H₀ + mk + 1) * (2 + mk + ℓ * be) * ag := by
        have h1 : H₀ + mk ≤ H₀ + mk + 1 := by linarith
        have h2 : 0 ≤ 2 + mk + ℓ * be := by positivity
        exact mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_right h1 h2) hag0
      have k2 : (1 / 2 + mk / 4 + ℓ * be / 4) * sg ≤ (H₀ + mk + 1) * (2 + mk + ℓ * be) * sg := by
        have : 1 / 2 + mk / 4 + ℓ * be / 4 ≤ 1 * (2 + mk + ℓ * be) := by linarith
        have h2 : 0 ≤ 2 + mk + ℓ * be := by positivity
        have h3 : 1 ≤ H₀ + mk + 1 := by linarith
        calc (1 / 2 + mk / 4 + ℓ * be / 4) * sg ≤ (1 * (2 + mk + ℓ * be)) * sg :=
              mul_le_mul_of_nonneg_right this hsg0
          _ ≤ (H₀ + mk + 1) * (2 + mk + ℓ * be) * sg :=
              mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_right h3 h2) hsg0
      nlinarith
    have hpos : 0 < |γ - gC g| := abs_pos.mpr (sub_ne_zero.mpr (hne γ hγ g hg))
    have hΦ : (H₀ + mk + 1) * (2 + mk + ℓ * be) * (ag + sg) ≤ Real.sqrt (X g) * Φ g := by
      simp only [Φ, phiE]
      have : 0 ≤ (H₀ + mk + 1) * (2 + mk + ℓ * be) := by positivity
      calc (H₀ + mk + 1) * (2 + mk + ℓ * be) * (ag + sg)
          ≤ (H₀ + mk + 1) * (2 + mk + ℓ * be) * (3 * Real.sqrt (4 * E / ℓ) * Real.sqrt (X g)) := by gcongr
        _ = _ := by ring
    calc ‖oubil ℓ (f γ) (f' γ) (h g) (h' g)‖
        ≤ (2 * ((H₀ + mk) * Sg γ * ag + Sg γ * sg / 4) +
          ℓ * (mk / ℓ * ((H₀ + mk) * Sg γ) * ag + (H₀ + mk) * Sg γ * (be * ag) +
            (mk / ℓ * Sg γ * sg + Sg γ * (be * sg)) / 4)) / |γ - gC g| := hcb
      _ ≤ Sg γ * (Real.sqrt (X g) * Φ g) / |γ - gC g| := by
          gcongr
          exact hnum.trans (mul_le_mul_of_nonneg_left hΦ hS0)
      _ = _ := by ring
  -- summation
  have hΦ0 : ∀ g, 0 ≤ Φ g := fun g => by
    simp only [Φ, phiE]
    have := bernC_nonneg (by linarith : (0 : ℝ) < κ) (rg g).M (w := gW g) hη
    positivity
  set K : ℝ → ℝ := fun γ => ∑ g ∈ 𝒢, Φ g ^ 2 / (γ - gC g) ^ 2
  have hK0 : ∀ γ, 0 ≤ K γ := fun γ => Finset.sum_nonneg fun g _ => by positivity
  have hinner : ∀ γ ∈ S, ∑ g ∈ 𝒢, Real.sqrt (X g) * (Φ g / |γ - gC g|) ≤
      Real.sqrt (∑ g ∈ 𝒢, X g) * Real.sqrt (K γ) := by
    intro γ _
    have h1 := sum_mul_le_sqrt 𝒢 (fun g => Real.sqrt (X g)) (fun g => Φ g / |γ - gC g|)
      (fun g _ => Real.sqrt_nonneg _) (fun g _ => div_nonneg (hΦ0 g) (abs_nonneg _))
    have e1 : ∑ g ∈ 𝒢, Real.sqrt (X g) ^ 2 = ∑ g ∈ 𝒢, X g :=
      Finset.sum_congr rfl fun g _ => Real.sq_sqrt (hX0 g)
    have e2 : ∑ g ∈ 𝒢, (Φ g / |γ - gC g|) ^ 2 = K γ := by
      simp only [K]
      exact Finset.sum_congr rfl fun g _ => by rw [div_pow, sq_abs]
    rw [e1, e2] at h1
    exact h1
  have hsum1 : ‖∑ γ ∈ S, ∑ g ∈ 𝒢, oubil ℓ (f γ) (f' γ) (h g) (h' g)‖ ≤
      ∑ γ ∈ S, Sg γ * (Real.sqrt (∑ g ∈ 𝒢, X g) * Real.sqrt (K γ)) := by
    refine (norm_sum_le _ _).trans (Finset.sum_le_sum fun γ hγ => ?_)
    refine (norm_sum_le _ _).trans ?_
    calc ∑ g ∈ 𝒢, ‖oubil ℓ (f γ) (f' γ) (h g) (h' g)‖
        ≤ ∑ g ∈ 𝒢, Sg γ * (Real.sqrt (X g) * (Φ g / |γ - gC g|)) :=
          Finset.sum_le_sum fun g hg => hpair γ hγ g hg
      _ = Sg γ * ∑ g ∈ 𝒢, Real.sqrt (X g) * (Φ g / |γ - gC g|) := by rw [Finset.mul_sum]
      _ ≤ Sg γ * (Real.sqrt (∑ g ∈ 𝒢, X g) * Real.sqrt (K γ)) :=
          mul_le_mul_of_nonneg_left (hinner γ hγ) (hSg0 γ)
  have hsum2 : ∑ γ ∈ S, Sg γ * Real.sqrt (K γ) ≤
      Real.sqrt (∑ γ ∈ S, Sg γ ^ 2) * Real.sqrt (∑ γ ∈ S, K γ) := by
    have h1 := sum_mul_le_sqrt S Sg (fun γ => Real.sqrt (K γ)) (fun γ _ => hSg0 γ)
      (fun γ _ => Real.sqrt_nonneg _)
    have e : ∑ γ ∈ S, Real.sqrt (K γ) ^ 2 = ∑ γ ∈ S, K γ :=
      Finset.sum_congr rfl fun γ _ => Real.sq_sqrt (hK0 γ)
    rw [e] at h1; exact h1
  set D := ∑ γ ∈ S, (γ ^ 2 + 1 / 4) * Ig γ
  have hsumI : ∑ γ ∈ S, Ig γ ≤ 4 * D := by
    simp only [D, Finset.mul_sum]
    refine Finset.sum_le_sum fun γ _ => ?_
    have := hIg0 γ
    nlinarith [sq_nonneg γ]
  have hS2 : Real.sqrt (∑ γ ∈ S, Sg γ ^ 2) ≤ 2 * Real.sqrt (σ / ℓ) * Real.sqrt D := by
    have h1 : ∑ γ ∈ S, Sg γ ^ 2 ≤ σ / ℓ * (4 * D) := by
      calc ∑ γ ∈ S, Sg γ ^ 2 ≤ ∑ γ ∈ S, σ / ℓ * Ig γ := Finset.sum_le_sum hSg2
        _ = σ / ℓ * ∑ γ ∈ S, Ig γ := by rw [Finset.mul_sum]
        _ ≤ σ / ℓ * (4 * D) := by
            have : 0 ≤ σ / ℓ := by positivity
            gcongr
    have hD0 : 0 ≤ D := Finset.sum_nonneg fun γ _ => by have := hIg0 γ; positivity
    calc Real.sqrt (∑ γ ∈ S, Sg γ ^ 2) ≤ Real.sqrt (σ / ℓ * (4 * D)) := Real.sqrt_le_sqrt h1
      _ = 2 * Real.sqrt (σ / ℓ) * Real.sqrt D := by
          rw [Real.sqrt_mul (by positivity), Real.sqrt_mul (by norm_num),
            show Real.sqrt 4 = 2 by rw [show (4 : ℝ) = 2 ^ 2 by norm_num, Real.sqrt_sq (by norm_num)]]
          ring
  have hY0 : 0 ≤ Real.sqrt (∑ g ∈ 𝒢, X g) := Real.sqrt_nonneg _
  calc ‖∑ γ ∈ S, ∑ g ∈ 𝒢, oubil ℓ (f γ) (f' γ) (h g) (h' g)‖
      ≤ ∑ γ ∈ S, Sg γ * (Real.sqrt (∑ g ∈ 𝒢, X g) * Real.sqrt (K γ)) := hsum1
    _ = Real.sqrt (∑ g ∈ 𝒢, X g) * ∑ γ ∈ S, Sg γ * Real.sqrt (K γ) := by
        rw [Finset.mul_sum]; refine Finset.sum_congr rfl fun γ _ => ?_; ring
    _ ≤ Real.sqrt (∑ g ∈ 𝒢, X g) * (Real.sqrt (∑ γ ∈ S, Sg γ ^ 2) * Real.sqrt (∑ γ ∈ S, K γ)) :=
        mul_le_mul_of_nonneg_left hsum2 hY0
    _ ≤ Real.sqrt (∑ g ∈ 𝒢, X g) * ((2 * Real.sqrt (σ / ℓ) * Real.sqrt D) * Real.sqrt (∑ γ ∈ S, K γ)) := by
        gcongr
    _ = _ := by simp only [K]; ring

end Riemann.Deep.LeanNBRate


