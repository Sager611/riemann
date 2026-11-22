import Riemann.Deep.LeanNBRate.GroupSums

/-!
# U2: group sums of a zero configuration
-/

open MeasureTheory Filter Topology Set Complex ComplexConjugate

noncomputable section

namespace Riemann.Deep.LeanNBRate

theorem sum_groups_sum (P : ExpSum) (T : Set ℝ) (s₀ : ℝ) (f : ℝ → ℝ) :
    ∑ g ∈ groups P T s₀, ∑ γ ∈ g, f γ = ∑ γ ∈ nodeSet P T, f γ := by
  classical
  have hU : (groups P T s₀).biUnion id = nodeSet P T := by
    ext x
    simp only [Finset.mem_biUnion, id]
    constructor
    · rintro ⟨g, hg, hx⟩; exact groups_subset hg hx
    · intro hx; exact groups_cover hx
  rw [← hU, Finset.sum_biUnion]
  · rfl
  · intro g hg g' hg' hne
    show Disjoint g g'
    rw [Finset.disjoint_left]
    intro x hx hx'
    exact hne (groups_disj hg hg' hx hx')

/-- (2κ)^{4M} ≤ √X from M ≤ ε₁ log X. -/
theorem pow_le_sqrt_of_le_log {κ ε₁ X : ℝ} {M : ℕ} (hκ : 1 ≤ κ) (hε₁ : 0 ≤ ε₁) (hX : 1 ≤ X)
    (hsmall : 4 * ε₁ * Real.log (2 * κ) ≤ 1 / 2) (hM : (M : ℝ) ≤ ε₁ * Real.log X) :
    (2 * κ) ^ (4 * M) ≤ Real.sqrt X := by
  have h2κ : 0 < 2 * κ := by linarith
  have hl0 : 0 ≤ Real.log (2 * κ) := Real.log_nonneg (by linarith)
  have hlX : 0 ≤ Real.log X := Real.log_nonneg hX
  rw [Real.sqrt_eq_rpow, ← Real.exp_log (pow_pos h2κ _), ← Real.exp_log (by linarith : (0 : ℝ) < X),
    ← Real.exp_mul]
  apply Real.exp_le_exp.mpr
  rw [Real.log_pow]
  push_cast
  calc (4 * (M : ℝ)) * Real.log (2 * κ) ≤ 4 * (ε₁ * Real.log X) * Real.log (2 * κ) := by gcongr
    _ = (4 * ε₁ * Real.log (2 * κ)) * Real.log X := by ring
    _ ≤ 1 / 2 * Real.log X := mul_le_mul_of_nonneg_right hsmall hlX
    _ = Real.log X * (1 / 2) := by ring

theorem F1_of_ge_one {t : ℝ} (ht : 1 ≤ t) : F1 t = 4 * Real.sqrt 2 * Real.sqrt (t + 1) / t ^ 2 := by
  unfold F1; rw [max_eq_left ht]

/-- A group centre near a node: √(|c|+2)/c² ≤ F₁(|γ|). -/
theorem sqrt_div_sq_le_F1 {c γ : ℝ} (hγ : 2 ≤ |γ|) (hc : |γ - c| ≤ 1) :
    Real.sqrt (|c| + 2) / c ^ 2 ≤ F1 |γ| := by
  rw [F1_of_ge_one (by linarith)]
  have hc1 : |c| ≤ |γ| + 1 := by
    have := abs_sub_abs_le_abs_sub c γ
    rw [abs_sub_comm] at hc
    linarith
  have hc2 : |γ| - 1 ≤ |c| := by
    have := abs_sub_abs_le_abs_sub γ c
    linarith
  have hcpos : 0 < |c| := by linarith
  have hγpos : 0 < |γ| := by linarith
  have e1 : Real.sqrt (|c| + 2) ≤ Real.sqrt 2 * Real.sqrt (|γ| + 1) := by
    rw [← Real.sqrt_mul (by norm_num)]
    apply Real.sqrt_le_sqrt
    linarith
  have e2 : 1 / c ^ 2 ≤ 4 / γ ^ 2 := by
    rw [← sq_abs c, ← sq_abs γ, div_le_div_iff₀ (by positivity) (by positivity)]
    have : |γ| / 2 ≤ |c| := by linarith
    nlinarith
  have hs0 : 0 ≤ Real.sqrt (|c| + 2) := Real.sqrt_nonneg _
  calc Real.sqrt (|c| + 2) / c ^ 2 = Real.sqrt (|c| + 2) * (1 / c ^ 2) := by ring
    _ ≤ (Real.sqrt 2 * Real.sqrt (|γ| + 1)) * (4 / γ ^ 2) := by
        apply mul_le_mul e1 e2 (by positivity) (by positivity)
    _ = 4 * Real.sqrt 2 * Real.sqrt (|γ| + 1) / |γ| ^ 2 := by rw [sq_abs]; ring

/-- One group: (2κ)^{4(M_g-1)}/c_g² ≤ Σ_{γ∈g} nw(γ)·F₁(|γ|). -/
theorem group_weight_le {P : ExpSum} {T : ℝ} (hP : ZeroConfig P T) {Th : Set ℝ} {s₀ κ ε₁ : ℝ}
    {g : Finset ℝ} (hg : g ∈ groups P Th s₀) (hκ : 1 ≤ κ) (hε₁ : 0 ≤ ε₁)
    (hsmall : 4 * ε₁ * Real.log (2 * κ) ≤ 1 / 2)
    (hM : ((P.part ↑g (gC g)).M : ℝ) ≤ ε₁ * Real.log (|gC g| + 2)) (hW : gW g ≤ 1)
    (hTh : ∀ γ ∈ nodeSet P Th, 2 ≤ |γ|) :
    (2 * κ) ^ (4 * ((P.part ↑g (gC g)).M - 1)) / (gC g) ^ 2 ≤ ∑ γ ∈ g, (nw P γ : ℝ) * F1 |γ| := by
  have hne := groups_nonempty hg
  obtain ⟨γ₀, hγ₀⟩ := hne
  have hγ₀N := groups_subset hg hγ₀
  have hγ₀2 := hTh γ₀ hγ₀N
  have hd : |γ₀ - gC g| ≤ 1 := (abs_sub_gC_le hγ₀).trans hW
  have h2κ : 1 ≤ 2 * κ := by linarith
  have hp1 : (2 * κ) ^ (4 * ((P.part ↑g (gC g)).M - 1)) ≤ (2 * κ) ^ (4 * (P.part ↑g (gC g)).M) :=
    pow_le_pow_right₀ h2κ (by omega)
  have hp2 := pow_le_sqrt_of_le_log hκ hε₁ (by have := abs_nonneg (gC g); linarith : (1 : ℝ) ≤ |gC g| + 2)
    hsmall hM
  have hc0 : 0 < (gC g) ^ 2 := by
    have : 1 ≤ |gC g| := by
      have := abs_sub_abs_le_abs_sub γ₀ (gC g)
      linarith
    have : 0 < |gC g| := by linarith
    rw [← sq_abs]; positivity
  have hF := sqrt_div_sq_le_F1 hγ₀2 hd
  have hnw : (1 : ℝ) ≤ nw P γ₀ := by exact_mod_cast nw_pos hP hγ₀N
  have hF0 : 0 ≤ F1 |γ₀| := F1_nonneg _
  calc (2 * κ) ^ (4 * ((P.part ↑g (gC g)).M - 1)) / (gC g) ^ 2
      ≤ Real.sqrt (|gC g| + 2) / (gC g) ^ 2 := by gcongr; exact hp1.trans hp2
    _ ≤ F1 |γ₀| := hF
    _ ≤ (nw P γ₀ : ℝ) * F1 |γ₀| := le_mul_of_one_le_left hF0 hnw
    _ ≤ ∑ γ ∈ g, (nw P γ : ℝ) * F1 |γ| :=
        Finset.single_le_sum (f := fun γ => (nw P γ : ℝ) * F1 |γ|)
          (fun γ _ => mul_nonneg (Nat.cast_nonneg _) (F1_nonneg _)) hγ₀

/-- All groups: Σ_g (2κ)^{4(M_g-1)}/c_g² ≤ Σ_{nodes} nw·F₁. -/
theorem groups_weight_le {P : ExpSum} {T : ℝ} (hP : ZeroConfig P T) {Th : Set ℝ} {s₀ κ ε₁ : ℝ}
    (hκ : 1 ≤ κ) (hε₁ : 0 ≤ ε₁) (hsmall : 4 * ε₁ * Real.log (2 * κ) ≤ 1 / 2)
    (hM : ∀ g ∈ groups P Th s₀, ((P.part ↑g (gC g)).M : ℝ) ≤ ε₁ * Real.log (|gC g| + 2))
    (hW : ∀ g ∈ groups P Th s₀, gW g ≤ 1) (hTh : ∀ γ ∈ nodeSet P Th, 2 ≤ |γ|) :
    ∑ g ∈ groups P Th s₀, (2 * κ) ^ (4 * ((P.part ↑g (gC g)).M - 1)) / (gC g) ^ 2 ≤
      ∑ γ ∈ nodeSet P Th, (nw P γ : ℝ) * F1 |γ| := by
  rw [← sum_groups_sum P Th s₀]
  exact Finset.sum_le_sum fun g hg => group_weight_le hP hg hκ hε₁ hsmall (hM g hg) (hW g hg) hTh

/-- The number of nodes is at most (Σ nw F₁)/F₁(T). -/
theorem card_nodes_le {P : ExpSum} {T : ℝ} (hP : ZeroConfig P T) (hT : 1 ≤ T) (Th : Set ℝ) :
    ((nodeSet P Th).card : ℝ) ≤ (∑ γ ∈ nodeSet P Th, (nw P γ : ℝ) * F1 |γ|) / F1 T := by
  have hF1T : 0 < F1 T := by rw [F1_of_ge_one hT]; positivity
  rw [le_div_iff₀ hF1T, Finset.card_eq_sum_ones, Nat.cast_sum, Finset.sum_mul]
  refine Finset.sum_le_sum fun γ hγ => ?_
  have hγT : |γ| ≤ T := node_height hP hγ
  have hnw : (1 : ℝ) ≤ nw P γ := by exact_mod_cast nw_pos hP hγ
  have hmono : F1 T ≤ F1 |γ| := F1_antitone (mem_Ici.mpr (abs_nonneg γ)) (mem_Ici.mpr (by linarith)) hγT
  push_cast
  calc 1 * F1 T ≤ 1 * F1 |γ| := by linarith
    _ ≤ (nw P γ : ℝ) * F1 |γ| := mul_le_mul_of_nonneg_right hnw (F1_nonneg _)

end Riemann.Deep.LeanNBRate

