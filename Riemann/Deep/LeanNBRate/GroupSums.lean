import Riemann.Deep.LeanNBRate.ZeroConfig

/-!
# Sums over the zeros: weights F₁, tails, and the group sums
-/

open MeasureTheory Filter Topology Set Complex ComplexConjugate

noncomputable section

namespace Riemann.Deep.LeanNBRate

/-- F₁(t) = 4√2 (max t 1 + 1)^{1/2} / (max t 1)². -/
def F1 (t : ℝ) : ℝ := 4 * Real.sqrt 2 * Real.sqrt (max t 1 + 1) / (max t 1) ^ 2

theorem F1_nonneg (t : ℝ) : 0 ≤ F1 t := by unfold F1; positivity

theorem F1_antitone : AntitoneOn F1 (Ici 0) := by
  intro x _ y _ hxy
  unfold F1
  set a := max x 1
  set b := max y 1
  have ha : 1 ≤ a := le_max_right _ _
  have hab : a ≤ b := max_le_max hxy le_rfl
  -- √(b+1)/b² ≤ √(a+1)/a²  ⟸  (b+1) a⁴ ≤ (a+1) b⁴
  have hpos : 0 < a := by linarith
  have key : Real.sqrt (b + 1) * a ^ 2 ≤ Real.sqrt (a + 1) * b ^ 2 := by
    have h1 : Real.sqrt (b + 1) ≤ Real.sqrt (a + 1) * (b / a) := by
      rw [show Real.sqrt (a + 1) * (b / a) = Real.sqrt ((a + 1) * (b / a) ^ 2) by
        rw [Real.sqrt_mul (by linarith), Real.sqrt_sq (by positivity)]]
      apply Real.sqrt_le_sqrt
      have : b / a ≥ 1 := by rw [ge_iff_le, le_div_iff₀ hpos]; linarith
      have hba : (b / a) ^ 2 ≥ b / a := by nlinarith
      have : (a + 1) * (b / a) = b + b / a := by field_simp
      nlinarith
    have h2 : Real.sqrt (a + 1) * (b / a) * a ^ 2 = Real.sqrt (a + 1) * (b * a) := by
      field_simp
    have h3 : Real.sqrt (a + 1) * (b * a) ≤ Real.sqrt (a + 1) * b ^ 2 := by
      have : b * a ≤ b ^ 2 := by nlinarith
      exact mul_le_mul_of_nonneg_left this (Real.sqrt_nonneg _)
    calc Real.sqrt (b + 1) * a ^ 2 ≤ Real.sqrt (a + 1) * (b / a) * a ^ 2 := by gcongr
      _ = Real.sqrt (a + 1) * (b * a) := h2
      _ ≤ _ := h3
  have hb : 0 < b := by linarith
  rw [div_le_div_iff₀ (by positivity) (by positivity)]
  have : 0 ≤ 4 * Real.sqrt 2 := by positivity
  nlinarith

theorem summable_F1_log : Summable fun k : ℕ => F1 k * Real.log (k + 1) := by
  have hs : Summable fun k : ℕ => (k : ℝ) ^ (-(5 / 4 : ℝ)) :=
    Real.summable_nat_rpow.mpr (by norm_num)
  refine Summable.of_norm_bounded_eventually_nat
    (g := fun k : ℕ => 16 * Real.sqrt 2 * 2 * (k : ℝ) ^ (-(5 / 4 : ℝ))) (hs.mul_left _) ?_
  filter_upwards [Filter.eventually_ge_atTop 1] with k hk
  have hk1 : (1 : ℝ) ≤ k := by exact_mod_cast hk
  have hk0 : (0 : ℝ) < k := by linarith
  have hmax : max (k : ℝ) 1 = k := max_eq_left hk1
  unfold F1
  rw [hmax]
  have hlog : Real.log (k + 1) ≤ 4 * (k + 1 : ℝ) ^ (1 / 4 : ℝ) := by
    have := Real.log_le_rpow_div (by positivity : (0 : ℝ) ≤ k + 1) (by norm_num : (0 : ℝ) < 1 / 4)
    rw [div_eq_mul_inv] at this
    calc Real.log (k + 1) ≤ (k + 1 : ℝ) ^ (1 / 4 : ℝ) * (1 / 4)⁻¹ := this
      _ = 4 * (k + 1 : ℝ) ^ (1 / 4 : ℝ) := by norm_num; ring
  have hlog0 : 0 ≤ Real.log (k + 1) := Real.log_nonneg (by linarith)
  rw [Real.norm_eq_abs, abs_of_nonneg (mul_nonneg (by positivity) hlog0)]
  -- (k+1) ≤ 2k
  have hk2 : (k + 1 : ℝ) ≤ 2 * k := by linarith
  have hsq : Real.sqrt (k + 1) = (k + 1 : ℝ) ^ (1 / 2 : ℝ) := Real.sqrt_eq_rpow _
  rw [hsq]
  have e1 : (k + 1 : ℝ) ^ (1 / 2 : ℝ) * (k + 1 : ℝ) ^ (1 / 4 : ℝ) = (k + 1 : ℝ) ^ (3 / 4 : ℝ) := by
    rw [← Real.rpow_add (by positivity)]; norm_num
  have e2 : (k + 1 : ℝ) ^ (3 / 4 : ℝ) ≤ (2 * k : ℝ) ^ (3 / 4 : ℝ) :=
    Real.rpow_le_rpow (by positivity) hk2 (by norm_num)
  have e3 : (2 * k : ℝ) ^ (3 / 4 : ℝ) = (2 : ℝ) ^ (3 / 4 : ℝ) * (k : ℝ) ^ (3 / 4 : ℝ) :=
    Real.mul_rpow (by norm_num) hk0.le
  have e4 : (2 : ℝ) ^ (3 / 4 : ℝ) ≤ 2 := by
    calc (2 : ℝ) ^ (3 / 4 : ℝ) ≤ (2 : ℝ) ^ (1 : ℝ) := Real.rpow_le_rpow_of_exponent_le (by norm_num) (by norm_num)
      _ = 2 := Real.rpow_one 2
  have e5 : (k : ℝ) ^ (3 / 4 : ℝ) / (k : ℝ) ^ 2 = (k : ℝ) ^ (-(5 / 4 : ℝ)) := by
    rw [show ((k : ℝ) ^ 2) = (k : ℝ) ^ (2 : ℝ) by norm_cast, ← Real.rpow_sub hk0]; norm_num
  calc 4 * Real.sqrt 2 * (k + 1 : ℝ) ^ (1 / 2 : ℝ) / (k : ℝ) ^ 2 * Real.log (k + 1)
      ≤ 4 * Real.sqrt 2 * (k + 1 : ℝ) ^ (1 / 2 : ℝ) / (k : ℝ) ^ 2 * (4 * (k + 1 : ℝ) ^ (1 / 4 : ℝ)) := by
        gcongr
    _ = 16 * Real.sqrt 2 * ((k + 1 : ℝ) ^ (1 / 2 : ℝ) * (k + 1 : ℝ) ^ (1 / 4 : ℝ)) / (k : ℝ) ^ 2 := by ring
    _ = 16 * Real.sqrt 2 * (k + 1 : ℝ) ^ (3 / 4 : ℝ) / (k : ℝ) ^ 2 := by rw [e1]
    _ ≤ 16 * Real.sqrt 2 * ((2 : ℝ) ^ (3 / 4 : ℝ) * (k : ℝ) ^ (3 / 4 : ℝ)) / (k : ℝ) ^ 2 := by
        rw [← e3]; gcongr
    _ ≤ 16 * Real.sqrt 2 * (2 * (k : ℝ) ^ (3 / 4 : ℝ)) / (k : ℝ) ^ 2 := by gcongr
    _ = 16 * Real.sqrt 2 * 2 * ((k : ℝ) ^ (3 / 4 : ℝ) / (k : ℝ) ^ 2) := by ring
    _ = _ := by rw [e5]

/-- The zero weights mult ρ · F₁(|γ|) are summable. -/
theorem summable_wF1 (hc : Titchmarsh_Thm9_2) :
    Summable fun ρ : NTZ => (mult (ρ : ℂ) : ℝ) * F1 |(ρ : ℂ).im| :=
  summable_NTZ_of_shells F1 F1_nonneg F1_antitone (summable_shell hc F1 F1_nonneg summable_F1_log)

/-- A finite sum over the nodes of a zero configuration is at most the corresponding tsum. -/
theorem node_sum_le_tsum {P : ExpSum} {T : ℝ} (hP : ZeroConfig P T) (T' : Set ℝ) (F : ℝ → ℝ)
    (hF0 : ∀ t, 0 ≤ F t) (hs : Summable fun ρ : NTZ => (mult (ρ : ℂ) : ℝ) * F |(ρ : ℂ).im|)
    (A : Finset ℝ) (hA : A ⊆ nodeSet P T') :
    ∑ γ ∈ A, (nw P γ : ℝ) * F |γ| ≤ ∑' ρ : NTZ, (mult (ρ : ℂ) : ℝ) * F |(ρ : ℂ).im| := by
  classical
  set e : A → NTZ := fun γ => ⟨(1 / 2 : ℂ) + (γ.1 : ℂ) * I, node_mem_NTZ hP (hA γ.2)⟩
  have he : Function.Injective e := by
    intro x y hxy
    have := congrArg (fun ρ : NTZ => (ρ : ℂ).im) hxy
    simp [e] at this
    exact this
  have h1 : ∑ γ ∈ A, (nw P γ : ℝ) * F |γ| ≤ ∑ γ ∈ A.attach, (mult ((e γ : NTZ) : ℂ) : ℝ) * F |((e γ : NTZ) : ℂ).im| := by
    rw [← Finset.sum_attach A]
    refine Finset.sum_le_sum fun γ _ => ?_
    have h2 : ((e γ : NTZ) : ℂ).im = γ.1 := by simp [e]
    rw [h2]
    exact mul_le_mul_of_nonneg_right (by exact_mod_cast nw_le_mult hP (hA γ.2)) (hF0 _)
  refine h1.trans ?_
  rw [← Finset.sum_image (f := fun ρ : NTZ => (mult (ρ : ℂ) : ℝ) * F |(ρ : ℂ).im|) (g := e)
    (s := A.attach) (fun x _ y _ hxy => he hxy)]
  exact hs.sum_le_tsum _ fun ρ _ => mul_nonneg (Nat.cast_nonneg _) (hF0 _)

/-- Tails: node sums above a large height are small. -/
theorem node_sum_tail (F : ℝ → ℝ) (hF0 : ∀ t, 0 ≤ F t)
    (hs : Summable fun ρ : NTZ => (mult (ρ : ℂ) : ℝ) * F |(ρ : ℂ).im|) (ε : ℝ) (hε : 0 < ε) :
    ∃ H₀ : ℝ, ∀ H ≥ H₀, ∀ (P : ExpSum) (T : ℝ), ZeroConfig P T → ∀ (T' : Set ℝ) (A : Finset ℝ),
      A ⊆ nodeSet P T' → (∀ γ ∈ A, H < |γ|) → ∑ γ ∈ A, (nw P γ : ℝ) * F |γ| ≤ ε := by
  classical
  obtain ⟨s, hs'⟩ := (summable_iff_vanishing_norm.mp hs) ε hε
  refine ⟨1 + ∑ ρ ∈ s, |(ρ : ℂ).im|, fun H hH P T hP T' A hA hAH => ?_⟩
  have hHs : ∀ ρ ∈ s, |(ρ : ℂ).im| < H := by
    intro ρ hρ
    have h1 : |(ρ : ℂ).im| ≤ ∑ ρ ∈ s, |(ρ : ℂ).im| :=
      Finset.single_le_sum (f := fun ρ : NTZ => |(ρ : ℂ).im|) (fun _ _ => abs_nonneg _) hρ
    linarith
  set e : A → NTZ := fun γ => ⟨(1 / 2 : ℂ) + (γ.1 : ℂ) * I, node_mem_NTZ hP (hA γ.2)⟩
  have he : Function.Injective e := by
    intro x y hxy
    have := congrArg (fun ρ : NTZ => (ρ : ℂ).im) hxy
    simp [e] at this
    exact this
  set t := A.attach.image e
  have hdisj : Disjoint t s := by
    rw [Finset.disjoint_left]
    intro ρ hρt hρs
    obtain ⟨γ, _, rfl⟩ := Finset.mem_image.mp hρt
    have h1 := hHs _ hρs
    have h2 : ((e γ : NTZ) : ℂ).im = γ.1 := by simp [e]
    rw [h2] at h1
    linarith [hAH γ.1 γ.2]
  have hv := hs' t hdisj
  rw [Real.norm_eq_abs] at hv
  have hsum0 : 0 ≤ ∑ ρ ∈ t, (mult (ρ : ℂ) : ℝ) * F |(ρ : ℂ).im| :=
    Finset.sum_nonneg fun ρ _ => mul_nonneg (Nat.cast_nonneg _) (hF0 _)
  rw [abs_of_nonneg hsum0] at hv
  have h1 : ∑ γ ∈ A, (nw P γ : ℝ) * F |γ| ≤ ∑ ρ ∈ t, (mult (ρ : ℂ) : ℝ) * F |(ρ : ℂ).im| := by
    rw [← Finset.sum_attach A, Finset.sum_image (f := fun ρ : NTZ => (mult (ρ : ℂ) : ℝ) * F |(ρ : ℂ).im|)
      (g := e) (s := A.attach) (fun x _ y _ hxy => he hxy)]
    refine Finset.sum_le_sum fun γ _ => ?_
    have h2 : ((e γ : NTZ) : ℂ).im = γ.1 := by simp [e]
    rw [h2]
    exact mul_le_mul_of_nonneg_right (by exact_mod_cast nw_le_mult hP (hA γ.2)) (hF0 _)
  linarith

end Riemann.Deep.LeanNBRate

