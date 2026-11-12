import Riemann.Deep.DPCUncond.Lean.DLam
import Riemann.Deep.DPCUncond.Lean.Hairpin

/-!
# dpcuncond: Steps 5 and 6 from finer named inputs (density theorem, depth classes, Montgomery–Vaughan)

Step 6 (the hairpin bound) is derived, not assumed: Density (Jutila/Conrey, uniform for σ ≥ 1/2, Lz l.221–230) and
Step 2 (card_off_le) bound the depth classes S_a by min(3ε, C_d e^{−κa}) T L (6a); Step6Classes (6a's decomposition
of 𝓗_T into classes, 6b–6e's class bounds A₁, A₂, A₃ and Hölder (3, 3, 3)) is the named analytic input;
holder_bookkeeping (6e) and hairpin_sum (6f) are proved. Step 5a is derived from Montgomery–Vaughan at length
T^{2θ} (MV5a) using θ ≤ 1/2; Step 5c stays named. theoremU_fine, D_pos_two_thirds_fine: the main results from these.
-/

namespace DPCLean

open Finset PPLean

variable {ex lg : ℝ → ℝ}

/-- Off-line zeros with T ≤ γ ≤ 2T and depth (β − 1/2) log T ≥ a. -/
noncomputable def deepSet (Z : Zeros) (β : ℕ → ℝ) (lg : ℝ → ℝ) (T : ℝ) (a : ℕ) : Finset ℕ :=
  (offSet Z β T).filter (fun n => (a : ℝ) ≤ (β n - 1 / 2) * lg T)

/-- The depth class S_a of 6a: off-line zeros with T ≤ γ ≤ 2T and (β − 1/2) log T ∈ [a, a + 1). -/
noncomputable def classSet (Z : Zeros) (β : ℕ → ℝ) (lg : ℝ → ℝ) (T : ℝ) (a : ℕ) : Finset ℕ :=
  (offSet Z β T).filter (fun n => (a : ℝ) ≤ (β n - 1 / 2) * lg T ∧ (β n - 1 / 2) * lg T < a + 1)

/-- The zero-density theorem, uniform for σ ≥ 1/2 (Jutila 1983, slope κ = 1 − δ for every δ > 0; Conrey 8/7 − δ;
as quoted Lz 2611.14486 l.221–230): N(σ, 2T) ≪_δ T^{1−κ(σ−1/2)} log T at σ = 1/2 + a/log T gives
#{T ≤ γ ≤ 2T, β ≥ 1/2 + a/log T} ≤ C_d T log T e^{−κa}. -/
def Density (Z : Zeros) (β : ℕ → ℝ) (ex lg : ℝ → ℝ) : Prop :=
  ∀ δ : ℝ, 0 < δ → δ < 1 → ∃ Cd T₀ : ℝ, 0 < Cd ∧ ∀ T : ℝ, T₀ ≤ T → ∀ a : ℕ,
    ((deepSet Z β lg T a).card : ℝ) ≤ Cd * T * lg T * ex (-((1 - δ) * a))

/-- Steps 6a–6e at the level of the depth classes S_a: for all w, r > 1 there are C, K and class totals
Ha T a = ∫|R_T||R*|(1 + |V_L| + 2Ñ)^D dμ_a and A₁ = ∫|R_T|³dμ_a, A₂ = ∫|R*|³dμ_a, A₃ = ∫(1 + |V_L| + 2Ñ)^{3D}dμ_a
such that, for large T, (6a) |𝓗_T| ≤ Σ_{a≤M(T)} Ha T a + o(T) (each merged cut belongs to the class of its deepest
zero; zeros outside [T, 2T] give o(1), Step 3), (6e) Ha³ ≤ A₁A₂A₃ (Hölder) and, whenever |S_a| ≤ s T log T,
A₁ ≤ K T (a+2)^C s^{(1−1/w)/4}, A₂ ≤ K T (a+2)^C e^{3θ(a+1)} s^{(1−1/w)/4}, A₃ ≤ K T (a+2)^C s^{(1−1/w)(1−1/r)}
(6b–6d: Poisson representations, sparse kernels, Fujii's moments, Lz Lemma DoubleSumZeros l.251–283). -/
def Step6Classes (Z : Zeros) (β : ℕ → ℝ) (ex lg : ℝ → ℝ) (θ : ℝ) (D : Data) : Prop :=
  ∀ w r : ℝ, 1 < w → 1 < r → ∃ (C : ℕ) (K : ℝ), 0 ≤ K ∧
    ∃ (Ha A₁ A₂ A₃ : ℝ → ℕ → ℝ) (Mt : ℝ → ℕ), ∀ η : ℝ, 0 < η → ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T →
      |D.Hp T| ≤ ∑ a ∈ range (Mt T + 1), Ha T a + η * T ∧
      ∀ a : ℕ, 0 ≤ A₂ T a ∧ 0 ≤ A₃ T a ∧ Ha T a ^ 3 ≤ A₁ T a * A₂ T a * A₃ T a ∧
        ∀ s : ℝ, 0 < s → ((classSet Z β lg T a).card : ℝ) ≤ s * T * lg T →
          A₁ T a ≤ K * T * ((a : ℝ) + 2) ^ C * ex (lg s * ((1 - 1 / w) / 4)) ∧
          A₂ T a ≤ K * T * ((a : ℝ) + 2) ^ C * (ex θ ^ (a + 1)) ^ 3 * ex (lg s * ((1 - 1 / w) / 4)) ∧
          A₃ T a ≤ K * T * ((a : ℝ) + 2) ^ C * ex (lg s * ((1 - 1 / w) * (1 - 1 / r)))

/-- **Step 6 from the density theorem and the class bounds** (6a counts, 6e bookkeeping and 6f sum proved). -/
theorem step6_of_classes (hEL : ExpLog ex lg) {Z : Zeros} {β : ℕ → ℝ} {θ : ℝ} (hθ : 0 ≤ θ) {D : Data}
    (hden : Density Z β ex lg) (h6c : Step6Classes Z β ex lg θ D) : Step6 Z β ex lg θ D := by
  intro δ w r hδ0 hδ1 hw hr hκE
  obtain ⟨C, K, hK, Ha, A₁, A₂, A₃, Mt, hcl⟩ := h6c w r hw hr
  obtain ⟨Cd, Td, hCd, hdT⟩ := hden δ hδ0 hδ1
  have hE0 : 0 ≤ Eexp w r := by nlinarith
  refine ⟨C + 1, K * max 1 (ex (lg Cd * Eexp w r)),
    mul_nonneg hK (le_trans zero_le_one (le_max_left _ _)), fun c _ _ ε hε _ η hη => ?_⟩
  obtain ⟨T₀, hT₀⟩ := hcl η hη
  refine ⟨max (max T₀ Td) (ex 1), fun T hT hoff => ?_⟩
  have hT1 : T₀ ≤ T := le_trans (le_trans (le_max_left _ _) (le_max_left _ _)) hT
  have hT2 : Td ≤ T := le_trans (le_trans (le_max_right _ _) (le_max_left _ _)) hT
  have hT3 : ex 1 ≤ T := le_trans (le_max_right _ _) hT
  have hT0 : 0 < T := lt_of_lt_of_le (hEL.pos 1) hT3
  have hlg : 1 ≤ lg T := hEL.le_lg hT3
  have hTL : 0 ≤ T * lg T := mul_nonneg hT0.le (by linarith)
  obtain ⟨hdec, hcls⟩ := hT₀ T hT1
  refine ⟨Mt T, ?_⟩
  have hs0 : ∀ a : ℕ, 0 < min (3 * ε) (Cd * ex (-((1 - δ) * a))) :=
    fun a => lt_min (by positivity) (mul_pos hCd (hEL.pos _))
  have hbound : ∀ a : ℕ, Ha T a ≤ K * T * ((a : ℝ) + 2) ^ C * ex θ ^ (a + 1) *
      ex (lg (min (3 * ε) (Cd * ex (-((1 - δ) * a)))) * Eexp w r) := by
    intro a
    obtain ⟨hA2, hA3, hH3, hb⟩ := hcls a
    have hc1 : ((classSet Z β lg T a).card : ℝ) ≤ 3 * ε * T * lg T := by
      have h2 : ((classSet Z β lg T a).card : ℝ) ≤ ((offSet Z β T).card : ℝ) := by
        exact_mod_cast card_le_card (filter_subset _ _)
      have h3 : 0 ≤ ε * (T * lg T) := mul_nonneg hε.le hTL
      linarith
    have hc2 : ((classSet Z β lg T a).card : ℝ) ≤ Cd * ex (-((1 - δ) * a)) * T * lg T := by
      have h1 : classSet Z β lg T a ⊆ deepSet Z β lg T a := by
        intro n hn
        simp only [classSet, deepSet, mem_filter] at hn ⊢
        exact ⟨hn.1, hn.2.1⟩
      have h2 : ((classSet Z β lg T a).card : ℝ) ≤ ((deepSet Z β lg T a).card : ℝ) := by
        exact_mod_cast card_le_card h1
      have h3 := hdT T hT2 a
      have e : Cd * T * lg T * ex (-((1 - δ) * a)) = Cd * ex (-((1 - δ) * a)) * T * lg T := by ring
      linarith
    have hcnt : ((classSet Z β lg T a).card : ℝ) ≤
        min (3 * ε) (Cd * ex (-((1 - δ) * a))) * T * lg T := by
      rcases le_total (3 * ε) (Cd * ex (-((1 - δ) * a))) with h | h
      · rw [min_eq_left h]
        exact hc1
      · rw [min_eq_right h]
        exact hc2
    obtain ⟨hb1, hb2, hb3⟩ := hb _ (hs0 a) hcnt
    exact holder_bookkeeping (s := fun e => ex (lg (min (3 * ε) (Cd * ex (-((1 - δ) * a)))) * e))
      (fun x y => by simp only [mul_add, hEL.add]) (fun x => (hEL.pos _).le) hK hT0.le
      (by positivity) (by exact pow_nonneg (hEL.pos _).le _) hA2 hA3 hb1 hb2 hb3 hH3
  have hsum := hairpin_sum hEL (C := C) (θ := θ) (κ := 1 - δ) (E := Eexp w r) (c := c) (ε := ε)
    (Mt T) (Ha T) (fun a => min (3 * ε) (Cd * ex (-((1 - δ) * a)))) hθ hE0 hK hT0.le hCd hs0
    (fun a => min_le_left _ _) (fun a => min_le_right _ _) hbound
  linarith

/-- Step 5a, Montgomery–Vaughan's mean value theorem for R_T² = Σ b_n n^{−it} (length T^{2θ}) with Σ|b_n|² ≤ B
(Mertens, the T^{η/2}-rough support): ∫|R_T|⁴ ≤ (2T + K T^{2θ}) B, T^{2θ} = e^{2θ log T}. -/
def MV5a (ex lg : ℝ → ℝ) (θ : ℝ) (D : Data) : Prop :=
  ∃ Kmv B T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → 0 ≤ D.R4 T ∧ D.R4 T ≤ (2 * T + Kmv * ex (lg T * (2 * θ))) * B

/-- Step 5c: ∫|P(V)|⁴ ≪ T (M. Riesz DPC l.641, Fujii Lz l.242–249, density; 5b–5c). -/
def Step5c (D : Data) : Prop :=
  ∃ K₂ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → 0 ≤ D.P4 T ∧ D.P4 T ≤ K₂ * T

/-- **Step 5 from Montgomery–Vaughan**: T^{2θ} ≤ T for θ ≤ 1/2, so ∫|R_T|⁴ ≪ T. -/
theorem step5_of (hEL : ExpLog ex lg) {θ : ℝ} (hθ : θ ≤ 1 / 2) {D : Data} (h5a : MV5a ex lg θ D)
    (h5c : Step5c D) : Step5 D := by
  obtain ⟨Kmv, B, T₁, hT₁⟩ := h5a
  obtain ⟨K₂, T₂, hT₂⟩ := h5c
  refine ⟨(2 + |Kmv|) * |B|, K₂, max (max T₁ T₂) (ex 0), fun T hT => ?_⟩
  have hT1 : T₁ ≤ T := le_trans (le_trans (le_max_left _ _) (le_max_left _ _)) hT
  have hT2 : T₂ ≤ T := le_trans (le_trans (le_max_right _ _) (le_max_left _ _)) hT
  have hT3 : ex 0 ≤ T := le_trans (le_max_right _ _) hT
  have hT0 : 0 < T := lt_of_lt_of_le (hEL.pos 0) hT3
  have hlg : 0 ≤ lg T := hEL.le_lg hT3
  have hX : ex (lg T * (2 * θ)) ≤ T := by
    have h1 : lg T * (2 * θ) ≤ lg T := by
      nlinarith [mul_nonneg hlg (show (0 : ℝ) ≤ 1 - 2 * θ by linarith)]
    calc ex (lg T * (2 * θ)) ≤ ex (lg T) := hEL.mono.monotone h1
      _ = T := hEL.ex_lg T hT0
  have hX0 : 0 ≤ ex (lg T * (2 * θ)) := (hEL.pos _).le
  obtain ⟨hR0, hR1⟩ := hT₁ T hT1
  obtain ⟨hP0, hP1⟩ := hT₂ T hT2
  refine ⟨hR0, ?_, hP0, hP1⟩
  have habs : |2 * T + Kmv * ex (lg T * (2 * θ))| ≤ (2 + |Kmv|) * T := by
    calc |2 * T + Kmv * ex (lg T * (2 * θ))| ≤ |2 * T| + |Kmv * ex (lg T * (2 * θ))| := abs_add_le _ _
      _ = 2 * T + |Kmv| * ex (lg T * (2 * θ)) := by
          rw [abs_of_nonneg (show (0 : ℝ) ≤ 2 * T by linarith), abs_mul, abs_of_nonneg hX0]
      _ ≤ 2 * T + |Kmv| * T := by
          have h1 := mul_le_mul_of_nonneg_left hX (abs_nonneg Kmv)
          linarith
      _ = (2 + |Kmv|) * T := by ring
  calc D.R4 T ≤ (2 * T + Kmv * ex (lg T * (2 * θ))) * B := hR1
    _ ≤ |(2 * T + Kmv * ex (lg T * (2 * θ))) * B| := le_abs_self _
    _ = |2 * T + Kmv * ex (lg T * (2 * θ))| * |B| := abs_mul _ _
    _ ≤ (2 + |Kmv|) * T * |B| := mul_le_mul_of_nonneg_right habs (abs_nonneg B)
    _ = (2 + |Kmv|) * |B| * T := by ring

/-- The finer named inputs for one Step-0 construction. -/
structure InputsFine (Z : Zeros) (β : ℕ → ℝ) (ex lg : ℝ → ℝ) (φ θ : ℝ) (D : Data) : Prop where
  step3 : Step3 D
  main : MainTerm D
  step6c : Step6Classes Z β ex lg θ D
  step1 : Step1Pos D
  holder : Holder7 D
  mv5a : MV5a ex lg θ D
  step5c : Step5c D
  cover : Covering Z φ D

/-- DPCInputs with Step 6 replaced by the class-level inputs and Step 5a by Montgomery–Vaughan. -/
def DPCInputsFine (Z : Zeros) (β : ℕ → ℝ) (ex lg : ℝ → ℝ) : Prop :=
  ∀ φ θ : ℝ, 0 < φ → φ < 1 → 0 < θ → θ < 1 / 2 → 0 < φ * θ → φ * θ < 1 → 1 < φ * (1 + θ) →
    φ * (1 + θ) < 2 → ∃ D : Data, InputsFine Z β ex lg φ θ D

theorem inputs_of_fine (hEL : ExpLog ex lg) {Z : Zeros} {β : ℕ → ℝ} (hden : Density Z β ex lg)
    (hin : DPCInputsFine Z β ex lg) : DPCInputs Z β ex lg := by
  intro φ θ h1 h2 h3 h4 h5 h6 h7 h8
  obtain ⟨D, inp⟩ := hin φ θ h1 h2 h3 h4 h5 h6 h7 h8
  exact ⟨D, ⟨inp.step3, inp.main, step6_of_classes hEL h3.le hden inp.step6c, inp.step1, inp.holder,
    step5_of hEL h4.le inp.mv5a inp.step5c, inp.cover⟩⟩

/-- **Theorem U** from the finer inputs (density theorem, class bounds, Montgomery–Vaughan). -/
theorem theoremU_fine {Z : Zeros} (hZ : Z.Ordered) {β : ℕ → ℝ} (hβ : Reflect Z β) (hEL : ExpLog ex lg)
    {tp : ℝ} (htp : 0 < tp) (hsc : ∀ T : ℝ, 0 < T → Z.sc T * tp = lg T) (hden : Density Z β ex lg)
    (hin : DPCInputsFine Z β ex lg) : TheoremU Z :=
  theoremU hZ hβ hEL htp hsc (inputs_of_fine hEL hden hin)

/-- **D(λ) > 0 for every λ > 2/3** from the finer inputs. -/
theorem D_pos_two_thirds_fine {Z : Zeros} (hZ : Z.Ordered) {β : ℕ → ℝ} (hβ : Reflect Z β)
    (hEL : ExpLog ex lg) {tp : ℝ} (htp : 0 < tp) (hsc : ∀ T : ℝ, 0 < T → Z.sc T * tp = lg T)
    (hR : RvMUp Z) (hden : Density Z β ex lg) (hin : DPCInputsFine Z β ex lg) :
    ∀ lam : ℝ, 2 / 3 < lam → 0 < Dlam Z lam :=
  D_pos_two_thirds hZ hβ hEL htp hsc hR (inputs_of_fine hEL hden hin)

/-- Theorem U along any sequence T_k → ∞ (T → ∞ indexed by ℕ). -/
theorem theoremU_seq {Z : Zeros} (hU : TheoremU Z) (t : ℕ → ℝ)
    (ht : ∀ M : ℝ, ∃ k₀ : ℕ, ∀ k : ℕ, k₀ ≤ k → M ≤ t k) {φ : ℝ} (h1 : 2 / 3 < φ) (h2 : φ < 1) :
    ∃ c : ℝ, 0 < c ∧ ∃ k₀ : ℕ, ∀ k : ℕ, k₀ ≤ k → c * t k * Z.sc (t k) ≤ ((gapSet Z φ (t k)).card : ℝ) := by
  obtain ⟨c, hc, T₀, hT₀⟩ := hU φ h1 h2
  obtain ⟨k₀, hk₀⟩ := ht T₀
  exact ⟨c, hc, k₀, fun k hk => hT₀ (t k) (hk₀ k hk)⟩

end DPCLean
