import Riemann.Deep.DPCUncond.Lean.Gaps
import Riemann.Deep.DPCUncond.Lean.Params
import Riemann.Deep.DPCUncond.Lean.Phi
import Riemann.Deep.DPCUncond.Lean.Step7

/-!
# dpcuncond: Theorem U from the named analytic inputs (certificate-interface pattern)

TheoremU Z: for every φ ∈ (2/3, 1) there is c > 0 with 𝒩(T) ≥ c T sc T (sc T = log T/(2π), so c T sc T =
(c/2π) T log T) for all large T. theoremU: TheoremU from the zero model (Zeros.Ordered, Reflect), the abstract
exp/log (ExpLog, sc T · 2π = log T) and DPCInputs: for every admissible (φ, θ) the analytic quantities Data of
DPC's Step-0 construction satisfy the named hypotheses Step3, MainTerm, Step6, Step1Pos, Holder7, Step5, Covering.
Proved here, not assumed: the choice of θ (Step 0), of δ, w, r (6e), of c = E/(2θ) (Repair 6), Φ(ε) → 0 (6f),
#Z_off ≤ 2𝒩 (Step 2) and the contradiction (Step 7).
-/

namespace DPCLean

open Finset PPLean

/-- **Theorem U** (PROOF.md), in pplean's zero model, sc T = log T/(2π): 𝒩(T) ≥ c(φ) T log T, unfolded. -/
def TheoremU (Z : Zeros) : Prop :=
  ∀ φ : ℝ, 2 / 3 < φ → φ < 1 → ∃ c : ℝ, 0 < c ∧ ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T →
    c * T * Z.sc T ≤ ((gapSet Z φ T).card : ℝ)

/-- The analytic quantities at fixed Step-0 data (φ, θ, q, b₀, P, f, η; PROOF Step 0), as functions of T. -/
structure Data where
  /-- Re ∫_T^{2T} G_T(t)|R_T(t)|² P(V_φ(t)) dt (DPC l.690, 769; Step 3's ∫_ℝ differs by O(e^{−cL²}), l.698) -/
  I : ℝ → ℝ
  /-- Re M_T, the Im t = −1 term after the shift (Step 3; DPC Prop 1, l.341–388) -/
  M : ℝ → ℝ
  /-- Re 𝓗_T, the hairpin total, a sum over Z_off (Steps 3, 6) -/
  Hp : ℝ → ℝ
  /-- Σ_n |c_n|² -/
  S : ℝ → ℝ
  /-- ∫_T^{2T} |R_T|⁴ (Step 5a) -/
  R4 : ℝ → ℝ
  /-- ∫_T^{2T} |P(V_φ)|⁴ (Step 5c) -/
  P4 : ℝ → ℝ
  /-- ∫_{E_T} G_T |R_T|² |P(V_φ)| (Step 7) -/
  Bad : ℝ → ℝ
  /-- meas E_T, E_T = {t ∈ [T, 2T] : the window (t − h/2, t + h/2) holds ≥ 2 zeros} (Step 7) -/
  mE : ℝ → ℝ

/-- Step 3, the exact moment identity (H) (DPC l.335–388 after the shift): ∫ = M_T + 𝓗_T + o(T), real parts. -/
def Step3 (D : Data) : Prop :=
  ∀ η : ℝ, 0 < η → ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → |D.I T - (D.M T + D.Hp T)| ≤ η * T

/-- Main term: DPC Prop 1 (l.262–388) and Lemma 6 (l.686–741: the limiting quotient l.689–693; l.695: the
denominator divided by T tends to a positive constant; l.700: ∫G_T|R_T|² ~ (T/2)Σ|c_n|²) for Step 0's construction
(l.743–771): Re M_T ≤ −3δ₀(T/2)Σ|c_n|²(1 + o(1)), and Σ|c_n|² ≥ S₀ > 0 (Step 7's δ₁ > 0). -/
def MainTerm (D : Data) : Prop :=
  ∃ δ₀ S₀ : ℝ, 0 < δ₀ ∧ 0 < S₀ ∧ ∀ η : ℝ, 0 < η → ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T →
    S₀ ≤ D.S T ∧ D.M T ≤ -3 * δ₀ * (T / 2) * D.S T * (1 - η)

/-- Step 6 (6a–6f; density N(σ, T) ≪ T^{1−κ(σ−1/2)} log T with Jutila's κ = 1 − δ, Lz l.223–225): for all
admissible δ, w, r (θ < (1 − δ)E) there are C, K, independent of ε, such that for every c ∈ (0, E/θ): if
#(Z_off ∩ [T, 2T]) ≤ 2εTL (Step 2's output), then |𝓗_T| ≤ K Φ(ε) T + o(T), L = log T. -/
def Step6 (Z : Zeros) (β : ℕ → ℝ) (ex lg : ℝ → ℝ) (θ : ℝ) (D : Data) : Prop :=
  ∀ δ w r : ℝ, 0 < δ → δ < 1 → 1 < w → 1 < r → θ < (1 - δ) * Eexp w r →
  ∃ (C : ℕ) (K : ℝ), 0 ≤ K ∧ ∀ c : ℝ, 0 < c → c * θ < Eexp w r →
  ∀ ε : ℝ, 0 < ε → ε < 1 → ∀ η : ℝ, 0 < η → ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T →
    ((offSet Z β T).card : ℝ) ≤ 2 * ε * T * lg T →
    ∃ m : ℕ, |D.Hp T| ≤ K * Phi ex lg C θ (1 - δ) (Eexp w r) c ε m * T + η * T

/-- Step 1 (Im V counts all zeros, so Im V ∈ {0, 1} outside E_T, where Re P(V) ≥ 0 by the two-line signs
l.760–765; DPC l.774): Re ∫_T^{2T} G_T|R_T|²P(V) ≥ −∫_{E_T} G_T|R_T|²|P(V)| − o(T). -/
def Step1Pos (D : Data) : Prop :=
  ∀ η : ℝ, 0 < η → ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → -D.Bad T - η * T ≤ D.I T

/-- Step 7, Hölder with exponents (2, 4, 4), G_T bounded (DPC l.776–780):
(∫_{E_T} G_T|R_T|²|P(V)|)⁴ ≤ KG (∫_T^{2T}|R_T|⁴)² ∫_T^{2T}|P(V)|⁴ meas E_T. -/
def Holder7 (D : Data) : Prop :=
  ∃ KG T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → 0 ≤ D.Bad T ∧ D.Bad T ^ 4 ≤ KG * D.R4 T ^ 2 * D.P4 T * D.mE T

/-- Step 5 (5a Montgomery–Vaughan at length T^{2θ} < T; 5b Fujii; 5c M. Riesz + density): ∫|R_T|⁴ ≪ T and
∫|P(V)|⁴ ≪ T, constants independent of ε. -/
def Step5 (D : Data) : Prop :=
  ∃ K₁ K₂ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → 0 ≤ D.R4 T ∧ D.R4 T ≤ K₁ * T ∧ 0 ≤ D.P4 T ∧ D.P4 T ≤ K₂ * T

/-- Step 7, covering (DPC l.774, 782): meas E_T ≤ h𝒩(T) + h with h = 2πφ/log T = φ/sc T (a consecutive pair with
gap < h covers centres of measure ≤ h; pairs not counted by 𝒩(T) cover only [T, T + h/2) ∪ (2T − h/2, 2T]). -/
def Covering (Z : Zeros) (φ : ℝ) (D : Data) : Prop :=
  ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → 0 ≤ D.mE T ∧
    D.mE T ≤ φ / Z.sc T * ((gapSet Z φ T).card : ℝ) + φ / Z.sc T

/-- The named analytic inputs for one Step-0 construction. -/
structure Inputs (Z : Zeros) (β : ℕ → ℝ) (ex lg : ℝ → ℝ) (φ θ : ℝ) (D : Data) : Prop where
  step3 : Step3 D
  main : MainTerm D
  step6 : Step6 Z β ex lg θ D
  step1 : Step1Pos D
  holder : Holder7 D
  step5 : Step5 D
  cover : Covering Z φ D

/-- Step 0: for every (φ, θ) with DPC's sign condition 0 < φθ < 1 < φ(1 + θ) < 2 (l.753) and θ < 1/2 (Step 5),
DPC's construction (l.743–771) has analytic quantities satisfying the named inputs. -/
def DPCInputs (Z : Zeros) (β : ℕ → ℝ) (ex lg : ℝ → ℝ) : Prop :=
  ∀ φ θ : ℝ, 0 < φ → φ < 1 → 0 < θ → θ < 1 / 2 → 0 < φ * θ → φ * θ < 1 → 1 < φ * (1 + θ) →
    φ * (1 + θ) < 2 → ∃ D : Data, Inputs Z β ex lg φ θ D

/-- **Theorem U** from the named inputs. -/
theorem theoremU {Z : Zeros} (hZ : Z.Ordered) {β : ℕ → ℝ} (hβ : Reflect Z β) {ex lg : ℝ → ℝ}
    (hEL : ExpLog ex lg) {tp : ℝ} (htp : 0 < tp) (hsc : ∀ T : ℝ, 0 < T → Z.sc T * tp = lg T)
    (hin : DPCInputs Z β ex lg) : TheoremU Z := by
  intro φ hφ1 hφ2
  have hφ0 : 0 < φ := by linarith
  obtain ⟨θ, hθ1, hθ2⟩ := step0_nonempty hφ1
  obtain ⟨hs1, hs2, hs3, hs4⟩ := step0_signs hφ1 hφ2 hθ1 hθ2
  have hθ0 : 0 < θ := by nlinarith
  obtain ⟨D, inp⟩ := hin φ θ hφ0 hφ2 hθ0 hθ2 hs1 hs2 hs3 hs4
  obtain ⟨δ, w, r, hδ0, hδ1, hw, hr, hκE⟩ := exists_params hθ2
  have hE0 : 0 < Eexp w r := by nlinarith
  obtain ⟨C, K, hK, h6⟩ := inp.step6 δ w r hδ0 hδ1 hw hr hκE
  have hc : 0 < Eexp w r / (2 * θ) := by positivity
  have hcθ : Eexp w r / (2 * θ) * θ < Eexp w r := by
    rw [div_mul_eq_mul_div, div_lt_iff₀ (by positivity)]
    nlinarith [mul_pos hE0 hθ0]
  have hΦ : ∀ τ : ℝ, 0 < τ → ∃ ε₁ : ℝ, 0 < ε₁ ∧ ∀ ε : ℝ, 0 < ε → ε ≤ ε₁ → ∀ m : ℕ,
      K * Phi ex lg C θ (1 - δ) (Eexp w r) (Eexp w r / (2 * θ)) ε m ≤ τ := by
    intro τ hτ
    obtain ⟨ε₁, hε₁, hb⟩ := Phi_small hEL C hE0.le hκE hc hcθ (τ / (K + 1)) (by positivity)
    refine ⟨ε₁, hε₁, fun ε hε hεε m => ?_⟩
    have h2 : K * Phi ex lg C θ (1 - δ) (Eexp w r) (Eexp w r / (2 * θ)) ε m ≤ K * (τ / (K + 1)) :=
      mul_le_mul_of_nonneg_left (hb ε hε hεε m) hK
    have h3 : K * (τ / (K + 1)) ≤ τ := by
      rw [mul_div_assoc', div_le_iff₀ (by linarith)]
      linarith
    linarith
  have hL : ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → 1 ≤ lg T := ⟨ex 1, fun T hT => hEL.le_lg hT⟩
  have h2 : ∀ T : ℝ, ((offSet Z β T).card : ℝ) ≤ 2 * ((gapSet Z φ T).card : ℝ) := fun T => by
    exact_mod_cast card_off_le hZ hβ hφ0 T
  obtain ⟨Tc, hTc⟩ := inp.cover
  have hcov : ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → 0 ≤ D.mE T ∧ 0 ≤ φ / Z.sc T ∧ φ / Z.sc T * lg T ≤ tp * φ ∧
      D.mE T ≤ φ / Z.sc T * ((gapSet Z φ T).card : ℝ) + φ / Z.sc T := by
    refine ⟨max Tc (ex 1), fun T hT => ?_⟩
    have hT1 : Tc ≤ T := le_trans (le_max_left _ _) hT
    have hT2 : ex 1 ≤ T := le_trans (le_max_right _ _) hT
    have hT0 : 0 < T := lt_of_lt_of_le (hEL.pos 1) hT2
    have hlg : 1 ≤ lg T := hEL.le_lg hT2
    have hsT : Z.sc T * tp = lg T := hsc T hT0
    have hs : 0 < Z.sc T := by
      by_contra hn
      have hn' : Z.sc T ≤ 0 := not_lt.1 hn
      nlinarith
    obtain ⟨hm0, hmE⟩ := hTc T hT1
    refine ⟨hm0, div_nonneg hφ0.le hs.le, ?_, hmE⟩
    rw [← hsT, ← mul_assoc, div_mul_cancel₀ _ hs.ne']
    linarith
  obtain ⟨ε₀, hε₀, T₀, hT₀⟩ := step7 (I := D.I) (M := D.M) (Hp := D.Hp) (S := D.S) (R4 := D.R4)
    (P4 := D.P4) (Bad := D.Bad) (mE := D.mE) (L := lg) (h := fun T => φ / Z.sc T)
    (Nc := fun T => ((gapSet Z φ T).card : ℝ)) (Off := fun T => ((offSet Z β T).card : ℝ))
    (Φ := fun ε m => K * Phi ex lg C θ (1 - δ) (Eexp w r) (Eexp w r / (2 * θ)) ε m) (H₀ := tp * φ)
    (mul_pos htp hφ0) hL (fun T => Nat.cast_nonneg _) inp.step3 inp.main (h6 _ hc hcθ) hΦ h2
    inp.step1 inp.holder inp.step5 hcov
  refine ⟨ε₀ * tp, mul_pos hε₀ htp, max T₀ 1, fun T hT => ?_⟩
  have hT1 : T₀ ≤ T := le_trans (le_max_left _ _) hT
  have hT0 : 0 < T := by linarith [le_max_right T₀ 1]
  have h := hT₀ T hT1
  have e : ε₀ * tp * T * Z.sc T = ε₀ * T * lg T := by rw [← hsc T hT0]; ring
  rw [e]
  exact h.le

end DPCLean

