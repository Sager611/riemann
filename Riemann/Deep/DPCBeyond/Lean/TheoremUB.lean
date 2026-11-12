import Riemann.Deep.DPCBeyond.Lean.DLamB
import Riemann.Deep.DPCBeyond.Lean.SumB
import Riemann.Deep.DPCBeyond.Lean.Step7B
import Riemann.Deep.DPCBeyond.Lean.ParamsB

/-!
# dpcbeyond: Theorem U′ from the named analytic inputs (certificate-interface pattern)

TheoremUB Z κ: for every φ ∈ (max{2/(2+κ), 1/2}, 1) there is c > 0 with 𝒩(T) ≥ c T sc T (sc T = log T/(2π)) for
large T. theoremUB: TheoremUB from the zero model (Zeros.Ordered, Reflect), ExpLog with sc T · 2π = log T, the
density slope κ (DensityK) and DPCInputsB: for every admissible (φ, θ) with θ < min{κ/2, 1}, the analytic quantities
DataB of DPC's Step-0 construction satisfy Theorem U's Step3, MainTerm, Step1Pos, Covering (dpclean, reused
verbatim on DataB.base) and the new HolderB (Step 7′), MomentM (Lemma M3), Step5cB (5c) and Step6ClassesB (6e′).
Proved here: Steps 0′/8′ and the order of choices (ParamsB), 6e′'s bookkeeping and Φ(ε) → 0 (HairpinB, SumB),
Step 7′ (Step7B), Step 2 (dpclean's card_off_le).
-/

namespace DPCBLean

open Finset PPLean DPCLean

variable {ex lg : ℝ → ℝ}

/-- **Theorem U′** (PROOF.md) in pplean's zero model: 𝒩(T) ≥ c(φ) T log T for φ ∈ (max{2/(2+κ), 1/2}, 1). -/
def TheoremUB (Z : Zeros) (κ : ℝ) : Prop := UFrom Z (max (2 / (2 + κ)) (1 / 2))

/-- The analytic quantities at fixed Step-0 data (φ, θ, P, f, η; PROOF Step 0), as functions of T: Theorem U's
Data with ∫|R_T|⁴, ∫|P(V)|⁴ replaced by the moment families Rq q T = ∫_T^{2T}|R_T|^q, Ps s T = ∫_T^{2T}|P(V_φ)|^s. -/
structure DataB where
  /-- Re ∫_T^{2T} G_T|R_T|² P(V_φ) -/
  I : ℝ → ℝ
  /-- Re M_T -/
  M : ℝ → ℝ
  /-- Re 𝓗_T, the hairpin total -/
  Hp : ℝ → ℝ
  /-- Σ_n |c_n|² -/
  S : ℝ → ℝ
  /-- Rq q T = ∫_T^{2T} |R_T|^q (Lemma M3) -/
  Rq : ℝ → ℝ → ℝ
  /-- Ps s T = ∫_T^{2T} |P(V_φ)|^s (5c) -/
  Ps : ℝ → ℝ → ℝ
  /-- ∫_{E_T} G_T |R_T|² |P(V_φ)| -/
  Bad : ℝ → ℝ
  /-- meas E_T -/
  mE : ℝ → ℝ

/-- Theorem U's Data of the same construction (R4 = Rq 4, P4 = Ps 4), for reusing Step3, MainTerm, Step1Pos,
Covering verbatim. -/
def DataB.base (D : DataB) : Data := ⟨D.I, D.M, D.Hp, D.S, D.Rq 4, D.Ps 4, D.Bad, D.mE⟩

/-- Step 7′, Hölder with exponents (q/2, s, 1/α), α = 1 − 2/q − 1/s, G_T bounded:
∫_{E_T} G_T|R_T|²|P(V)| ≤ K_G ‖R_T‖_q² ‖P(V)‖_s meas(E_T)^α. -/
def HolderB (ex lg : ℝ → ℝ) (D : DataB) : Prop :=
  ∀ q s : ℝ, 2 < q → 1 < s → 2 / q + 1 / s < 1 → ∃ KG T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T →
    D.Bad T ≤ KG * rp ex lg (D.Rq q T) (2 / q) * rp ex lg (D.Ps s T) (1 / s) *
      rp ex lg (D.mE T) (1 - 2 / q - 1 / s)

/-- Lemma M (M3): ∫_T^{2T}|R_T|^q ≪ T for every q < κ/θ (here q > 2), constant independent of ε. -/
def MomentM (κ θ : ℝ) (D : DataB) : Prop :=
  ∀ q : ℝ, 2 < q → q < κ / θ → ∃ K₁ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → D.Rq q T ≤ K₁ * T

/-- Step 5c: ∫_T^{2T}|P(V_φ)|^s ≪ T, i.e. ‖P(V)‖_s ≪ T^{1/s}, for every s > 1 (M. Riesz, Fujii, density). -/
def Step5cB (D : DataB) : Prop :=
  ∀ s : ℝ, 1 < s → ∃ K₂ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → D.Ps s T ≤ K₂ * T

/-- The named analytic inputs for one Step-0 construction (θ < min{κ/2, 1}). -/
structure InputsB (Z : Zeros) (β : ℕ → ℝ) (ex lg : ℝ → ℝ) (κ φ θ : ℝ) (D : DataB) : Prop where
  step3 : Step3 D.base
  main : MainTerm D.base
  step6 : Step6ClassesB Z β ex lg κ θ D.Hp
  step1 : Step1Pos D.base
  holder : HolderB ex lg D
  momM : MomentM κ θ D
  step5c : Step5cB D
  cover : Covering Z φ D.base

/-- Step 0′: for every (φ, θ) with DPC's sign condition 0 < φθ < 1 < φ(1 + θ) < 2 (l.753) and θ < min{κ/2, 1},
DPC's construction (l.743–771) has analytic quantities satisfying the named inputs. -/
def DPCInputsB (Z : Zeros) (β : ℕ → ℝ) (ex lg : ℝ → ℝ) (κ : ℝ) : Prop :=
  ∀ φ θ : ℝ, 0 < φ → φ < 1 → 0 < θ → θ < min (κ / 2) 1 → 0 < φ * θ → φ * θ < 1 →
    1 < φ * (1 + θ) → φ * (1 + θ) < 2 → ∃ D : DataB, InputsB Z β ex lg κ φ θ D

/-- Step 7's covering bound in step7's form (h = φ/sc T, h · log T = 2πφ); copied from dpclean's theoremU. -/
lemma cover_of {Z : Zeros} (hEL : ExpLog ex lg) {tp φ : ℝ} (htp : 0 < tp) (hφ0 : 0 < φ)
    (hsc : ∀ T : ℝ, 0 < T → Z.sc T * tp = lg T) {mE : ℝ → ℝ}
    (hcov : ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → 0 ≤ mE T ∧
      mE T ≤ φ / Z.sc T * ((gapSet Z φ T).card : ℝ) + φ / Z.sc T) :
    ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → 0 ≤ mE T ∧ 0 ≤ φ / Z.sc T ∧ φ / Z.sc T * lg T ≤ tp * φ ∧
      mE T ≤ φ / Z.sc T * ((gapSet Z φ T).card : ℝ) + φ / Z.sc T := by
  obtain ⟨Tc, hTc⟩ := hcov
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

/-- Lemma M3 and 5c at the chosen q, s, with one threshold. -/
lemma moments_of {D : DataB} {κ θ q s : ℝ} (hM : MomentM κ θ D) (h5 : Step5cB D) (hq : 2 < q)
    (hqκ : q < κ / θ) (hs : 1 < s) :
    ∃ K₁ K₂ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → D.Rq q T ≤ K₁ * T ∧ D.Ps s T ≤ K₂ * T := by
  obtain ⟨K₁, T₁, h₁⟩ := hM q hq hqκ
  obtain ⟨K₂, T₂, h₂⟩ := h5 s hs
  exact ⟨K₁, K₂, max T₁ T₂, fun T hT =>
    ⟨h₁ T (le_trans (le_max_left _ _) hT), h₂ T (le_trans (le_max_right _ _) hT)⟩⟩

/-- 6e′ (e)'s split parameter: c = (1 − 1/w)/(4θ) has c · 2θ < 1 − 1/w. -/
lemma split_c {θ w : ℝ} (hθ0 : 0 < θ) (hw : 1 < w) :
    0 < (1 - 1 / w) / (4 * θ) ∧ (1 - 1 / w) / (4 * θ) * (2 * θ) < 1 - 1 / w := by
  have hw0 : 0 < w := by linarith
  have hE1 : 0 < 1 - 1 / w := by
    have h1 : 1 / w < 1 := (div_lt_one hw0).2 hw
    linarith
  have e : (1 - 1 / w) / (4 * θ) * (2 * θ) = (1 - 1 / w) / 2 := by
    rw [div_mul_eq_mul_div, div_eq_div_iff (by positivity) (by norm_num)]
    ring
  exact ⟨by positivity, by rw [e]; linarith⟩

/-- **Theorem U′** from the named inputs. -/
theorem theoremUB {Z : Zeros} (hZ : Z.Ordered) {β : ℕ → ℝ} (hβ : Reflect Z β)
    (hEL : ExpLog ex lg) {tp : ℝ} (htp : 0 < tp) (hsc : ∀ T : ℝ, 0 < T → Z.sc T * tp = lg T)
    {κ : ℝ} (hκ : 0 < κ) (hden : DensityK Z β ex lg κ) (hin : DPCInputsB Z β ex lg κ) :
    TheoremUB Z κ := by
  intro φ hφ1 hφ2
  have hφ0 : 0 < φ := by
    have h := le_max_right (2 / (2 + κ)) (1 / 2)
    linarith
  obtain ⟨θ, hθ1, hθ2⟩ := step0B_nonempty hκ hφ1
  obtain ⟨hθ0, hθκ, -, hs1, hs2, hs3, hs4⟩ := step0B_signs hφ1 hφ2 hθ1 hθ2
  obtain ⟨D, inp⟩ := hin φ θ hφ0 hφ2 hθ0 hθ2 hs1 hs2 hs3 hs4
  obtain ⟨q, s, hq2, hqκ, hs1', hqs⟩ := exists_q_s hθ0 hθκ
  obtain ⟨u, w, hu1, -, -, hw1, huθ, hwθ, hv⟩ := exists_u_w hθ0 hθκ
  obtain ⟨C, K₁, K₂, hK₁, hK₂, h6⟩ := step6B_of_classes hEL hθ0.le hden inp.step6 hu1 hw1 hv huθ
  obtain ⟨hc0, hcθ⟩ := split_c hθ0 hw1
  have hΦ := PhiB_small hEL C hK₁ hK₂ hκ hw1 hv hwθ hc0 hcθ
  have hL : ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → 1 ≤ lg T := ⟨ex 1, fun T hT => hEL.le_lg hT⟩
  have h2 : ∀ T : ℝ, ((offSet Z β T).card : ℝ) ≤ 2 * ((gapSet Z φ T).card : ℝ) := fun T => by
    exact_mod_cast card_off_le hZ hβ hφ0 T
  obtain ⟨ε₀, hε₀, T₀, hT₀⟩ := step7B hEL (I := D.I) (M := D.M) (Hp := D.Hp) (S := D.S)
    (Rq := D.Rq q) (Ps := D.Ps s) (Bad := D.Bad) (mE := D.mE) (L := lg)
    (h := fun T => φ / Z.sc T) (Nc := fun T => ((gapSet Z φ T).card : ℝ))
    (Off := fun T => ((offSet Z β T).card : ℝ))
    (Φ := fun ε m => K₁ * Phi ex lg C (2 * θ) κ (1 - 1 / w) ((1 - 1 / w) / (4 * θ)) ε m +
      K₂ * Phi ex lg C 0 κ (1 - 1 / u - 1 / w) ((1 - 1 / w) / (4 * θ)) ε m) (H₀ := tp * φ)
    (mul_pos htp hφ0) hq2 (by linarith only [hs1']) hqs hL (fun T => Nat.cast_nonneg _) inp.step3
    inp.main (h6 ((1 - 1 / w) / (4 * θ))) hΦ h2 inp.step1 (inp.holder q s hq2 hs1' hqs)
    (moments_of inp.momM inp.step5c hq2 hqκ hs1') (cover_of hEL htp hφ0 hsc (mE := D.mE) inp.cover)
  refine ⟨ε₀ * tp, mul_pos hε₀ htp, max T₀ 1, fun T hT => ?_⟩
  have hT1 : T₀ ≤ T := le_trans (le_max_left _ _) hT
  have hT0 : 0 < T := lt_of_lt_of_le one_pos (le_trans (le_max_right _ _) hT)
  have h := hT₀ T hT1
  have e : ε₀ * tp * T * Z.sc T = ε₀ * T * lg T := by rw [← hsc T hT0]; ring
  rw [e]
  exact h.le

end DPCBLean

