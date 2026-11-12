import Riemann.Deep.DPCUncond.Lean.Fine
import Riemann.Deep.DPCBeyond.Lean.Rpow

/-!
# dpcbeyond Step 6e′ (new step (d)): hairpins through Lemma Z, |𝓗_T| ≤ Φ(ε)T + o(T)

PROOF.md 6e′: per depth class S_a, (a)–(b) 𝓗_{T,a} ≪ (a+2)^C L⁻¹ (e^{2θ(a+1)} Σ_{S_a} Y_ρ + Σ_{S_a} Off(x_ρ)² Y_ρ),
(c) Σ_ρ Y_ρ^w ≪ (a+2)^C TL and Σ_ρ Off(x_ρ)^{2u} ≪ TL (2uθ < κ, M4), (d) Hölder over S_a, giving
𝓗_{T,a} ≪ T(a+2)^C [e^{2θ(a+1)} σ_a^{1−1/w} + σ_a^{1/v}]; (e) the class sum is Φ(ε)T with Φ(ε) → 0 when
2θ < κ(1 − 1/w). Named here: DensityK (slope κ) and Step6ClassesB ((6a), (a)–(d) as inequalities); proved:
class_count (Step 2's output and the density bound |S_a|), class_bound (the exponent bookkeeping of (d)),
step6B_of_classes (dpclean's hairpin_sum twice: θ ↦ 2θ, E = 1 − 1/w and θ ↦ 0, E = 1/v) and PhiB_small
(dpclean's Phi_small twice). deepSet, classSet, offSet, Phi, hairpin_sum, Phi_small: dpclean (DPCUncond/Lean).
-/

namespace DPCBLean

open Finset PPLean DPCLean

variable {ex lg : ℝ → ℝ}

/-- Uniform zero density with slope κ (PROOF.md Theorem U′; Conrey κ = 8/7 − δ via Baluyot, Jutila κ = 1 − δ;
as quoted Lz 2611.14486 l.221–230): N(σ, 2T) ≪ T^{1−κ(σ−1/2)} log T at σ = 1/2 + a/log T gives
#{T ≤ γ ≤ 2T, β ≥ 1/2 + a/log T} ≤ C_d T log T e^{−κa}. -/
def DensityK (Z : Zeros) (β : ℕ → ℝ) (ex lg : ℝ → ℝ) (κ : ℝ) : Prop :=
  ∃ Cd T₀ : ℝ, 0 < Cd ∧ ∀ T : ℝ, T₀ ≤ T → ∀ a : ℕ,
    ((deepSet Z β lg T a).card : ℝ) ≤ Cd * T * lg T * ex (-(κ * a))

/-- Step 6e′ (a)–(d) at the level of the depth classes S_a = classSet (PROOF.md 6e′, notation of Theorem U's
6a–6d): for all u, w > 1 with 1/v := 1 − 1/u − 1/w > 0 and 2uθ < κ there are C, K and class quantities
Ha T a = 𝓗_{T,a}, Y₁ T a = Σ_{ρ∈S_a} Y_ρ, Y₂ T a = Σ_{ρ∈S_a} Off(x_ρ)² Y_ρ, MY T a = Σ_ρ Y_ρ^w (all zeros with
γ ∈ [T, 2T], segments of length Y_a = (a+1)/L), MO T = Σ_ρ Off(x_ρ)^{2u}, and a truncation M(T), such that for
large T: (6a) |𝓗_T| ≤ Σ_{a ≤ M(T)} 𝓗_{T,a} + o(T); (c) MO ≤ K T L and MY ≤ K (a+2)^C T L;
(a)–(b) 𝓗_{T,a} ≤ K (a+2)^C L⁻¹ (e^{2θ(a+1)} Y₁ + Y₂) (Lemma Z with W_{−y} ≤ On + X^{−y}Off and W ≤ On + Off,
the y-integral, summed over S_a); (d) Hölder over S_a: Y₁ ≤ |S_a|^{1−1/w} MY^{1/w},
Y₂ ≤ |S_a|^{1/v} MO^{1/u} MY^{1/w}. -/
def Step6ClassesB (Z : Zeros) (β : ℕ → ℝ) (ex lg : ℝ → ℝ) (κ θ : ℝ) (Hp : ℝ → ℝ) : Prop :=
  ∀ u w : ℝ, 1 < u → 1 < w → 0 < 1 - 1 / u - 1 / w → 2 * u * θ < κ → ∃ (C : ℕ) (K : ℝ), 0 ≤ K ∧
    ∃ (Ha Y₁ Y₂ MY : ℝ → ℕ → ℝ) (MO : ℝ → ℝ) (Mt : ℝ → ℕ), ∀ η : ℝ, 0 < η → ∃ T₀ : ℝ, ∀ T : ℝ,
      T₀ ≤ T → |Hp T| ≤ ∑ a ∈ range (Mt T + 1), Ha T a + η * T ∧ MO T ≤ K * (T * lg T) ∧
      ∀ a : ℕ, Ha T a ≤ K * ((a : ℝ) + 2) ^ C / lg T * (ex (2 * θ) ^ (a + 1) * Y₁ T a + Y₂ T a) ∧
        MY T a ≤ K * ((a : ℝ) + 2) ^ C * (T * lg T) ∧
        Y₁ T a ≤ rp ex lg ((classSet Z β lg T a).card : ℝ) (1 - 1 / w) * rp ex lg (MY T a) (1 / w) ∧
        Y₂ T a ≤ rp ex lg ((classSet Z β lg T a).card : ℝ) (1 - 1 / u - 1 / w) *
          rp ex lg (MO T) (1 / u) * rp ex lg (MY T a) (1 / w)

/-- 6a's count: |S_a| ≤ min(3ε, C_d e^{−κa}) T log T from #Z_off ≤ 2εTL (Step 2) and the density. -/
lemma class_count {Z : Zeros} {β : ℕ → ℝ} {κ Cd ε T : ℝ} (hε : 0 < ε) (hTL : 0 ≤ T * lg T)
    (hoff : ((offSet Z β T).card : ℝ) ≤ 2 * ε * T * lg T)
    (hd : ∀ a : ℕ, ((deepSet Z β lg T a).card : ℝ) ≤ Cd * T * lg T * ex (-(κ * a))) (a : ℕ) :
    ((classSet Z β lg T a).card : ℝ) ≤ min (3 * ε) (Cd * ex (-(κ * a))) * (T * lg T) := by
  have hc1 : ((classSet Z β lg T a).card : ℝ) ≤ 3 * ε * (T * lg T) := by
    have h2 : ((classSet Z β lg T a).card : ℝ) ≤ ((offSet Z β T).card : ℝ) := by
      exact_mod_cast card_le_card (filter_subset _ _)
    have h3 : 0 ≤ ε * (T * lg T) := mul_nonneg hε.le hTL
    have e1 : 2 * ε * T * lg T = 2 * (ε * (T * lg T)) := by ring
    have e2 : 3 * ε * (T * lg T) = 3 * (ε * (T * lg T)) := by ring
    linarith
  have hc2 : ((classSet Z β lg T a).card : ℝ) ≤ Cd * ex (-(κ * a)) * (T * lg T) := by
    have h1 : classSet Z β lg T a ⊆ deepSet Z β lg T a := by
      intro n hn
      simp only [classSet, deepSet, mem_filter] at hn ⊢
      exact ⟨hn.1, hn.2.1⟩
    have h2 : ((classSet Z β lg T a).card : ℝ) ≤ ((deepSet Z β lg T a).card : ℝ) := by
      exact_mod_cast card_le_card h1
    have h3 := hd a
    have e : Cd * T * lg T * ex (-(κ * a)) = Cd * ex (-(κ * a)) * (T * lg T) := by ring
    linarith
  rcases le_total (3 * ε) (Cd * ex (-(κ * a))) with h | h
  · rw [min_eq_left h]
    exact hc1
  · rw [min_eq_right h]
    exact hc2

/-- 6e′ (d), one class: Hölder over S_a, the moments (c) and |S_a| ≤ σTL give
𝓗_{T,a} ≤ K(K+1)² T (a+2)^{2C} [e^{2θ(a+1)} σ^{1−1/w} + σ^{1/v}] (B = (a+2)^C, G = e^{2θ(a+1)}). -/
lemma class_bound (h : ExpLog ex lg) {Ha Y₁ Y₂ MY MO n σ T L K B G w u : ℝ} (hT : 0 < T)
    (hL : 1 ≤ L) (hK : 0 ≤ K) (hB : 1 ≤ B) (hG : 0 ≤ G) (hσ : 0 < σ) (hw : 1 < w) (hu : 1 < u)
    (hv : 0 < 1 - 1 / u - 1 / w) (hn : n ≤ σ * (T * L)) (hMO : MO ≤ K * (T * L))
    (hMY : MY ≤ K * B * (T * L)) (hY₁ : Y₁ ≤ rp ex lg n (1 - 1 / w) * rp ex lg MY (1 / w))
    (hY₂ : Y₂ ≤ rp ex lg n (1 - 1 / u - 1 / w) * rp ex lg MO (1 / u) * rp ex lg MY (1 / w))
    (hHa : Ha ≤ K * B / L * (G * Y₁ + Y₂)) :
    Ha ≤ K * (K + 1) ^ 2 * T * B ^ 2 * G * rp ex lg σ (1 - 1 / w) +
      K * (K + 1) ^ 2 * T * B ^ 2 * rp ex lg σ (1 - 1 / u - 1 / w) := by
  have hL0 : 0 < L := by linarith only [hL]
  have hTL : 0 < T * L := mul_pos hT hL0
  have hK1 : 1 ≤ K + 1 := by linarith only [hK]
  have hK10 : 0 < K + 1 := by linarith only [hK]
  have hB0 : 0 < B := by linarith only [hB]
  have hw0 : 0 < w := by linarith only [hw]
  have hu0 : 0 < u := by linarith only [hu]
  have hiw : 0 ≤ 1 / w := by positivity
  have hiw1 : 1 / w ≤ 1 := by rw [div_le_one hw0]; exact hw.le
  have hiu : 0 ≤ 1 / u := by positivity
  have hiu1 : 1 / u ≤ 1 := by rw [div_le_one hu0]; exact hu.le
  have he1 : 0 ≤ 1 - 1 / w := by linarith only [hiw1]
  have hKB1 : 1 ≤ (K + 1) * B := one_le_mul_of_one_le_of_one_le hK1 hB
  have r1 : rp ex lg n (1 - 1 / w) ≤ rp ex lg σ (1 - 1 / w) * rp ex lg (T * L) (1 - 1 / w) := by
    rw [← rp_mul h hσ hTL]
    exact rp_le_rp h hn he1
  have r3 : rp ex lg n (1 - 1 / u - 1 / w) ≤
      rp ex lg σ (1 - 1 / u - 1 / w) * rp ex lg (T * L) (1 - 1 / u - 1 / w) := by
    rw [← rp_mul h hσ hTL]
    exact rp_le_rp h hn hv.le
  have r2 : rp ex lg MY (1 / w) ≤ (K + 1) * B * rp ex lg (T * L) (1 / w) := by
    have h1 : MY ≤ (K + 1) * B * (T * L) := by
      have h2 : K * B * (T * L) ≤ (K + 1) * B * (T * L) :=
        mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_right (by linarith only [hK]) hB0.le) hTL.le
      linarith only [hMY, h2]
    calc rp ex lg MY (1 / w) ≤ rp ex lg ((K + 1) * B * (T * L)) (1 / w) := rp_le_rp h h1 hiw
      _ = rp ex lg ((K + 1) * B) (1 / w) * rp ex lg (T * L) (1 / w) := rp_mul h (by positivity) hTL _
      _ ≤ (K + 1) * B * rp ex lg (T * L) (1 / w) := by
          apply mul_le_mul_of_nonneg_right _ (rp_nonneg h _ _)
          have h2 := rp_le_max h (x := (K + 1) * B) hiw hiw1
          rwa [max_eq_right hKB1] at h2
  have r4 : rp ex lg MO (1 / u) ≤ (K + 1) * rp ex lg (T * L) (1 / u) := by
    have h1 : MO ≤ (K + 1) * (T * L) := by
      have h2 : K * (T * L) ≤ (K + 1) * (T * L) :=
        mul_le_mul_of_nonneg_right (by linarith only [hK]) hTL.le
      linarith only [hMO, h2]
    calc rp ex lg MO (1 / u) ≤ rp ex lg ((K + 1) * (T * L)) (1 / u) := rp_le_rp h h1 hiu
      _ = rp ex lg (K + 1) (1 / u) * rp ex lg (T * L) (1 / u) := rp_mul h hK10 hTL _
      _ ≤ (K + 1) * rp ex lg (T * L) (1 / u) := by
          apply mul_le_mul_of_nonneg_right _ (rp_nonneg h _ _)
          have h2 := rp_le_max h (x := K + 1) hiu hiu1
          rwa [max_eq_right hK1] at h2
  have eTL1 : rp ex lg (T * L) (1 - 1 / w) * rp ex lg (T * L) (1 / w) = T * L := by
    rw [← rp_add h hTL, show 1 - 1 / w + 1 / w = 1 by ring, rp_one h hTL]
  have eTL2 : rp ex lg (T * L) (1 - 1 / u - 1 / w) * rp ex lg (T * L) (1 / u) *
      rp ex lg (T * L) (1 / w) = T * L := by
    rw [← rp_add h hTL, ← rp_add h hTL, show 1 - 1 / u - 1 / w + 1 / u + 1 / w = 1 by ring,
      rp_one h hTL]
  obtain ⟨P₁, hP₁⟩ : ∃ x : ℝ, x = rp ex lg σ (1 - 1 / w) := ⟨_, rfl⟩
  obtain ⟨P₃, hP₃⟩ : ∃ x : ℝ, x = rp ex lg σ (1 - 1 / u - 1 / w) := ⟨_, rfl⟩
  rw [← hP₁] at r1
  rw [← hP₃] at r3
  rw [← hP₁, ← hP₃]
  have hP₁0 : 0 ≤ P₁ := by rw [hP₁]; exact rp_nonneg h _ _
  have hP₃0 : 0 ≤ P₃ := by rw [hP₃]; exact rp_nonneg h _ _
  have hY₁' : Y₁ ≤ (K + 1) * B * P₁ * (T * L) := by
    calc Y₁ ≤ rp ex lg n (1 - 1 / w) * rp ex lg MY (1 / w) := hY₁
      _ ≤ (P₁ * rp ex lg (T * L) (1 - 1 / w)) * ((K + 1) * B * rp ex lg (T * L) (1 / w)) :=
          mul_le_mul r1 r2 (rp_nonneg h _ _) (mul_nonneg hP₁0 (rp_nonneg h _ _))
      _ = (K + 1) * B * P₁ * (rp ex lg (T * L) (1 - 1 / w) * rp ex lg (T * L) (1 / w)) := by ring
      _ = (K + 1) * B * P₁ * (T * L) := by rw [eTL1]
  have hY₂' : Y₂ ≤ (K + 1) ^ 2 * B * P₃ * (T * L) := by
    calc Y₂ ≤ rp ex lg n (1 - 1 / u - 1 / w) * rp ex lg MO (1 / u) * rp ex lg MY (1 / w) := hY₂
      _ ≤ (P₃ * rp ex lg (T * L) (1 - 1 / u - 1 / w)) * ((K + 1) * rp ex lg (T * L) (1 / u)) *
          ((K + 1) * B * rp ex lg (T * L) (1 / w)) := by
          apply mul_le_mul (mul_le_mul r3 r4 (rp_nonneg h _ _)
            (mul_nonneg hP₃0 (rp_nonneg h _ _))) r2 (rp_nonneg h _ _)
          exact mul_nonneg (mul_nonneg hP₃0 (rp_nonneg h _ _)) (mul_nonneg hK10.le (rp_nonneg h _ _))
      _ = (K + 1) ^ 2 * B * P₃ * (rp ex lg (T * L) (1 - 1 / u - 1 / w) * rp ex lg (T * L) (1 / u) *
          rp ex lg (T * L) (1 / w)) := by ring
      _ = (K + 1) ^ 2 * B * P₃ * (T * L) := by rw [eTL2]
  have hGY : G * Y₁ ≤ G * ((K + 1) * B * P₁ * (T * L)) := mul_le_mul_of_nonneg_left hY₁' hG
  have hsum : G * Y₁ + Y₂ ≤ ((K + 1) * B * G * P₁ * T + (K + 1) ^ 2 * B * P₃ * T) * L := by
    have e : ((K + 1) * B * G * P₁ * T + (K + 1) ^ 2 * B * P₃ * T) * L =
        G * ((K + 1) * B * P₁ * (T * L)) + (K + 1) ^ 2 * B * P₃ * (T * L) := by ring
    rw [e]
    linarith only [hGY, hY₂']
  have hKBL : 0 ≤ K * B / L := by positivity
  have h1 := le_trans hHa (mul_le_mul_of_nonneg_left hsum hKBL)
  have e : K * B / L * (((K + 1) * B * G * P₁ * T + (K + 1) ^ 2 * B * P₃ * T) * L) =
      K * B * ((K + 1) * B * G * P₁ * T + (K + 1) ^ 2 * B * P₃ * T) := by
    rw [← mul_assoc, div_mul_eq_mul_div, div_mul_cancel₀ _ hL0.ne']
  rw [e] at h1
  have hK2 : K + 1 ≤ (K + 1) ^ 2 := by
    rw [sq]
    exact le_mul_of_one_le_right hK10.le hK1
  have hc : 0 ≤ K * B * B * G * P₁ * T := by positivity
  have h2 : K * B * B * G * P₁ * T * (K + 1) ≤ K * B * B * G * P₁ * T * (K + 1) ^ 2 :=
    mul_le_mul_of_nonneg_left hK2 hc
  have e2 : K * B * ((K + 1) * B * G * P₁ * T + (K + 1) ^ 2 * B * P₃ * T) =
      K * B * B * G * P₁ * T * (K + 1) + K * (K + 1) ^ 2 * T * B ^ 2 * P₃ := by ring
  have e3 : K * (K + 1) ^ 2 * T * B ^ 2 * G * P₁ = K * B * B * G * P₁ * T * (K + 1) ^ 2 := by ring
  linarith only [h1, h2, e2, e3]

end DPCBLean
