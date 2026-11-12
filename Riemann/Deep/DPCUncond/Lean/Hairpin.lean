import Riemann.Deep.DPCUncond.Lean.Phi

/-!
# dpcuncond Step 6f: the hairpin sum over depth classes (proved)

hairpin_sum: per-class bounds 𝓗_{T,a} ≤ K T (a + 2)^C e^{θ(a+1)} s_a^E (the output of 6e, holder_bookkeeping)
with s_a ≤ 3ε (Step 2 and the contradiction hypothesis) and s_a ≤ C_d e^{−κa} (density) give
Σ_{a≤M} 𝓗_{T,a} ≤ K max(1, C_d^E) Φ(ε) T, Φ = Phi with exponent C + 1 (the a ≤ A(ε) terms use 3ε, the others the
density). This is PROOF 6f's first display.
-/

namespace DPCLean

open Finset

variable {ex lg : ℝ → ℝ}

lemma ExpLog.lg_ex (h : ExpLog ex lg) (x : ℝ) : lg (ex x) = x :=
  h.mono.injective (h.ex_lg _ (h.pos x))

/-- **6f summation**. -/
theorem hairpin_sum (h : ExpLog ex lg) {C : ℕ} {θ κ E c ε K T Cd : ℝ} (M : ℕ) (H s : ℕ → ℝ)
    (hθ : 0 ≤ θ) (hE : 0 ≤ E) (hK : 0 ≤ K) (hT : 0 ≤ T) (hCd : 0 < Cd)
    (hs0 : ∀ a, 0 < s a) (hsε : ∀ a, s a ≤ 3 * ε) (hsd : ∀ a : ℕ, s a ≤ Cd * ex (-(κ * a)))
    (hH : ∀ a : ℕ, H a ≤ K * T * ((a : ℝ) + 2) ^ C * ex θ ^ (a + 1) * ex (lg (s a) * E)) :
    ∑ a ∈ range (M + 1), H a ≤ K * max 1 (ex (lg Cd * E)) * Phi ex lg (C + 1) θ κ E c ε M * T := by
  have hq1 : 1 ≤ ex θ := by
    have h1 := h.mono.monotone hθ
    rwa [h.zero] at h1
  have hq0 : 0 ≤ ex θ := (h.pos θ).le
  have hKT : 0 ≤ K * T := mul_nonneg hK hT
  have hm1 : (1 : ℝ) ≤ max 1 (ex (lg Cd * E)) := le_max_left _ _
  have hρ : ex θ * ex (-(κ * E)) = ex (θ - κ * E) := by rw [sub_eq_add_neg, h.add]
  obtain ⟨A, hA⟩ : ∃ A : ℕ, A = Aeps lg c ε := ⟨_, rfl⟩
  obtain ⟨X, hX⟩ : ∃ X : ℝ, X = K * T * ((A : ℝ) + 2) ^ C * ex θ ^ (A + 1) * ex (lg (3 * ε) * E) :=
    ⟨_, rfl⟩
  have hX0 : 0 ≤ X := by
    rw [hX]
    exact mul_nonneg (mul_nonneg (mul_nonneg hKT (by positivity)) (pow_nonneg hq0 _)) (h.pos _).le
  -- the classes a ≤ A: s_a ≤ 3ε
  have hlow : ∀ a : ℕ, a ≤ A → H a ≤ X := by
    intro a ha
    have ha' : (a : ℝ) ≤ A := by exact_mod_cast ha
    have h1 : ((a : ℝ) + 2) ^ C ≤ ((A : ℝ) + 2) ^ C := pow_le_pow_left₀ (by positivity) (by linarith) C
    have h2 : ex θ ^ (a + 1) ≤ ex θ ^ (A + 1) := pow_le_pow_right₀ hq1 (by omega)
    have h3 : ex (lg (s a) * E) ≤ ex (lg (3 * ε) * E) := by
      apply h.mono.monotone
      apply mul_le_mul_of_nonneg_right _ hE
      exact h.le_lg (by rw [h.ex_lg _ (hs0 a)]; exact hsε a)
    have hb1 : 0 ≤ K * T * ((A : ℝ) + 2) ^ C := mul_nonneg hKT (by positivity)
    have hb2 : 0 ≤ K * T * ((A : ℝ) + 2) ^ C * ex θ ^ (A + 1) := mul_nonneg hb1 (pow_nonneg hq0 _)
    rw [hX]
    calc H a ≤ K * T * ((a : ℝ) + 2) ^ C * ex θ ^ (a + 1) * ex (lg (s a) * E) := hH a
      _ ≤ K * T * ((A : ℝ) + 2) ^ C * ex θ ^ (A + 1) * ex (lg (3 * ε) * E) :=
        mul_le_mul (mul_le_mul (mul_le_mul_of_nonneg_left h1 hKT) h2 (pow_nonneg hq0 _) hb1) h3
          (h.pos _).le hb2
  -- the classes a > A: s_a ≤ C_d e^{−κa}
  have hhigh : ∀ a : ℕ, H a ≤ K * T * max 1 (ex (lg Cd * E)) *
      (((a : ℝ) + 2) ^ (C + 1) * ex θ * ex (θ - κ * E) ^ a) := by
    intro a
    have h3 : ex (lg (s a) * E) ≤ ex (lg Cd * E) * ex (-(κ * E)) ^ a := by
      have hl1 : lg (s a) ≤ lg (Cd * ex (-(κ * a))) :=
        h.le_lg (by rw [h.ex_lg _ (hs0 a)]; exact hsd a)
      rw [h.lg_mul hCd (h.pos _), h.lg_ex] at hl1
      have hl2 : lg (s a) * E ≤ lg Cd * E + (a : ℝ) * (-(κ * E)) := by
        have h4 := mul_le_mul_of_nonneg_right hl1 hE
        linarith
      calc ex (lg (s a) * E) ≤ ex (lg Cd * E + (a : ℝ) * (-(κ * E))) := h.mono.monotone hl2
        _ = ex (lg Cd * E) * ex (-(κ * E)) ^ a := by rw [h.add, h.nat_mul]
    have hm : ex (lg Cd * E) ≤ max 1 (ex (lg Cd * E)) := le_max_right _ _
    have ha0 : (1 : ℝ) ≤ (a : ℝ) + 2 := by
      have h4 : (0 : ℝ) ≤ a := Nat.cast_nonneg a
      linarith
    have hpow : ((a : ℝ) + 2) ^ C ≤ ((a : ℝ) + 2) ^ (C + 1) := pow_le_pow_right₀ ha0 (by omega)
    have hb0 : 0 ≤ K * T * ((a : ℝ) + 2) ^ C * ex θ ^ (a + 1) :=
      mul_nonneg (mul_nonneg hKT (by positivity)) (pow_nonneg hq0 _)
    have hρ0 : 0 ≤ ex (θ - κ * E) ^ a := pow_nonneg (h.pos _).le a
    calc H a ≤ K * T * ((a : ℝ) + 2) ^ C * ex θ ^ (a + 1) * ex (lg (s a) * E) := hH a
      _ ≤ K * T * ((a : ℝ) + 2) ^ C * ex θ ^ (a + 1) * (ex (lg Cd * E) * ex (-(κ * E)) ^ a) :=
          mul_le_mul_of_nonneg_left h3 hb0
      _ = K * T * ex (lg Cd * E) * (((a : ℝ) + 2) ^ C * ex θ * (ex θ * ex (-(κ * E))) ^ a) := by
          rw [pow_succ, mul_pow]; ring
      _ = K * T * ex (lg Cd * E) * (((a : ℝ) + 2) ^ C * ex θ * ex (θ - κ * E) ^ a) := by rw [hρ]
      _ ≤ K * T * max 1 (ex (lg Cd * E)) * (((a : ℝ) + 2) ^ (C + 1) * ex θ * ex (θ - κ * E) ^ a) := by
          apply mul_le_mul (mul_le_mul_of_nonneg_left hm hKT) _ (by positivity) (by positivity)
          exact mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_right hpow hq0) hρ0
  -- split the a-sum at A
  have hsplit := (sum_filter_add_sum_filter_not (range (M + 1)) (fun a => a ≤ A) H).symm
  have hIoc : (range (M + 1)).filter (fun a => ¬ a ≤ A) = Ioc A M := by
    ext a
    simp only [mem_filter, mem_range, mem_Ioc]
    omega
  rw [hsplit, hIoc]
  have hfirst : ∑ a ∈ (range (M + 1)).filter (fun a => a ≤ A), H a ≤ ((A : ℝ) + 1) * X := by
    have h1 := sum_le_card_nsmul ((range (M + 1)).filter (fun a => a ≤ A)) H X
      (fun a ha => hlow a (mem_filter.1 ha).2)
    have h2 : ((range (M + 1)).filter (fun a => a ≤ A)).card ≤ A + 1 := by
      calc ((range (M + 1)).filter (fun a => a ≤ A)).card ≤ (range (A + 1)).card := by
            apply card_le_card
            intro a ha
            simp only [mem_filter, mem_range] at ha ⊢
            omega
        _ = A + 1 := card_range _
    rw [nsmul_eq_mul] at h1
    have h3 : (((range (M + 1)).filter (fun a => a ≤ A)).card : ℝ) ≤ (A : ℝ) + 1 := by
      exact_mod_cast h2
    calc ∑ a ∈ (range (M + 1)).filter (fun a => a ≤ A), H a
        ≤ (((range (M + 1)).filter (fun a => a ≤ A)).card : ℝ) * X := h1
      _ ≤ ((A : ℝ) + 1) * X := mul_le_mul_of_nonneg_right h3 hX0
  have hsecond : ∑ a ∈ Ioc A M, H a ≤ K * T * max 1 (ex (lg Cd * E)) *
      ∑ a ∈ Ioc A M, ((a : ℝ) + 2) ^ (C + 1) * ex θ * ex (θ - κ * E) ^ a := by
    rw [mul_sum]
    exact sum_le_sum (fun a _ => hhigh a)
  have hQ0 : 0 ≤ ex θ ^ (A + 1) * ex (lg (3 * ε) * E) := mul_nonneg (pow_nonneg hq0 _) (h.pos _).le
  have hY0 : 0 ≤ ((A : ℝ) + 2) ^ (C + 1) * ex θ ^ (A + 1) * ex (lg (3 * ε) * E) :=
    mul_nonneg (mul_nonneg (by positivity) (pow_nonneg hq0 _)) (h.pos _).le
  have hX' : ((A : ℝ) + 1) * X ≤ K * T * max 1 (ex (lg Cd * E)) *
      (((A : ℝ) + 2) ^ (C + 1) * ex θ ^ (A + 1) * ex (lg (3 * ε) * E)) := by
    have hAA : ((A : ℝ) + 1) * ((A : ℝ) + 2) ^ C ≤ ((A : ℝ) + 2) ^ (C + 1) := by
      rw [pow_succ]
      calc ((A : ℝ) + 1) * ((A : ℝ) + 2) ^ C ≤ ((A : ℝ) + 2) * ((A : ℝ) + 2) ^ C :=
            mul_le_mul_of_nonneg_right (by linarith) (by positivity)
        _ = ((A : ℝ) + 2) ^ C * ((A : ℝ) + 2) := by ring
    calc ((A : ℝ) + 1) * X
        = K * T * (((A : ℝ) + 1) * ((A : ℝ) + 2) ^ C) * (ex θ ^ (A + 1) * ex (lg (3 * ε) * E)) := by
          rw [hX]; ring
      _ ≤ K * T * ((A : ℝ) + 2) ^ (C + 1) * (ex θ ^ (A + 1) * ex (lg (3 * ε) * E)) :=
          mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_left hAA hKT) hQ0
      _ = (K * T) * (((A : ℝ) + 2) ^ (C + 1) * ex θ ^ (A + 1) * ex (lg (3 * ε) * E)) := by ring
      _ ≤ (K * T * max 1 (ex (lg Cd * E))) *
            (((A : ℝ) + 2) ^ (C + 1) * ex θ ^ (A + 1) * ex (lg (3 * ε) * E)) :=
          mul_le_mul_of_nonneg_right (le_mul_of_one_le_right hKT hm1) hY0
  unfold Phi
  rw [← hA]
  have e : K * max 1 (ex (lg Cd * E)) * (((A : ℝ) + 2) ^ (C + 1) * ex θ ^ (A + 1) * ex (lg (3 * ε) * E)
      + ∑ a ∈ Ioc A M, ((a : ℝ) + 2) ^ (C + 1) * ex θ * ex (θ - κ * E) ^ a) * T
      = K * T * max 1 (ex (lg Cd * E)) * (((A : ℝ) + 2) ^ (C + 1) * ex θ ^ (A + 1) * ex (lg (3 * ε) * E))
        + K * T * max 1 (ex (lg Cd * E)) *
          ∑ a ∈ Ioc A M, ((a : ℝ) + 2) ^ (C + 1) * ex θ * ex (θ - κ * E) ^ a := by ring
  rw [e]
  linarith

end DPCLean

