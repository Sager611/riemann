/-
E1Lean.Parity — SL Thm 5.1 (Zenodo 23134787 v0.03, §5) deduced from its inputs: the parity product
[SL ref 4, Prop 2.5] per window, Wang's pair-sum limit (arXiv 2609.07918) and the optimal-density limit
inf C(f) = 2 − c(θ); with R(2 − c, 1 − k) = h(θ; k).
-/
import Riemann.Deep.SantiAudit.Lean.Arch
import Mathlib.Tactic.LinearCombination

namespace E1Lean

variable (F : Fns)

/-- SL §5: R(A, b) = (4 − b − √(b(b + 8(A − 1))))/4, the smaller root of 2(1 − s)² = (A − s)b. -/
noncomputable def Rroot (A b : ℝ) : ℝ := (4 - b - F.sqrt (b * (b + 8 * (A - 1)))) / 4

/-- SL §5: R(2 − c(θ), 1 − k) = h(θ; k). -/
theorem Rroot_eq (c k : ℝ) : Rroot F (2 - c) (1 - k) = hc F c k := by
  unfold Rroot hc
  have e : (1 - k) * ((1 - k) + 8 * ((2 - c) - 1)) = (1 - k) * (9 - k - 8 * c) := by ring
  rw [e]; ring

/-- Root step: 2u² ≤ (A′ − 1 + u)b′ with u ≤ 1 and a small perturbation (A′, b′) of (A, b) ⟹ u < 1 − R(A, b) + η. -/
theorem root_step {A b A' b' u η : ℝ} (hb : 0 < b) (hD0 : 0 ≤ b * (b + 8 * (A - 1))) (hη : 0 < η)
    (hu1 : u ≤ 1) (hbb : b ≤ b') (hpert : (b' - b) + (b' * (A' - 1) - b * (A - 1)) < 2 * η ^ 2)
    (hq : 2 * u ^ 2 ≤ (A' - 1 + u) * b') : u < 1 - Rroot F A b + η := by
  have hr := F.sqrt_nonneg (b * (b + 8 * (A - 1)))
  have hss := F.mul_self_sqrt _ hD0
  have hu0 : 1 - Rroot F A b = (b + F.sqrt (b * (b + 8 * (A - 1)))) / 4 := by unfold Rroot; ring
  rw [hu0]
  by_contra hcon
  rw [not_lt] at hcon
  have hfac : 2 * u ^ 2 - b * u - b * (A - 1) = 2 * (u - (b + F.sqrt (b * (b + 8 * (A - 1)))) / 4) *
      (u - (b - F.sqrt (b * (b + 8 * (A - 1)))) / 4) := by linear_combination (1 / 8 : ℝ) * hss
  have h1 : η ≤ u - (b + F.sqrt (b * (b + 8 * (A - 1)))) / 4 := by linarith
  have h2 : η ≤ u - (b - F.sqrt (b * (b + 8 * (A - 1)))) / 4 := by linarith
  have hg : 2 * η ^ 2 ≤ 2 * u ^ 2 - b * u - b * (A - 1) := by
    rw [hfac]; nlinarith [mul_le_mul h1 h2 hη.le (hη.le.trans h1)]
  have hu0' : 0 ≤ u := by nlinarith
  have hbu : (b' - b) * u ≤ b' - b := by nlinarith [mul_nonneg (sub_nonneg.2 hbb) (sub_nonneg.2 hu1)]
  linarith

/-- Per window: parity product (ratio form), q ≤ A′, 1 − o ≤ b′, s ≥ 0, o ≤ 1 ⟹ s > R(A, b) − η. -/
theorem per_T {A b A' b' q s o η : ℝ} (hb : 0 < b) (hD0 : 0 ≤ b * (b + 8 * (A - 1))) (hη : 0 < η)
    (hA1 : 1 ≤ A') (hs0 : 0 ≤ s) (ho1 : o ≤ 1) (hqA : q ≤ A') (hob : 1 - o ≤ b') (hbb : b ≤ b')
    (hpert : (b' - b) + (b' * (A' - 1) - b * (A - 1)) < 2 * η ^ 2)
    (hpar : 2 * (1 - s) ^ 2 ≤ (q - s) * (1 - o)) : Rroot F A b - η < s := by
  have hR1 : Rroot F A b ≤ 1 := by
    unfold Rroot; have := F.sqrt_nonneg (b * (b + 8 * (A - 1))); linarith
  by_cases hs1 : 1 ≤ s
  · linarith
  · rw [not_le] at hs1
    have hu : 0 < 1 - s := by linarith
    have hpos : 0 < (q - s) * (1 - o) := lt_of_lt_of_le (mul_pos two_pos (pow_pos hu 2)) hpar
    have ho : 0 < 1 - o := by
      rcases (sub_nonneg.2 ho1).lt_or_eq with hlt | heq
      · exact hlt
      · rw [← heq, mul_zero] at hpos; exact absurd hpos (lt_irrefl 0)
    have hqs : 0 < q - s := by nlinarith
    have hq : 2 * (1 - s) ^ 2 ≤ (A' - 1 + (1 - s)) * b' := by
      have m1 : (q - s) * (1 - o) ≤ (A' - s) * (1 - o) := mul_le_mul_of_nonneg_right (by linarith) ho.le
      have m2 : (A' - s) * (1 - o) ≤ (A' - s) * b' := mul_le_mul_of_nonneg_left hob (by linarith)
      nlinarith
    have := root_step F hb hD0 hη (by linarith) hbb hpert hq
    linarith

/-- Prop 2.5 of SL ref 4 in counts ⟹ ratio form (N > 0). -/
theorem ratio_parity {Q S O N : ℝ} (hN : 0 < N) (hp : 2 * (N - S) ^ 2 ≤ (Q - S) * (N - O)) :
    2 * (1 - S / N) ^ 2 ≤ (Q / N - S / N) * (1 - O / N) := by
  have hN' : N ≠ 0 := hN.ne'
  have l : 2 * (1 - S / N) ^ 2 = 2 * (N - S) ^ 2 / N ^ 2 := by field_simp
  have r : (Q / N - S / N) * (1 - O / N) = (Q - S) * (N - O) / N ^ 2 := by field_simp
  rw [l, r]
  exact div_le_div_of_nonneg_right hp (sq_nonneg N)

/-- **SL Thm 5.1** (§5) from: window counts (N > 0, S ≥ 0, O ≤ N eventually); the parity product
  (Q_f − S)(N − O) ≥ 2(N − S)² eventually, for every test function f of the family ι [SL ref 4, Prop 2.5];
  Wang's limsup Q_f/N ≤ C(f) (arXiv 2609.07918, Thm 2.2 + (eq:unweighted)); inf_f C(f) = 2 − c(θ)
  (cosine density at λ, then λ ↑ θ). Conclusion: liminf O/N ≥ k ⟹ liminf S/N ≥ h(θ; k). -/
theorem sl51 (I : Inputs) {ι : Type} (Qf : ι → ℝ → ℝ → ℝ) (C : ι → ℝ) {θ k : ℝ}
    (hθ0 : 0 < θ) (hθ1 : θ ≤ 1) (hk1 : k < 1)
    (hcount : ∃ T₀ : ℝ, ∀ T ≥ T₀, 0 < I.N T (F.rpow T θ) ∧ 0 ≤ I.S T (F.rpow T θ) ∧
      I.O T (F.rpow T θ) ≤ I.N T (F.rpow T θ))
    (hpar : ∀ f, ∃ T₀ : ℝ, ∀ T ≥ T₀, 2 * (I.N T (F.rpow T θ) - I.S T (F.rpow T θ)) ^ 2 ≤
      (Qf f T (F.rpow T θ) - I.S T (F.rpow T θ)) * (I.N T (F.rpow T θ) - I.O T (F.rpow T θ)))
    (hwang : ∀ f, ∀ ε > 0, ∃ T₀ : ℝ, ∀ T ≥ T₀, ratio F (Qf f) I.N θ T ≤ C f + ε)
    (hopt : ∀ δ > 0, ∃ f, C f ≤ 2 - cW F θ + δ)
    (hO : LiminfGe (ratio F I.O I.N θ) k) : LiminfGe (ratio F I.S I.N θ) (h F θ k) := by
  intro η hη
  have hc1 := cW_lt_one F hθ0 hθ1
  set A := 2 - cW F θ with hAdef
  set b := 1 - k with hbdef
  have hb : 0 < b := by linarith
  have hA : 1 < A := by linarith
  have hD0 : 0 ≤ b * (b + 8 * (A - 1)) := mul_nonneg hb.le (by linarith)
  set M := 3 + 2 * b + (A - 1) with hM
  have hM0 : 0 < M := by linarith
  set ε := min 1 (η ^ 2 / M) with hε
  have hε0 : 0 < ε := lt_min one_pos (div_pos (pow_pos hη 2) hM0)
  have hε1 : ε ≤ 1 := min_le_left _ _
  have hεM : ε * M ≤ η ^ 2 := by
    have := min_le_right 1 (η ^ 2 / M)
    rw [← hε] at this
    calc ε * M ≤ η ^ 2 / M * M := mul_le_mul_of_nonneg_right this hM0.le
      _ = η ^ 2 := div_mul_cancel₀ _ hM0.ne'
  obtain ⟨f, hf⟩ := hopt ε hε0
  obtain ⟨T1, hT1⟩ := hcount
  obtain ⟨T2, hT2⟩ := hwang f ε hε0
  obtain ⟨T3, hT3⟩ := hO ε hε0
  obtain ⟨T4, hT4⟩ := hpar f
  refine ⟨max (max (max T1 T2) T3) T4, fun T hT => ?_⟩
  have hT' : max (max T1 T2) T3 ≤ T := le_trans (le_max_left _ _) hT
  have h1 := hT1 T (le_trans (le_trans (le_max_left _ _) (le_max_left _ _)) hT')
  have h2 := hT2 T (le_trans (le_trans (le_max_right _ _) (le_max_left _ _)) hT')
  have h3 := hT3 T (le_trans (le_max_right _ _) hT')
  obtain ⟨hN, hS, hON⟩ := h1
  have hpr := ratio_parity hN (hT4 T (le_trans (le_max_right _ _) hT))
  have hpert : (b + ε - b) + ((b + ε) * (A + 2 * ε - 1) - b * (A - 1)) < 2 * η ^ 2 := by
    nlinarith [mul_le_mul_of_nonneg_left hε1 hε0.le]
  have key := per_T F (A' := A + 2 * ε) (b' := b + ε) (q := ratio F (Qf f) I.N θ T)
    (s := ratio F I.S I.N θ T) (o := ratio F I.O I.N θ T) hb hD0 hη (by linarith)
    (div_nonneg hS hN.le) ((div_le_one hN).2 hON) (by linarith) (by linarith) (by linarith) hpert hpr
  rw [hAdef, hbdef, Rroot_eq] at key
  unfold h
  linarith

/-- The inputs of SL Thm 5.1 at θ: a test-function family ι (even η ∈ C_c^∞(−λ/2, λ/2), λ < θ, f = η²) with
  window pair sums Q_f(T, H) and Wang limits C(f), satisfying counts, parity product, Wang, optimum. -/
def SL51Inputs (I : Inputs) (θ : ℝ) : Prop :=
  ∃ (ι : Type) (Qf : ι → ℝ → ℝ → ℝ) (C : ι → ℝ),
    (∃ T₀ : ℝ, ∀ T ≥ T₀, 0 < I.N T (F.rpow T θ) ∧ 0 ≤ I.S T (F.rpow T θ) ∧
      I.O T (F.rpow T θ) ≤ I.N T (F.rpow T θ)) ∧
    (∀ f, ∃ T₀ : ℝ, ∀ T ≥ T₀, 2 * (I.N T (F.rpow T θ) - I.S T (F.rpow T θ)) ^ 2 ≤
      (Qf f T (F.rpow T θ) - I.S T (F.rpow T θ)) * (I.N T (F.rpow T θ) - I.O T (F.rpow T θ))) ∧
    (∀ f, ∀ ε > 0, ∃ T₀ : ℝ, ∀ T ≥ T₀, ratio F (Qf f) I.N θ T ≤ C f + ε) ∧
    (∀ δ > 0, ∃ f, C f ≤ 2 - cW F θ + δ)

/-- H3 (SL Thm 5.1) follows from its inputs at every θ ∈ (1/2, 1). -/
theorem H3_of_inputs (I : Inputs) (hin : ∀ θ : ℝ, 1 / 2 < θ → θ < 1 → SL51Inputs F I θ) : H3 F I := by
  intro θ k hθ1 hθ2 _ hk1 hO
  obtain ⟨ι, Qf, C, hcount, hpar, hwang, hopt⟩ := hin θ hθ1 hθ2
  exact sl51 F I Qf C (by linarith) hθ2.le hk1 hcount hpar hwang hopt hO

/-- **Theorem E1** with SL Thm 5.1 replaced by its inputs (parity product, Wang, optimum, counts). -/
theorem theoremE1_parity (I : Inputs) (h1 : H1 I) (h2 : H2 F I) (hin : ∀ θ : ℝ, 1 / 2 < θ → θ < 1 → SL51Inputs F I θ)
    (h4 : H4 I) (θ : ℝ) (hθ1 : 1 / 2 < θ) (hθ2 : θ < 1) :
    LiminfGe (ratio F I.S I.N θ) (h F θ (κ₁ F I)) ∧ h F (1 / 2) (κ₁ F I) ≤ h F θ (κ₁ F I) ∧
      (0.0055 : ℝ) < h F (1 / 2) (κ₁ F I) :=
  theoremE1 F I h1 h2 (H3_of_inputs F I hin) h4 θ hθ1 hθ2

/-- Last step of [SL ref 4, Prop 2.5] (FULLPROOF §3c): the Lamzouri-type bound Q ≥ 2tN − (2t − 1)S − t²d for
  every t ≥ 1, with d ≤ N − S and 2d ≤ N − O, gives the parity product (Q − S)(N − O) ≥ 2(N − S)². -/
theorem parity_of_lamzouri {Q N S O d : ℝ} (hd : 0 ≤ d) (hNS : d ≤ N - S) (hNO : 2 * d ≤ N - O)
    (hL : ∀ t : ℝ, 1 ≤ t → 2 * t * N - (2 * t - 1) * S - t ^ 2 * d ≤ Q) :
    2 * (N - S) ^ 2 ≤ (Q - S) * (N - O) := by
  rcases hd.lt_or_eq with hd0 | hd0
  · set t := (N - S) / d with htdef
    have ht1 : 1 ≤ t := by rw [htdef, le_div_iff₀ hd0]; linarith
    have htd : t * d = N - S := by rw [htdef]; field_simp
    have hQ := hL t ht1
    have hQS : (N - S) * t ≤ Q - S := by nlinarith
    have hQS0 : 0 ≤ Q - S := le_trans (by nlinarith) hQS
    calc 2 * (N - S) ^ 2 = ((N - S) * t) * (2 * d) := by rw [← htd]; ring
      _ ≤ (Q - S) * (N - O) := mul_le_mul hQS hNO (by linarith) hQS0
  · rw [← hd0] at hNS hNO hL
    have hQ1 := hL 1 le_rfl
    rcases (sub_nonneg.2 (show S ≤ N by linarith)).lt_or_eq with hpos | heq
    · exfalso
      have hL2 := hL (1 + |Q - S| / (N - S)) (by have := div_nonneg (abs_nonneg (Q - S)) hpos.le; linarith)
      have e : (1 + |Q - S| / (N - S)) * (N - S) = (N - S) + |Q - S| := by field_simp
      have ha := le_abs_self (Q - S)
      nlinarith [e, abs_nonneg (Q - S)]
    · rw [← heq]
      have : 0 ≤ Q - S := by linarith
      nlinarith

end E1Lean
