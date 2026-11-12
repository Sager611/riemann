import Riemann.Deep.ExplorePointProc.Lean.NoGood

/-!
# Theorem 1 (DERIVATION §4): the hard-core augmentation of Lamzouri's small-gap criterion

theorem1_Cc2: for an even F ≥ 0 in C_c² with F^ ≤ 0 on [U, ∞), U ≤ 3λ, A bounded on [λ, ∞), β ≥ 0 and (H),
J_θ(F) > F^(0) + 2c gives D(λ) > 0. theorem1: the same for Schwartz f through Lemma T'. thmMain2: Lamzouri's
Thm Main2 as the case A = β = c = 0, U = λ. Inputs: the named hypotheses of Model.lean.
-/

namespace PPLean

open Finset

/-! ## Theorem 1 -/

/-- **Theorem 1 for C_c² kernels** (DERIVATION §4): D(λ) > 0. -/
theorem theorem1_Cc2 {Z : Zeros} {C : TestCalc} (hZ : Z.Ordered) (hsc : ScaleToInf Z) (hRvM : RvM Z)
    (h12 : Steps12 Z C) (hFC : FourierCc2 C)
    {θ lam U β c : ℝ} {F A : ℝ → ℝ}
    (hθ1 : 1 / 2 < θ) (hθ2 : θ < 4 / 7) (hlam0 : 0 < lam) (hlam1 : lam ≤ 1) (hU3 : U ≤ 3 * lam)
    (hFc : C.Cc2 F) (hFe : ∀ x, F (-x) = F x) (hFn : ∀ x, 0 ≤ F x)
    (hFneg : ∀ u, U ≤ u → C.hat F u ≤ 0)
    (hA : ∃ K : ℝ, ∀ y, lam ≤ y → |A y| ≤ K) (hβ : 0 ≤ β) (hH : HH lam β c (C.hat F) A)
    (hJ : C.hat F 0 + 2 * c < C.J θ F) : DPos Z lam := by
  obtain ⟨hFhe, M, hM⟩ := hFC F hFc hFe
  obtain ⟨KA, hKA⟩ := hA
  obtain ⟨μ, hμdef⟩ : ∃ μ : ℝ, μ = C.J θ F - (C.hat F 0 + 2 * c) := ⟨_, rfl⟩
  have hμ : 0 < μ := by rw [hμdef]; linarith
  obtain ⟨δ₀, hδ₀, hδ⟩ := h12 θ lam hθ1 hθ2 hlam0 hlam1 F hFc hFe hFn (μ / 6) (by positivity)
  obtain ⟨δ, hδpos, hδ1, hδ2, hδ3⟩ := exists_delta hδ₀ hμ
    (by positivity : (0 : ℝ) ≤ |C.hat F 0 + 2 * c| + 2 * β)
  obtain ⟨T₀, hT₀⟩ := hδ δ hδpos hδ1
  obtain ⟨T₁, hT₁⟩ := hsc (6 * |M| / μ + 1)
  obtain ⟨T₂, hT₂⟩ := hRvM (μ / (12 * β + 1)) (by positivity)
  obtain ⟨T₃, hT₃⟩ := hRvM 1 one_pos
  obtain ⟨K₀, hK₀⟩ : ∃ K₀ : ℝ, K₀ = 4 * |c| + 4 * β + 4 * |KA| + 2 * (|M| / lam ^ 2) := ⟨_, rfl⟩
  have hK₀nn : 0 ≤ K₀ := by rw [hK₀]; positivity
  refine ⟨δ / 2, by positivity, max (max T₀ T₁) (max (max T₂ T₃) (2 * (6 * K₀ / μ + 4))), ?_⟩
  intro T hT
  obtain ⟨⟨hT0, hT1⟩, ⟨hT2, hT3⟩, hT4⟩ :
      (T₀ ≤ T ∧ T₁ ≤ T) ∧ (T₂ ≤ T ∧ T₃ ≤ T) ∧ 2 * (6 * K₀ / μ + 4) ≤ T := by
    simpa only [max_le_iff] using hT
  have hs1 := hT₁ T hT1
  have hs0 : 0 ≤ 6 * |M| / μ := by positivity
  have hs_ge1 : 1 ≤ Z.sc T := by linarith
  have hK0μ : 0 ≤ 6 * K₀ / μ := by positivity
  have hTpos : 0 < T := by linarith
  have hNT : T ≤ 2 * (Z.N T : ℝ) := by
    have h := hT₃ T hT3
    have h' := mul_le_mul_of_nonneg_left hs_ge1 hTpos.le
    linarith
  have hNbig : 6 * K₀ / μ + 4 ≤ (Z.N T : ℝ) := by linarith
  have hKN : K₀ ≤ μ * (Z.N T : ℝ) / 6 := by
    have h1 : 6 * K₀ / μ ≤ (Z.N T : ℝ) := by linarith
    rw [div_le_iff₀ hμ] at h1
    linarith
  have hE : δ * (Z.N T : ℝ) ≤ ((Elam Z lam T).card : ℝ) := by
    by_contra hcon
    have hcon' := not_le.1 hcon
    exact no_good_height hZ hlam0 hU3 hFhe hFneg hM hKA hβ hH hμ hμdef hδpos hδ2 hδ3 hs1
      (eps_small hβ hμ) (hT₂ T hT2) hK₀ hKN (by linarith) hcon' (hT₀ T hT0 hcon')
  have h2 := card_Elam_le Z lam T
  have h2' : ((Elam Z lam T).card : ℝ) ≤ 2 * ((smallGaps Z lam T).card : ℝ) := by exact_mod_cast h2
  linarith

/-- **Theorem 1** (DERIVATION §4): Schwartz f, through Lemma T'. -/
theorem theorem1 {Z : Zeros} {C : TestCalc} (hZ : Z.Ordered) (hsc : ScaleToInf Z) (hRvM : RvM Z)
    (h12 : Steps12 Z C) (hFC : FourierCc2 C) (hT : LemmaT' C)
    {θ lam U β c : ℝ} {f A : ℝ → ℝ}
    (hθ1 : 1 / 2 < θ) (hθ2 : θ < 4 / 7) (hlam0 : 0 < lam) (hlam1 : lam ≤ 1)
    (hU1 : lam < U) (hU3 : U ≤ 3 * lam)
    (hfS : C.Sch f) (hfe : ∀ x, f (-x) = f x) (hfn : ∀ x, 0 ≤ f x)
    (hfneg : ∀ u, U ≤ u → C.hat f u ≤ 0)
    (hA : ∃ K : ℝ, ∀ y, lam ≤ y → |A y| ≤ K) (hβ : 0 ≤ β) (hH : HH lam β c (C.hat f) A)
    (hJ : C.hat f 0 + 2 * c < C.J θ f) : DPos Z lam := by
  obtain ⟨F, hFc, hFe, hFn, hFneg, hFH, hFJ⟩ :=
    hT θ lam U β c f A hlam0 hU1.le (by linarith) hfS hfe hfn hfneg hH hJ
  exact theorem1_Cc2 hZ hsc hRvM h12 hFC hθ1 hθ2 hlam0 hlam1 hU3 hFc hFe hFn hFneg hA hβ hFH hFJ

/-- Lamzouri's Thm Main2 (TeX l.392–400) as Theorem 1's case A = β = c = 0, U = λ (sign condition on the
closed [λ, ∞), which Lamzouri's |u| > λ gives by continuity of f^; no normalisation f^(0) = 1 needed). -/
theorem thmMain2 {Z : Zeros} {C : TestCalc} (hZ : Z.Ordered) (hsc : ScaleToInf Z) (hRvM : RvM Z)
    (h12 : Steps12 Z C) (hFC : FourierCc2 C) {θ lam : ℝ} {f : ℝ → ℝ}
    (hθ1 : 1 / 2 < θ) (hθ2 : θ < 4 / 7) (hlam0 : 0 < lam) (hlam1 : lam ≤ 1)
    (hfc : C.Cc2 f) (hfe : ∀ x, f (-x) = f x) (hfn : ∀ x, 0 ≤ f x)
    (hfneg : ∀ u, lam ≤ u → C.hat f u ≤ 0) (hJ : C.hat f 0 < C.J θ f) : DPos Z lam := by
  refine theorem1_Cc2 (U := lam) (β := 0) (c := 0) (A := fun _ => 0) hZ hsc hRvM h12 hFC hθ1 hθ2
    hlam0 hlam1 (by linarith) hfc hfe hfn hfneg ⟨0, fun y _ => by simp⟩ le_rfl ?_ (by linarith)
  intro x y hx hy
  have h1 := hfneg x hx
  have h2 := hfneg (x + y) (by linarith)
  show C.hat f x + C.hat f (x + y) + 0 - 0 + 0 * (1 - x) ≤ 0
  linarith

/-- Homogeneity (DERIVATION §4, Normalisation): dividing the kernel, A, β and c by t > 0 (t = f^(0)) keeps (H). -/
lemma HH_div {lam β c t : ℝ} {F A : ℝ → ℝ} (ht : 0 < t) (h : HH lam β c F A) :
    HH lam (β / t) (c / t) (fun u => F u / t) (fun u => A u / t) := by
  intro x y hx hy
  have h1 := h x y hx hy
  have e : F x / t + F (x + y) / t + A x / t - A y / t + β / t * (1 - x)
      = (F x + F (x + y) + A x - A y + β * (1 - x)) / t := by ring
  show F x / t + F (x + y) / t + A x / t - A y / t + β / t * (1 - x) ≤ c / t
  rw [e]
  exact div_le_div_of_nonneg_right h1 ht.le

/-- ... and the margin J > F^(0) + 2c. -/
lemma margin_div {F0 c J t : ℝ} (ht : 0 < t) (h : F0 + 2 * c < J) :
    F0 / t + 2 * (c / t) < J / t := by
  have e : F0 / t + 2 * (c / t) = (F0 + 2 * c) / t := by ring
  rw [e]
  exact div_lt_div_of_pos_right h ht

end PPLean
