import Riemann.Deep.ExplorePointProc.Lean.FinSum
import Mathlib.Data.Finset.Sort
import Mathlib.Order.Interval.Finset.Nat
import Mathlib.Tactic.Positivity

/-!
# Theorem 1 of explore40 (DERIVATION §4) from named analytic hypotheses

Model of the zeros (ordinates only, Lamzouri 2611.14486 l.175–183, 404–410): γ n (n ≥ 1, with multiplicity,
nondecreasing), N T = N(T), and the unfolding scale sc T = log T/(2π) (abstract: Real.log and π are not
importable under the 290 MB cap). D(λ) > 0 is DPos, E_λ(T) is Elam, Z_0(T) = Z(T) ∖ E_λ(T) is Z0; the
hard-core property of Z_0, the count #Z_0 = N(T) − #E_λ(T), #E_λ ≤ 2·#(small gaps) and the reindexing of the
pair sum over Z_0 as a sum over the increasing sequence xs of its unfolded ordinates are proved here.
The analytic inputs are the named hypotheses Steps12 (Lamzouri's (Eq:ReductionGoodSet2), l.425–440, in the
on-line form of l.444–446), FourierCc2, LemmaT', ScaleToInf, RvM; the test-function calculus (Fourier
transform, J_θ, the classes C_c² and Schwartz) is abstract (TestCalc).
-/

namespace PPLean

open Finset

/-- Ordinates of the nontrivial zeros of ζ counted with multiplicity (γ n for n ≥ 1; γ 0 unused),
N T = N(T) (zeros with 0 < γ ≤ T) and the unfolding scale sc T = log T/(2π). -/
structure Zeros where
  γ : ℕ → ℝ
  N : ℝ → ℕ
  sc : ℝ → ℝ

/-- Structural facts of the model: γ_n nondecreasing for n ≥ 1, γ_1 > 0, N(T) = #{n ≥ 1 : γ_n ≤ T}. -/
structure Zeros.Ordered (Z : Zeros) : Prop where
  mono : ∀ k : ℕ, 1 ≤ k → Z.γ k ≤ Z.γ (k + 1)
  pos : 0 < Z.γ 1
  count : ∀ (T : ℝ) (k : ℕ), 1 ≤ k → (Z.γ k ≤ T ↔ k ≤ Z.N T)

lemma Zeros.Ordered.mono' {Z : Zeros} (hZ : Z.Ordered) {a b : ℕ} (ha : 1 ≤ a) (hab : a ≤ b) :
    Z.γ a ≤ Z.γ b := by
  induction b, hab using Nat.le_induction with
  | base => exact le_rfl
  | succ b hab ih => exact ih.trans (hZ.mono b (by omega))

/-- Lamzouri l.183: {n : 0 < γ_n ≤ T, γ_{n+1} − γ_n ≤ 2πλ/log T}, unfolded: (γ_{n+1} − γ_n)·log T/(2π) ≤ λ. -/
noncomputable def smallGaps (Z : Zeros) (lam T : ℝ) : Finset ℕ :=
  (Icc 1 (Z.N T)).filter (fun k => (Z.γ (k + 1) - Z.γ k) * Z.sc T ≤ lam)

/-- D(λ) > 0, where D(λ) = liminf_{T→∞} #smallGaps(λ, T)/N(T) (l.183): the ratio is eventually ≥ some δ > 0. -/
def DPos (Z : Zeros) (lam : ℝ) : Prop :=
  ∃ δ : ℝ, 0 < δ ∧ ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → δ * (Z.N T : ℝ) ≤ ((smallGaps Z lam T).card : ℝ)

/-- E_λ(T), Lamzouri (Eq:DefMultiSetE) l.404: ρ_n with 0 < γ_n ≤ T, n ≥ 2 and
min(γ_{n+1} − γ_n, γ_n − γ_{n−1}) ≤ 2πλ/log T. -/
noncomputable def Elam (Z : Zeros) (lam T : ℝ) : Finset ℕ :=
  (Icc 1 (Z.N T)).filter
    (fun k => 2 ≤ k ∧ min (Z.γ (k + 1) - Z.γ k) (Z.γ k - Z.γ (k - 1)) * Z.sc T ≤ lam)

/-- Z_0(T) = Z(T) ∖ E_λ(T), l.410. -/
noncomputable def Z0 (Z : Zeros) (lam T : ℝ) : Finset ℕ := Icc 1 (Z.N T) \ Elam Z lam T

/-- l.406–408: #E_λ(T) ≤ 2 #{small gaps} (each ρ_n ∈ E_λ has a small gap at index n or n − 1). -/
lemma card_Elam_le (Z : Zeros) (lam T : ℝ) :
    (Elam Z lam T).card ≤ 2 * (smallGaps Z lam T).card := by
  have hsub : Elam Z lam T ⊆ smallGaps Z lam T ∪ (smallGaps Z lam T).image (· + 1) := by
    intro k hk
    simp only [Elam, mem_filter, mem_Icc] at hk
    obtain ⟨⟨hk1, hkN⟩, hk2, hmin⟩ := hk
    rcases min_choice (Z.γ (k + 1) - Z.γ k) (Z.γ k - Z.γ (k - 1)) with h | h
    · rw [h] at hmin
      apply mem_union_left
      simp only [smallGaps, mem_filter, mem_Icc]
      exact ⟨⟨hk1, hkN⟩, hmin⟩
    · rw [h] at hmin
      apply mem_union_right
      rw [mem_image]
      refine ⟨k - 1, ?_, by show k - 1 + 1 = k; omega⟩
      simp only [smallGaps, mem_filter, mem_Icc]
      refine ⟨⟨by omega, by omega⟩, ?_⟩
      rw [show k - 1 + 1 = k by omega]
      exact hmin
  calc (Elam Z lam T).card
      ≤ (smallGaps Z lam T ∪ (smallGaps Z lam T).image (· + 1)).card := card_le_card hsub
    _ ≤ (smallGaps Z lam T).card + ((smallGaps Z lam T).image (· + 1)).card := card_union_le _ _
    _ ≤ (smallGaps Z lam T).card + (smallGaps Z lam T).card :=
        Nat.add_le_add_left card_image_le _
    _ = 2 * (smallGaps Z lam T).card := by ring

lemma smallGaps_mono (Z : Zeros) {lam lam' : ℝ} (h : lam ≤ lam') (T : ℝ) :
    smallGaps Z lam T ⊆ smallGaps Z lam' T := by
  intro k hk
  simp only [smallGaps, mem_filter] at hk ⊢
  exact ⟨hk.1, hk.2.trans h⟩

/-- D is monotone: D(λ) > 0 and λ ≤ λ' give D(λ') > 0. -/
theorem DPos.mono {Z : Zeros} {lam lam' : ℝ} (h : lam ≤ lam') (hD : DPos Z lam) : DPos Z lam' := by
  obtain ⟨δ, hδ, T₀, hT⟩ := hD
  refine ⟨δ, hδ, T₀, fun T hT' => (hT T hT').trans ?_⟩
  exact_mod_cast card_le_card (smallGaps_mono Z h T)

lemma card_Z0 (Z : Zeros) (lam T : ℝ) :
    (Z0 Z lam T).card + (Elam Z lam T).card = Z.N T := by
  have hsub : Elam Z lam T ⊆ Icc 1 (Z.N T) := filter_subset _ _
  have h := card_sdiff_add_card_eq_card hsub
  rw [Nat.card_Icc] at h
  unfold Z0
  omega

/-- l.441–443: distinct elements of Z_0 are more than λ apart (unfolded). -/
lemma Z0_sep {Z : Zeros} (hZ : Z.Ordered) {lam T : ℝ} (hsc : 0 ≤ Z.sc T) {a b : ℕ}
    (ha : a ∈ Z0 Z lam T) (hb : b ∈ Z0 Z lam T) (hab : a < b) :
    lam < (Z.γ b - Z.γ a) * Z.sc T := by
  obtain ⟨haI, -⟩ := mem_sdiff.1 ha
  obtain ⟨hbI, hbE⟩ := mem_sdiff.1 hb
  rw [mem_Icc] at haI
  have hnot : ¬ (2 ≤ b ∧ min (Z.γ (b + 1) - Z.γ b) (Z.γ b - Z.γ (b - 1)) * Z.sc T ≤ lam) :=
    fun h => hbE (mem_filter.2 ⟨hbI, h⟩)
  have h1 : lam < min (Z.γ (b + 1) - Z.γ b) (Z.γ b - Z.γ (b - 1)) * Z.sc T := by
    by_contra h
    exact hnot ⟨by omega, not_lt.1 h⟩
  have h2 : min (Z.γ (b + 1) - Z.γ b) (Z.γ b - Z.γ (b - 1)) ≤ Z.γ b - Z.γ (b - 1) :=
    min_le_right _ _
  have h3 : Z.γ a ≤ Z.γ (b - 1) := hZ.mono' haI.1 (by omega)
  have h4 := mul_le_mul_of_nonneg_right h2 hsc
  have h5 : (Z.γ b - Z.γ (b - 1)) * Z.sc T ≤ (Z.γ b - Z.γ a) * Z.sc T :=
    mul_le_mul_of_nonneg_right (by linarith) hsc
  linarith

/-- The unfolded ordinates of Z_0(T) in increasing order: x_i = γ_{e(i)} · log T/(2π), i < #Z_0. -/
noncomputable def xs (Z : Zeros) (lam T : ℝ) (i : ℕ) : ℝ :=
  if h : i < (Z0 Z lam T).card then Z.γ ((Z0 Z lam T).orderEmbOfFin rfl ⟨i, h⟩) * Z.sc T else 0

lemma sum_Z0_eq (Z : Zeros) (lam T : ℝ) (φ : ℕ → ℝ) :
    ∑ a ∈ Z0 Z lam T, φ a = ∑ i ∈ range (Z0 Z lam T).card,
      (if h : i < (Z0 Z lam T).card then φ ((Z0 Z lam T).orderEmbOfFin rfl ⟨i, h⟩) else 0) := by
  symm
  apply Finset.sum_bij (fun i hi => (Z0 Z lam T).orderEmbOfFin rfl ⟨i, mem_range.1 hi⟩)
  · intro i hi
    exact orderEmbOfFin_mem _ _ _
  · intro i₁ h₁ i₂ h₂ h
    have h' := ((Z0 Z lam T).orderEmbOfFin rfl).injective h
    exact congrArg Fin.val h'
  · intro b hb
    have hb' : b ∈ Finset.map ((Z0 Z lam T).orderEmbOfFin rfl).toEmbedding univ := by
      rw [map_orderEmbOfFin_univ]
      exact hb
    obtain ⟨i, _, hi⟩ := mem_map.1 hb'
    exact ⟨i.1, mem_range.2 i.2, hi⟩
  · intro i hi
    rw [dite_eq_left (mem_range.1 hi)]

/-- Lamzouri's pair sum over Z_0, as the pair sum of the increasing sequence xs. -/
lemma pairSum_Z0 (Z : Zeros) (lam T : ℝ) (Φ : ℝ → ℝ) :
    ∑ a ∈ Z0 Z lam T, ∑ b ∈ Z0 Z lam T, Φ ((Z.γ a - Z.γ b) * Z.sc T)
      = pairSum Φ (xs Z lam T) (Z0 Z lam T).card := by
  unfold pairSum
  rw [sum_Z0_eq Z lam T (fun a => ∑ b ∈ Z0 Z lam T, Φ ((Z.γ a - Z.γ b) * Z.sc T))]
  apply Finset.sum_congr rfl
  intro i hi
  rw [dite_eq_left (mem_range.1 hi), sum_Z0_eq Z lam T (fun b => Φ ((Z.γ _ - Z.γ b) * Z.sc T))]
  apply Finset.sum_congr rfl
  intro j hj
  rw [dite_eq_left (mem_range.1 hj)]
  simp only [xs, dite_eq_left (mem_range.1 hi), dite_eq_left (mem_range.1 hj)]
  congr 1
  ring

lemma xs_gap {Z : Zeros} (hZ : Z.Ordered) {lam T : ℝ} (hsc : 0 ≤ Z.sc T) {i : ℕ}
    (hi : i + 1 < (Z0 Z lam T).card) : lam < xs Z lam T (i + 1) - xs Z lam T i := by
  have hi' : i < (Z0 Z lam T).card := by omega
  simp only [xs, dite_eq_left hi, dite_eq_left hi']
  rw [← sub_mul]
  apply Z0_sep hZ hsc (orderEmbOfFin_mem _ _ _) (orderEmbOfFin_mem _ _ _)
  apply ((Z0 Z lam T).orderEmbOfFin rfl).strictMono
  exact Fin.mk_lt_mk.2 (Nat.lt_succ_self i)

lemma xs_bounds {Z : Zeros} (hZ : Z.Ordered) {lam T : ℝ} (hsc : 0 ≤ Z.sc T) {i : ℕ}
    (hi : i < (Z0 Z lam T).card) : 0 ≤ xs Z lam T i ∧ xs Z lam T i ≤ T * Z.sc T := by
  simp only [xs, dite_eq_left hi]
  have hm := orderEmbOfFin_mem (Z0 Z lam T) rfl ⟨i, hi⟩
  have hm' := (mem_sdiff.1 hm).1
  rw [mem_Icc] at hm'
  have h1 := hZ.mono' le_rfl hm'.1
  have h2 := (hZ.count T _ hm'.1).2 hm'.2
  have h0 := hZ.pos
  constructor
  · exact mul_nonneg (by linarith) hsc
  · exact mul_le_mul_of_nonneg_right h2 hsc

/-! ## The abstract test-function calculus and the named analytic hypotheses -/

/-- Abstract calculus. Intended: hat f u = f^(u) = ∫ f(α) e^{−2πiαu} dα (real for even real f);
J θ f = f(0) + 2∫_0^θ f(α) α dα; Cc2 f: f ∈ C_c²(ℝ); Sch f: f Schwartz. -/
structure TestCalc where
  hat : (ℝ → ℝ) → ℝ → ℝ
  J : ℝ → (ℝ → ℝ) → ℝ
  Cc2 : (ℝ → ℝ) → Prop
  Sch : (ℝ → ℝ) → Prop

/-- Lamzouri's steps (i)–(ii): (Eq:ReductionGoodSet2), TeX l.425–440 (Prop SmallSetZeros, l.331–340, for
g = ψ_ε F, then the positivity identity l.430–436, F ≥ g ≥ 0), with the summand in its on-line form
F^((γ − γ')log T/2π)·4/(4 + (γ − γ')²) (l.444–446: Z_0 lies on the critical line), and the error
O(δ^{1/16−7θ/64}) for T large written as an ε–δ–T₀ statement. -/
def Steps12 (Z : Zeros) (C : TestCalc) : Prop :=
  ∀ θ lam : ℝ, 1 / 2 < θ → θ < 4 / 7 → 0 < lam → lam ≤ 1 →
  ∀ F : ℝ → ℝ, C.Cc2 F → (∀ x, F (-x) = F x) → (∀ x, 0 ≤ F x) →
  ∀ ε : ℝ, 0 < ε → ∃ δ₀ : ℝ, 0 < δ₀ ∧ ∀ δ : ℝ, 0 < δ → δ < δ₀ → ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T →
    ((Elam Z lam T).card : ℝ) < δ * (Z.N T : ℝ) →
    (C.J θ F - ε) * (Z.N T : ℝ) ≤ ∑ a ∈ Z0 Z lam T, ∑ b ∈ Z0 Z lam T,
      C.hat F ((Z.γ a - Z.γ b) * Z.sc T) * (4 / (4 + (Z.γ a - Z.γ b) ^ 2))

/-- Fourier facts for an even F ∈ C_c²: F^ is even and F^(u) = O(u⁻²) (|F^(u)| ≤ ‖F''‖₁/(4π²u²)). -/
def FourierCc2 (C : TestCalc) : Prop :=
  ∀ F : ℝ → ℝ, C.Cc2 F → (∀ x, F (-x) = F x) →
    (∀ u, C.hat F (-u) = C.hat F u) ∧ ∃ M : ℝ, ∀ u, |C.hat F u| * u ^ 2 ≤ M

/-- Lemma T' (DERIVATION §4): F = f h_R + c_r G (C_c², even, ≥ 0) with F^ ≤ f^ on [λ, ∞), J_θ(F) ≥ J_θ(f),
F^(0) ≤ f^(0) + 2c_r: so F^ ≤ 0 on [U, ∞), (H) holds for F with the same A, β, c, and the margin stays > 0. -/
def LemmaT' (C : TestCalc) : Prop :=
  ∀ (θ lam U β c : ℝ) (f A : ℝ → ℝ), 0 < lam → lam ≤ U → 0 ≤ θ →
    C.Sch f → (∀ x, f (-x) = f x) → (∀ x, 0 ≤ f x) → (∀ u, U ≤ u → C.hat f u ≤ 0) →
    HH lam β c (C.hat f) A → C.hat f 0 + 2 * c < C.J θ f →
    ∃ F : ℝ → ℝ, C.Cc2 F ∧ (∀ x, F (-x) = F x) ∧ (∀ x, 0 ≤ F x) ∧
      (∀ u, U ≤ u → C.hat F u ≤ 0) ∧ HH lam β c (C.hat F) A ∧ C.hat F 0 + 2 * c < C.J θ F

/-- log T/(2π) → ∞. -/
def ScaleToInf (Z : Zeros) : Prop := ∀ K : ℝ, ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → K ≤ Z.sc T

/-- Riemann–von Mangoldt, N(T) = (T/2π) log(T/2πe) + O(log T): T log T/(2π) ≤ (1 + ε) N(T) for large T
(the unfolded ordinates of Z(T) lie in (0, T log T/2π]: DERIVATION's x_{N_0−1} − x_1 ≤ N(1 + o(1))). -/
def RvM (Z : Zeros) : Prop :=
  ∀ ε : ℝ, 0 < ε → ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → T * Z.sc T ≤ (1 + ε) * (Z.N T : ℝ)

/-! ## The contradiction (DERIVATION §4, end of the proof) -/

/-- The final bookkeeping: lower bound (J − μ/6)N against the (iii') upper bound. -/
lemma core_arith {N n a J μ δ η β ε₂ c X PS B1 B2 : ℝ}
    (hN : 0 < N) (hμ : 0 < μ) (hβ : 0 ≤ β) (hJ : J = a + μ)
    (hn1 : (1 - δ) * N ≤ n) (hn2 : n ≤ N) (hδ0 : 0 < δ)
    (hlow : (J - μ / 6) * N ≤ PS)
    (hup : PS ≤ n * a + 4 * n * η + 4 * |c| + 2 * β * X - 2 * β * n + 4 * β + B1 + B2)
    (hX : X ≤ (1 + ε₂) * N) (hη : 0 ≤ η)
    (hδ3 : δ * (|a| + 2 * β) ≤ μ / 6) (hη6 : 4 * η ≤ μ / 6) (hε₂ : 2 * β * ε₂ ≤ μ / 6)
    (hK : 4 * |c| + 4 * β + B1 + B2 ≤ μ * N / 6) : False := by
  subst hJ
  have h1 : n * a ≤ N * a + δ * N * |a| := by
    rcases le_or_gt 0 a with ha | ha
    · have h := mul_le_mul_of_nonneg_right hn2 ha
      have : 0 ≤ δ * N * |a| := by positivity
      linarith
    · have h := mul_le_mul_of_nonpos_right hn1 ha.le
      rw [abs_of_neg ha]
      linarith
  have h2 : n * η ≤ N * η := mul_le_mul_of_nonneg_right hn2 hη
  have h3 : β * X ≤ β * ((1 + ε₂) * N) := mul_le_mul_of_nonneg_left hX hβ
  have h4 : β * ((1 - δ) * N) ≤ β * n := mul_le_mul_of_nonneg_left hn1 hβ
  have h5 : δ * (|a| + 2 * β) * N ≤ μ / 6 * N := mul_le_mul_of_nonneg_right hδ3 hN.le
  have h6 : 4 * η * N ≤ μ / 6 * N := mul_le_mul_of_nonneg_right hη6 hN.le
  have h7 : 2 * β * ε₂ * N ≤ μ / 6 * N := mul_le_mul_of_nonneg_right hε₂ hN.le
  have h8 : 0 < μ * N := mul_pos hμ hN
  linarith

end PPLean
