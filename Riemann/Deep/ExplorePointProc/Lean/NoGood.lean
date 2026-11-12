import Riemann.Deep.ExplorePointProc.Lean.Model

/-!
# Theorem 1 (DERIVATION §4), part 1: real-arithmetic helpers and the contradiction at one height T

no_good_height: at a height T with #E_λ(T) < δN(T), Lamzouri's lower bound for the weighted pair sum over
Z_0(T) contradicts step (iii') (FinSum.pairSum_weighted_le) once T is large. (Split from Theorem1.lean for the
290 MB memory cap.)
-/

namespace PPLean

open Finset

/-! ## Real-arithmetic helpers -/

/-- Weight step: v w(u) ≤ v + |M|/(4s²) for Montgomery's weight w(u) = 4s²/(4s² + u²) when |v| u² ≤ M. -/
lemma weight_le {s M v u : ℝ} (hs : 0 < s) (hM : |v| * u ^ 2 ≤ M) :
    v * (4 * s ^ 2 / (4 * s ^ 2 + u ^ 2)) ≤ v + |M| / (4 * s ^ 2) := by
  have h4s : 0 < 4 * s ^ 2 := by positivity
  have hden : 0 < 4 * s ^ 2 + u ^ 2 := by positivity
  have hη4 : |M| / (4 * s ^ 2) * (4 * s ^ 2) = |M| := div_mul_cancel₀ _ h4s.ne'
  have hη0 : 0 ≤ |M| / (4 * s ^ 2) := by positivity
  generalize |M| / (4 * s ^ 2) = η at hη4 hη0 ⊢
  have key : -(v * u ^ 2) ≤ |M| := by
    have h2 := neg_abs_le v
    have h3 := mul_nonneg (by linarith : 0 ≤ v + |v|) (sq_nonneg u)
    linarith [le_abs_self M]
  have h5 := mul_nonneg hη0 (sq_nonneg u)
  rw [mul_div_assoc', div_le_iff₀ hden]
  linarith

/-- 4/(4 + d²) = w(d s) with w(u) = 4s²/(4s² + u²): Lamzouri's weight at γ − γ' = d in unfolded units. -/
lemma weight_unfold {d s : ℝ} (hs : 0 < s) :
    (4 : ℝ) / (4 + d ^ 2) = 4 * s ^ 2 / (4 * s ^ 2 + (d * s) ^ 2) := by
  have h1 : (4 : ℝ) + d ^ 2 ≠ 0 := by positivity
  have h2 : (4 : ℝ) * s ^ 2 + (d * s) ^ 2 ≠ 0 := by positivity
  rw [div_eq_div_iff h1 h2]
  ring

/-- The boundary pair term: v w ≤ |v| ≤ |M|/λ² when |v| g² ≤ M, g > λ > 0, 0 ≤ w ≤ 1. -/
lemma last_term {lam g v w M : ℝ} (hlam : 0 < lam) (hg : lam < g) (hM : |v| * g ^ 2 ≤ M)
    (hw0 : 0 ≤ w) (hw1 : w ≤ 1) : v * w ≤ |M| / lam ^ 2 := by
  have hsq : lam ^ 2 ≤ g ^ 2 := by
    rw [sq, sq]
    exact mul_le_mul hg.le hg.le hlam.le (by linarith)
  have hab : |v| * lam ^ 2 ≤ |M| := by
    have h2 := mul_le_mul_of_nonneg_left hsq (abs_nonneg v)
    linarith [le_abs_self M]
  have hab' : |v| ≤ |M| / lam ^ 2 := by
    rw [le_div_iff₀ (by positivity)]
    exact hab
  have h1 := mul_le_mul_of_nonneg_right (le_abs_self v) hw0
  have h2 := mul_le_mul_of_nonneg_left hw1 (abs_nonneg v)
  linarith

lemma eta_small {s M μ : ℝ} (hμ : 0 < μ) (hs1 : 6 * |M| / μ + 1 ≤ s) :
    4 * (|M| / (4 * s ^ 2)) ≤ μ / 6 := by
  have hs0 : 0 ≤ 6 * |M| / μ := by positivity
  have hspos : 0 < s := by linarith
  have hs6 : 6 * |M| ≤ μ * s := by
    have h' : 6 * |M| / μ ≤ s := by linarith
    rw [div_le_iff₀ hμ] at h'
    linarith
  have hss : s ≤ s ^ 2 := by
    have := mul_le_mul_of_nonneg_left (by linarith : (1 : ℝ) ≤ s) hspos.le
    rw [sq]
    linarith
  have hs6' : 6 * |M| ≤ μ * s ^ 2 := by
    have := mul_le_mul_of_nonneg_left hss hμ.le
    linarith
  have h4s : 0 < 4 * s ^ 2 := by positivity
  rw [mul_div_assoc', div_le_div_iff₀ h4s (by norm_num)]
  linarith

lemma eps_small {β μ : ℝ} (hβ : 0 ≤ β) (hμ : 0 < μ) : 2 * β * (μ / (12 * β + 1)) ≤ μ / 6 := by
  rw [mul_div_assoc', div_le_div_iff₀ (by positivity) (by norm_num)]
  linarith

lemma exists_delta {δ₀ μ X : ℝ} (hδ₀ : 0 < δ₀) (hμ : 0 < μ) (hX : 0 ≤ X) :
    ∃ δ : ℝ, 0 < δ ∧ δ < δ₀ ∧ δ ≤ 1 / 2 ∧ δ * X ≤ μ / 6 := by
  have hpos : 0 < 6 * (X + 1) := by positivity
  refine ⟨min (δ₀ / 2) (min (1 / 2) (μ / (6 * (X + 1)))),
    lt_min (by positivity) (lt_min (by norm_num) (by positivity)),
    lt_of_le_of_lt (min_le_left _ _) (by linarith),
    (min_le_right _ _).trans (min_le_left _ _), ?_⟩
  have hd : min (δ₀ / 2) (min (1 / 2) (μ / (6 * (X + 1)))) ≤ μ / (6 * (X + 1)) :=
    (min_le_right _ _).trans (min_le_right _ _)
  calc min (δ₀ / 2) (min (1 / 2) (μ / (6 * (X + 1)))) * X ≤ μ / (6 * (X + 1)) * X :=
        mul_le_mul_of_nonneg_right hd hX
    _ ≤ μ / 6 := by
        rw [div_mul_eq_mul_div, div_le_div_iff₀ hpos (by norm_num)]
        linarith

/-! ## One height T: the contradiction -/

/-- At a height T with #E_λ(T) < δN(T), Lamzouri's lower bound (J − μ/6)N(T) for the weighted pair sum over
Z_0(T) contradicts the (iii') upper bound once T is large (DERIVATION §4). -/
lemma no_good_height {Z : Zeros} {C : TestCalc} (hZ : Z.Ordered)
    {θ lam U β c μ δ M KA K₀ ε₂ T : ℝ} {F A : ℝ → ℝ}
    (hlam0 : 0 < lam) (hU3 : U ≤ 3 * lam) (hFhe : ∀ u, C.hat F (-u) = C.hat F u)
    (hFneg : ∀ u, U ≤ u → C.hat F u ≤ 0) (hM : ∀ u, |C.hat F u| * u ^ 2 ≤ M)
    (hKA : ∀ y, lam ≤ y → |A y| ≤ KA) (hβ : 0 ≤ β) (hH : HH lam β c (C.hat F) A)
    (hμ : 0 < μ) (hμdef : μ = C.J θ F - (C.hat F 0 + 2 * c))
    (hδpos : 0 < δ) (hδ2 : δ ≤ 1 / 2) (hδ3 : δ * (|C.hat F 0 + 2 * c| + 2 * β) ≤ μ / 6)
    (hs1 : 6 * |M| / μ + 1 ≤ Z.sc T) (hε₂ : 2 * β * ε₂ ≤ μ / 6)
    (hRvT : T * Z.sc T ≤ (1 + ε₂) * (Z.N T : ℝ))
    (hK₀ : K₀ = 4 * |c| + 4 * β + 4 * |KA| + 2 * (|M| / lam ^ 2))
    (hKN : K₀ ≤ μ * (Z.N T : ℝ) / 6) (hN4 : 4 ≤ (Z.N T : ℝ))
    (hcon : ((Elam Z lam T).card : ℝ) < δ * (Z.N T : ℝ))
    (hlow : (C.J θ F - μ / 6) * (Z.N T : ℝ) ≤ ∑ a ∈ Z0 Z lam T, ∑ b ∈ Z0 Z lam T,
      C.hat F ((Z.γ a - Z.γ b) * Z.sc T) * (4 / (4 + (Z.γ a - Z.γ b) ^ 2))) : False := by
  have hs0 : 0 ≤ 6 * |M| / μ := by positivity
  have hspos : 0 < Z.sc T := by linarith only [hs0, hs1]
  have hNpos : 0 < (Z.N T : ℝ) := by linarith only [hN4]
  -- Montgomery's weight in unfolded units
  obtain ⟨w, hw⟩ : ∃ w : ℝ → ℝ, ∀ u, w u = 4 * Z.sc T ^ 2 / (4 * Z.sc T ^ 2 + u ^ 2) :=
    ⟨fun u => 4 * Z.sc T ^ 2 / (4 * Z.sc T ^ 2 + u ^ 2), fun u => rfl⟩
  have h4s : 0 < 4 * Z.sc T ^ 2 := by positivity
  have hweven : ∀ u, w (-u) = w u := fun u => by rw [hw, hw, neg_sq]
  have hwnn : ∀ u, 0 ≤ w u := fun u => by rw [hw]; positivity
  have hw0 : w 0 = 1 := by
    rw [hw, show (0 : ℝ) ^ 2 = 0 by norm_num, add_zero, div_self h4s.ne']
  have hw1 : ∀ u, w u ≤ 1 := fun u => by
    rw [hw, div_le_one (by positivity)]
    linarith only [sq_nonneg u]
  obtain ⟨η, hηdef⟩ : ∃ η : ℝ, η = |M| / (4 * Z.sc T ^ 2) := ⟨_, rfl⟩
  have hη0 : 0 ≤ η := by rw [hηdef]; positivity
  have hη : ∀ u, C.hat F u * w u ≤ C.hat F u + η := fun u => by
    rw [hw, hηdef]
    exact weight_le hspos (hM u)
  have hη6 : 4 * η ≤ μ / 6 := by rw [hηdef]; exact eta_small hμ hs1
  -- the pair sum over Z_0 as the pair sum of the increasing sequence xs
  have hPS : ∑ a ∈ Z0 Z lam T, ∑ b ∈ Z0 Z lam T,
      C.hat F ((Z.γ a - Z.γ b) * Z.sc T) * (4 / (4 + (Z.γ a - Z.γ b) ^ 2))
      = pairSum (fun u => C.hat F u * w u) (xs Z lam T) (Z0 Z lam T).card := by
    rw [← pairSum_Z0]
    refine Finset.sum_congr rfl (fun a _ => Finset.sum_congr rfl (fun b _ => ?_))
    show _ = C.hat F ((Z.γ a - Z.γ b) * Z.sc T) * w ((Z.γ a - Z.γ b) * Z.sc T)
    rw [hw, weight_unfold hspos]
  -- counts: #Z_0 = N − #E_λ ∈ [(1 − δ)N, N], ≥ 2
  have hcard : ((Z0 Z lam T).card : ℝ) + ((Elam Z lam T).card : ℝ) = (Z.N T : ℝ) := by
    exact_mod_cast card_Z0 Z lam T
  have hn1 : (1 - δ) * (Z.N T : ℝ) ≤ ((Z0 Z lam T).card : ℝ) := by linarith only [hcard, hcon]
  have hn2 : ((Z0 Z lam T).card : ℝ) ≤ (Z.N T : ℝ) := by
    have : (0 : ℝ) ≤ ((Elam Z lam T).card : ℝ) := Nat.cast_nonneg _
    linarith only [hcard, this]
  have hn2' : 2 ≤ (Z0 Z lam T).card := by
    have h0 : (1 / 2) * (Z.N T : ℝ) ≤ (1 - δ) * (Z.N T : ℝ) :=
      mul_le_mul_of_nonneg_right (by linarith only [hδ2]) hNpos.le
    have h1 : (2 : ℝ) ≤ ((Z0 Z lam T).card : ℝ) := by linarith only [h0, hn1, hN4]
    exact_mod_cast h1
  -- (iii'), weighted
  have hgap' : ∀ i, i + 1 < (Z0 Z lam T).card → lam ≤ xs Z lam T (i + 1) - xs Z lam T i :=
    fun i hi => (xs_gap hZ hspos.le hi).le
  have hbound := pairSum_weighted_le hlam0.le hU3 hFhe hFneg hH hweven hwnn hw0 hη hn2'
    (xs Z lam T) hgap'
  -- boundary terms
  have hg1 : lam < xs Z lam T 1 - xs Z lam T 0 := xs_gap hZ hspos.le (i := 0) (by omega)
  have hg2 : lam < xs Z lam T ((Z0 Z lam T).card - 1) - xs Z lam T ((Z0 Z lam T).card - 2) := by
    have h := xs_gap hZ hspos.le (lam := lam) (T := T) (i := (Z0 Z lam T).card - 2) (by omega)
    rwa [show (Z0 Z lam T).card - 2 + 1 = (Z0 Z lam T).card - 1 by omega] at h
  have hA1 : A (xs Z lam T ((Z0 Z lam T).card - 1) - xs Z lam T ((Z0 Z lam T).card - 2))
      - A (xs Z lam T 1 - xs Z lam T 0) ≤ 2 * |KA| := by
    have h1 := abs_le.1 (hKA _ hg2.le)
    have h2 := abs_le.1 (hKA _ hg1.le)
    have h3 := le_abs_self KA
    linarith only [h1.2, h2.1, h3]
  have hG := last_term
    (w := w (xs Z lam T ((Z0 Z lam T).card - 1) - xs Z lam T ((Z0 Z lam T).card - 2)))
    hlam0 hg2 (hM _) (hwnn _) (hw1 _)
  have hXb : xs Z lam T ((Z0 Z lam T).card - 2) - xs Z lam T 0 ≤ (1 + ε₂) * (Z.N T : ℝ) := by
    have h1 := (xs_bounds hZ hspos.le (lam := lam) (i := (Z0 Z lam T).card - 2) (by omega)).2
    have h2 := (xs_bounds hZ hspos.le (lam := lam) (i := 0) (by omega)).1
    linarith only [h1, h2, hRvT]
  have hup : pairSum (fun u => C.hat F u * w u) (xs Z lam T) (Z0 Z lam T).card ≤
      ((Z0 Z lam T).card : ℝ) * (C.hat F 0 + 2 * c) + 4 * ((Z0 Z lam T).card : ℝ) * η + 4 * |c|
      + 2 * β * (xs Z lam T ((Z0 Z lam T).card - 2) - xs Z lam T 0)
      - 2 * β * ((Z0 Z lam T).card : ℝ) + 4 * β + 4 * |KA| + 2 * (|M| / lam ^ 2) := by
    have hc := neg_abs_le c
    linarith only [hbound, hc, hη0, hA1, hG]
  rw [hPS] at hlow
  exact core_arith hNpos hμ hβ (by rw [hμdef]; ring) hn1 hn2 hδpos hlow hup hXb hη0 hδ3 hη6 hε₂
    (by rw [hK₀] at hKN; linarith only [hKN])

end PPLean
