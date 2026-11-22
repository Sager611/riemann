module

public import Riemann.Deep.LeanLi.LemmaA2

/-!
# LeanLi: Lemma A (proofs.md), the exact decomposition of `R(n)` at height `T`

Under RH, for `0 < T ≤ X`:
`R(n) = -2 W_T(n) + 2n a_{T,X}(n) + Rest`, `W_T(n) = Σ_{γ≤T} cos nθ_γ`,
`a_{T,X}(n) = ∫_T^X S(t)/(t² + 1/4) sin nθ(t) dt`, and `Rest` a sum of explicitly bounded terms.
-/

@[expose] public section

open Real Filter Set MeasureTheory

namespace LeanLi

/-- `g_n(t) = 2(1 - cos nθ(t))`. -/
noncomputable def gfun (n : ℕ) (t : ℝ) : ℝ := 2 * (1 - Real.cos (n * liTheta t))

/-- `g_n'(t)`. -/
noncomputable def gder (n : ℕ) (t : ℝ) : ℝ := 2 * (Real.sin (n * liTheta t) * (n * -(1 / (t ^ 2 + 1 / 4))))

lemma hasDerivAt_gfun (n : ℕ) (t : ℝ) : HasDerivAt (gfun n) (gder n t) t := by
  have := (((hasDerivAt_liTheta t).const_mul (n : ℝ)).cos.const_sub 1).const_mul 2
  unfold gfun gder
  convert this using 1
  ring

lemma continuous_liTheta : Continuous liTheta :=
  continuous_iff_continuousAt.2 fun t => (hasDerivAt_liTheta t).continuousAt

lemma continuous_gder (n : ℕ) : Continuous (gder n) := by
  unfold gder
  have := continuous_inv_sq_add
  have := continuous_liTheta
  fun_prop

/-- `a_{T,X}(n)`. -/
noncomputable def aCoef (T X : ℝ) (n : ℕ) : ℝ :=
  ∫ t in T..X, zetaS t / (t ^ 2 + 1 / 4) * Real.sin (n * liTheta t)

/-- `W_T(n) = Σ_{γ≤T} cos nθ_γ` (with multiplicity). -/
noncomputable def Wv (T : ℝ) (n : ℕ) : ℝ := ∑ k ∈ ZK T, Real.cos (n * thK k)

/-- The integrand of `I(n)` in `t`. -/
noncomputable def hI (n : ℕ) (t : ℝ) : ℝ :=
  1 / Real.pi * ((1 - Real.cos (n * liTheta t)) * Real.log (t / (2 * Real.pi)))

/-- The remainder of Lemma A. -/
noncomputable def Rest (n : ℕ) (T X : ℝ) : ℝ :=
  2 * Real.cos (n * liTheta T) * ((zeroCount T : ℝ) - Afun T) + 2 + 2 * (∫ t in (0 : ℝ)..T, rfun t) +
    (∫ t in T..X, gfun n t * rfun t) + (Iclosed n - liTrend n) - (∫ t in Ioi X, hI n t) +
    1 / Real.pi * (∫ t in (0 : ℝ)..T, Real.cos (n * liTheta t) * Real.log (t / (2 * Real.pi))) +
    gfun n X * ((zeroCount X : ℝ) - Afun X) + ∑' k : ↑((ZK X : Set ZeroIdx)ᶜ), gfun n (ordK k)

/-- **Lemma A** (exact decomposition). -/
theorem decomp (L : LiteratureInputs) (hRH : RiemannHypothesis) (n : ℕ) {T X : ℝ} (hT : 0 < T)
    (hTX : T ≤ X) :
    R n = -2 * Wv T n + 2 * n * aCoef T X n + Rest n T X := by
  have hπ := Real.pi_pos
  have hX : 0 < X := lt_of_lt_of_le hT hTX
  -- Li's formula
  have hlam : HasSum (fun k : ZeroIdx => gfun n (ordK k)) (lam n) := hasSum_lam L hRH n
  have hsum : Summable fun k : ZeroIdx => gfun n (ordK k) := hlam.summable
  have h1 : lam n = (∑ k ∈ ZK X, gfun n (ordK k)) + ∑' k : ↑((ZK X : Set ZeroIdx)ᶜ), gfun n (ordK k) := by
    rw [← hlam.tsum_eq, ← hsum.sum_add_tsum_compl (s := ZK X)]
  -- split the zeros at T
  have h2 : ∑ k ∈ ZK X, gfun n (ordK k) =
      (∑ k ∈ ZK T, gfun n (ordK k)) + ∑ ρ ∈ (zset X).filter (fun ρ => T < ρ.im), (mult ρ : ℝ) * gfun n ρ.im := by
    rw [sum_ZK X, sum_ZK T, ← Finset.sum_filter_add_sum_filter_not (zset X) (fun ρ => ρ.im ≤ T),
      zset_filter_le hTX]
    congr 1
    refine Finset.sum_congr ?_ fun _ _ => rfl
    ext ρ; simp only [Finset.mem_filter, not_le]
  have h3 : ∑ k ∈ ZK T, gfun n (ordK k) = 2 * (zeroCount T : ℝ) - 2 * Wv T n := by
    rw [← card_ZK, Wv]
    simp only [gfun, thK]
    rw [Finset.mul_sum, Finset.card_eq_sum_ones, Nat.cast_sum, Finset.mul_sum, ← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl fun k _ => ?_
    push_cast; ring
  -- Abel on (T, X]
  have habel := zsum_abel_AS hTX (f := gfun n) (f' := gder n) (fun t _ => hasDerivAt_gfun n t)
    (continuous_gder n).continuousOn
  -- ∫ g' S = -2n a
  have h4 : ∫ t in T..X, gder n t * zetaS t = -(2 * n * aCoef T X n) := by
    unfold aCoef
    rw [← intervalIntegral.integral_const_mul, ← intervalIntegral.integral_neg]
    refine intervalIntegral.integral_congr fun t _ => ?_
    simp only [gder]; ring
  -- ∫ g A' = ∫ hI + ∫ g r on [T, X]
  have hgc : Continuous (gfun n) := by unfold gfun; have := continuous_liTheta; fun_prop
  have hiT : ∀ a b, 0 ≤ a → 0 ≤ b → IntervalIntegrable (hI n) volume a b := by
    intro a b ha hb
    have := (intervalIntegrable_log_div_two_pi ha hb).continuousOn_mul
      (g := fun t => 1 - Real.cos (n * liTheta t)) (Continuous.continuousOn (by have := continuous_liTheta; fun_prop))
    exact (this.const_mul (1 / Real.pi)).congr_uIoo fun t _ => by simp only [hI]
  have hgr : ∀ a b, 0 ≤ a → 0 ≤ b → IntervalIntegrable (fun t => gfun n t * rfun t) volume a b :=
    fun a b ha hb => (intervalIntegrable_rfun ha hb).continuousOn_mul hgc.continuousOn
  have h5 : ∫ t in T..X, gfun n t * Ad t = (∫ t in T..X, hI n t) + ∫ t in T..X, gfun n t * rfun t := by
    rw [← intervalIntegral.integral_add (hiT T X hT.le hX.le) (hgr T X hT.le hX.le)]
    refine intervalIntegral.integral_congr fun t _ => ?_
    simp only [hI, gfun, rfun]
    field_simp
    ring
  -- A(T) = 1 + ∫₀^T A'
  have h6 : Afun T = 1 + (1 / (2 * Real.pi)) * (∫ t in (0 : ℝ)..T, Real.log (t / (2 * Real.pi))) +
      ∫ t in (0 : ℝ)..T, rfun t := by
    have e := integral_Ad 0 T
    rw [Afun_zero] at e
    have e2 : ∫ t in (0 : ℝ)..T, Ad t =
        (1 / (2 * Real.pi)) * (∫ t in (0 : ℝ)..T, Real.log (t / (2 * Real.pi))) + ∫ t in (0 : ℝ)..T, rfun t := by
      rw [← intervalIntegral.integral_const_mul,
        ← intervalIntegral.integral_add ((intervalIntegrable_log_div_two_pi le_rfl hT.le).const_mul _)
          (intervalIntegrable_rfun le_rfl hT.le)]
      refine intervalIntegral.integral_congr fun t _ => ?_
      simp only [rfun]; ring
    linarith
  -- ∫₀^T log = π∫₀^T hI + ∫₀^T cos·log
  have hcl : IntervalIntegrable (fun t => Real.cos (n * liTheta t) * Real.log (t / (2 * Real.pi))) volume 0 T :=
    (intervalIntegrable_log_div_two_pi le_rfl hT.le).continuousOn_mul
      (Continuous.continuousOn (by have := continuous_liTheta; fun_prop))
  have h7 : (∫ t in (0 : ℝ)..T, Real.log (t / (2 * Real.pi))) =
      Real.pi * (∫ t in (0 : ℝ)..T, hI n t) +
        ∫ t in (0 : ℝ)..T, Real.cos (n * liTheta t) * Real.log (t / (2 * Real.pi)) := by
    rw [← intervalIntegral.integral_const_mul, ← intervalIntegral.integral_add
      ((hiT 0 T le_rfl hT.le).const_mul _) hcl]
    refine intervalIntegral.integral_congr fun t _ => ?_
    simp only [hI]
    field_simp
    ring
  -- ∫₀^X hI = Iclosed n - ∫_{Ioi X} hI
  have h8 : (∫ t in (0 : ℝ)..T, hI n t) + (∫ t in T..X, hI n t) = Iclosed n - ∫ t in Ioi X, hI n t := by
    rw [intervalIntegral.integral_add_adjacent_intervals (hiT 0 T le_rfl hT.le) (hiT T X hT.le hX.le),
      intervalIntegral.integral_of_le hX.le]
    have hint := integrableOn_hI n
    have hu : Ioc (0 : ℝ) X ∪ Ioi X = Ioi 0 := Ioc_union_Ioi_eq_Ioi hX.le
    have hsplit := setIntegral_union (f := hI n) (Ioc_disjoint_Ioi_same) measurableSet_Ioi
      (hint.mono_set (by intro t ht; exact ht.1)) (hint.mono_set (Ioi_subset_Ioi hX.le))
    rw [hu] at hsplit
    have hIc := integral_I_eq n
    have e : ∫ t in Ioi (0 : ℝ), hI n t = Iclosed n := by rw [← hIc]; rfl
    linarith
  -- assemble
  have hgT : gfun n T = 2 - 2 * Real.cos (n * liTheta T) := by unfold gfun; ring
  rw [hgT] at habel
  set H0 := ∫ t in (0 : ℝ)..T, hI n t with hH0
  set CLv := ∫ t in (0 : ℝ)..T, Real.cos (n * liTheta t) * Real.log (t / (2 * Real.pi)) with hCL
  set Lg := ∫ t in (0 : ℝ)..T, Real.log (t / (2 * Real.pi)) with hLg
  have e9 : (1 / (2 * Real.pi)) * Lg * 2 = H0 + 1 / Real.pi * CLv := by
    rw [h7]; field_simp
  unfold R Rest
  rw [h1, h2, h3, habel, h4, h5]
  linarith [h8, h6, e9]

end LeanLi
