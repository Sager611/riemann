import Riemann.Deep.FirstZero.Zeta0

/-!
# Order-8 Euler–Maclaurin for ζ with an explicit remainder

For 0 < Re s, s ≠ 1, N ≥ 1:
ζ(s) = Σ_{n ≤ N} n^{-s} - N^{1-s}/(1-s) - N^{-s}/2 + s·Φ(N) + E,
|E| ≤ |s|·|(s+1)⋯(s+7)|·10⁻⁶·N^{-σ-7}/(σ+7),
Φ(N) = g₀/12 - g₂/720 + g₄/30240 - g₆/1209600, g_j = (-1)^j (s+1)⋯(s+j) N^{-s-1-j}.
Built on the first-order formula PNT.Zeta0EqZeta by seven integrations by parts on each
[n, n+1] against the normalised Bernoulli polynomials B_k(u)/k!.
-/

open Complex Set Filter Topology MeasureTheory Interval

namespace Fayed.FirstZero

/-- Normalised Bernoulli polynomials BP k u = B_k(u) / k! for 1 ≤ k ≤ 8. -/
noncomputable def BP : ℕ → ℝ → ℝ
  | 1 => fun u => u - 1 / 2
  | 2 => fun u => u ^ 2 / 2 - u / 2 + 1 / 12
  | 3 => fun u => u ^ 3 / 6 - u ^ 2 / 4 + u / 12
  | 4 => fun u => u ^ 4 / 24 - u ^ 3 / 12 + u ^ 2 / 24 - 1 / 720
  | 5 => fun u => u ^ 5 / 120 - u ^ 4 / 48 + u ^ 3 / 72 - u / 720
  | 6 => fun u => u ^ 6 / 720 - u ^ 5 / 240 + u ^ 4 / 288 - u ^ 2 / 1440 + 1 / 30240
  | 7 => fun u => u ^ 7 / 5040 - u ^ 6 / 1440 + u ^ 5 / 1440 - u ^ 3 / 4320 + u / 30240
  | 8 => fun u => u ^ 8 / 40320 - u ^ 7 / 10080 + u ^ 6 / 8640 - u ^ 4 / 17280 + u ^ 2 / 60480
      - 1 / 1209600
  | _ => fun _ => 0

lemma BP_hasDerivAt {k : ℕ} (hk1 : 1 ≤ k) (hk : k ≤ 7) (u : ℝ) :
    HasDerivAt (BP (k + 1)) (BP k u) u := by
  have hid := hasDerivAt_id' u
  interval_cases k
  · convert (((hasDerivAt_pow 2 u).div_const 2).sub (hid.div_const 2)).add_const (1 / 12 : ℝ)
      using 1
    all_goals (try funext x)
    all_goals (simp only [BP, Pi.sub_apply]; (try norm_num))
  · convert (((hasDerivAt_pow 3 u).div_const 6).sub ((hasDerivAt_pow 2 u).div_const 4)).add
      (hid.div_const 12) using 1
    all_goals (try funext x)
    all_goals (simp only [BP, Pi.sub_apply, Pi.add_apply]; (try norm_num); (try ring))
  · convert ((((hasDerivAt_pow 4 u).div_const 24).sub ((hasDerivAt_pow 3 u).div_const 12)).add
      ((hasDerivAt_pow 2 u).div_const 24)).sub_const (1 / 720 : ℝ) using 1
    all_goals (try funext x)
    all_goals (simp only [BP, Pi.sub_apply, Pi.add_apply]; (try norm_num); (try ring))
  · convert ((((hasDerivAt_pow 5 u).div_const 120).sub ((hasDerivAt_pow 4 u).div_const 48)).add
      ((hasDerivAt_pow 3 u).div_const 72)).sub (hid.div_const 720) using 1
    all_goals (try funext x)
    all_goals (simp only [BP, Pi.sub_apply, Pi.add_apply]; (try norm_num); (try ring))
  · convert (((((hasDerivAt_pow 6 u).div_const 720).sub ((hasDerivAt_pow 5 u).div_const 240)).add
      ((hasDerivAt_pow 4 u).div_const 288)).sub ((hasDerivAt_pow 2 u).div_const 1440)).add_const
      (1 / 30240 : ℝ) using 1
    all_goals (try funext x)
    all_goals (simp only [BP, Pi.sub_apply, Pi.add_apply]; (try norm_num); (try ring))
  · convert (((((hasDerivAt_pow 7 u).div_const 5040).sub ((hasDerivAt_pow 6 u).div_const 1440)).add
      ((hasDerivAt_pow 5 u).div_const 1440)).sub ((hasDerivAt_pow 3 u).div_const 4320)).add
      (hid.div_const 30240) using 1
    all_goals (try funext x)
    all_goals (simp only [BP, Pi.sub_apply, Pi.add_apply]; (try norm_num); (try ring))
  · convert ((((((hasDerivAt_pow 8 u).div_const 40320).sub
      ((hasDerivAt_pow 7 u).div_const 10080)).add ((hasDerivAt_pow 6 u).div_const 8640)).sub
      ((hasDerivAt_pow 4 u).div_const 17280)).add ((hasDerivAt_pow 2 u).div_const 60480)).sub_const
      (1 / 1209600 : ℝ) using 1
    all_goals (try funext x)
    all_goals (simp only [BP, Pi.sub_apply, Pi.add_apply]; (try norm_num); (try ring))

lemma BP_continuous {k : ℕ} (hk1 : 1 ≤ k) (hk : k ≤ 8) : Continuous (BP k) := by
  interval_cases k <;> simp only [BP] <;> fun_prop

lemma BP_one {k : ℕ} (hk2 : 2 ≤ k) (hk : k ≤ 8) : BP k 1 = BP k 0 := by
  interval_cases k <;> norm_num [BP]

/-- Sup bound for B₈(u)/8! on [0,1] (true sup 127/154828800), via Bernstein coefficients in
x = 4 (u - 1/2)². -/
lemma abs_BP8_le {u : ℝ} (h0 : 0 ≤ u) (h1 : u ≤ 1) : |BP 8 u| ≤ 1 / 1000000 := by
  set x := 4 * (u - 1 / 2) ^ 2 with hxdef
  have hx0 : 0 ≤ x := by positivity
  have hx1 : 0 ≤ 1 - x := by rw [hxdef]; nlinarith
  have e1 : 1 / 1000000 - BP 8 u = 17393 / 96768000000 * (1 - x) ^ 4
      + 4 * (4081 / 3456000000) * x * (1 - x) ^ 3 + 6 * (30019 / 18144000000) * x ^ 2 * (1 - x) ^ 2
      + 4 * (1381 / 756000000) * x ^ 3 * (1 - x) + 1381 / 756000000 * x ^ 4 := by
    simp only [BP, hxdef]; ring
  have e2 : 1 / 1000000 + BP 8 u = 176143 / 96768000000 * (1 - x) ^ 4
      + 4 * (2831 / 3456000000) * x * (1 - x) ^ 3 + 6 * (6269 / 18144000000) * x ^ 2 * (1 - x) ^ 2
      + 4 * (131 / 756000000) * x ^ 3 * (1 - x) + 131 / 756000000 * x ^ 4 := by
    simp only [BP, hxdef]; ring
  have t0 := pow_nonneg hx1 4
  have t1 := mul_nonneg hx0 (pow_nonneg hx1 3)
  have t2 := mul_nonneg (pow_nonneg hx0 2) (pow_nonneg hx1 2)
  have t3 := mul_nonneg (pow_nonneg hx0 3) hx1
  have t4 := pow_nonneg hx0 4
  rw [abs_le]
  constructor
  · have : 0 ≤ 1 / 1000000 + BP 8 u := by
      rw [e2]
      have := mul_nonneg (by norm_num : (0:ℝ) ≤ 4 * (2831 / 3456000000)) t1
      have := mul_nonneg (by norm_num : (0:ℝ) ≤ 6 * (6269 / 18144000000)) t2
      have := mul_nonneg (by norm_num : (0:ℝ) ≤ 4 * (131 / 756000000)) t3
      nlinarith
    linarith
  · have : 0 ≤ 1 / 1000000 - BP 8 u := by
      rw [e1]
      have := mul_nonneg (by norm_num : (0:ℝ) ≤ 4 * (4081 / 3456000000)) t1
      have := mul_nonneg (by norm_num : (0:ℝ) ≤ 6 * (30019 / 18144000000)) t2
      have := mul_nonneg (by norm_num : (0:ℝ) ≤ 4 * (1381 / 756000000)) t3
      nlinarith
    linarith

/-- gc s j = (-1)^j (s+1)(s+2)⋯(s+j). -/
noncomputable def gc (s : ℂ) : ℕ → ℂ
  | 0 => 1
  | j + 1 => gc s j * (-s - 1 - j)

/-- gf s j x = (d/dx)^j x^{-s-1}. -/
noncomputable def gf (s : ℂ) (j : ℕ) (x : ℝ) : ℂ := gc s j * (x : ℂ) ^ (-s - 1 - (j : ℂ))

lemma gf_hasDerivAt (s : ℂ) (j : ℕ) {x : ℝ} (hx : 0 < x) :
    HasDerivAt (gf s j) (gf s (j + 1) x) x := by
  have h := ((Complex.hasStrictDerivAt_cpow_const (c := -s - 1 - (j : ℂ))
    (ofReal_mem_slitPlane.mpr hx)).hasDerivAt).comp_ofReal
  have h2 := h.const_mul (gc s j)
  have hfun : gf s j = fun y : ℝ => gc s j * (y : ℂ) ^ (-s - 1 - (j : ℂ)) := rfl
  rw [hfun]
  convert h2 using 1
  simp only [gf, gc]
  rw [show -s - 1 - ((j + 1 : ℕ) : ℂ) = (-s - 1 - (j : ℂ)) - 1 by push_cast; ring]
  ring

lemma norm_gf {s : ℂ} {j : ℕ} {x : ℝ} (hx : 0 < x) :
    ‖gf s j x‖ = ‖gc s j‖ * x ^ (-s.re - 1 - j) := by
  simp only [gf, norm_mul, norm_cpow_eq_rpow_re_of_pos hx, sub_re, neg_re, one_re, natCast_re]

/-- One integration by parts on [n, n+1]. -/
lemma ibp_step (s : ℂ) {n : ℝ} (hn : 0 < n) {k k' : ℕ} (hk' : k' = k + 1) (hk1 : 1 ≤ k)
    (hk : k ≤ 7) {j j' : ℕ} (hj' : j' = j + 1) :
    ∫ x in n..n + 1, (BP k (x - n) : ℂ) * gf s j x =
      (BP k' 0 : ℂ) * (gf s j (n + 1) - gf s j n)
        - ∫ x in n..n + 1, (BP k' (x - n) : ℂ) * gf s j' x := by
  subst hk' hj'
  have hpos : ∀ x ∈ [[n, n + 1]], 0 < x := by
    intro x hx; rw [uIcc_of_le (by linarith)] at hx; linarith [hx.1]
  have hu : ∀ x ∈ [[n, n + 1]], HasDerivAt (gf s j) (gf s (j + 1) x) x :=
    fun x hx => gf_hasDerivAt s j (hpos x hx)
  have hv : ∀ x ∈ [[n, n + 1]],
      HasDerivAt (fun y : ℝ => (BP (k + 1) (y - n) : ℂ)) (BP k (x - n) : ℂ) x := by
    intro x _
    have := (BP_hasDerivAt hk1 hk (x - n)).comp x ((hasDerivAt_id' x).sub_const n)
    simpa using this.ofReal_comp
  have hu' : IntervalIntegrable (gf s (j + 1)) volume n (n + 1) :=
    ContinuousOn.intervalIntegrable
      (fun x hx => (gf_hasDerivAt s (j + 1) (hpos x hx)).continuousAt.continuousWithinAt)
  have hv' : IntervalIntegrable (fun x : ℝ => (BP k (x - n) : ℂ)) volume n (n + 1) :=
    (Complex.continuous_ofReal.comp ((BP_continuous hk1 (by omega)).comp
      (continuous_id.sub continuous_const))).intervalIntegrable _ _
  have key := intervalIntegral.integral_mul_deriv_eq_deriv_mul hu hv hu' hv'
  have e1 : BP (k + 1) (n + 1 - n) = BP (k + 1) 0 := by
    rw [show n + 1 - n = (1 : ℝ) by ring]; exact BP_one (by omega) (by omega)
  rw [e1, sub_self] at key
  have l1 : (∫ x in n..n + 1, (BP k (x - n) : ℂ) * gf s j x)
      = ∫ x in n..n + 1, gf s j x * (BP k (x - n) : ℂ) := by
    congr 1; ext x; ring
  have l2 : (∫ x in n..n + 1, gf s (j + 1) x * (BP (k + 1) (x - n) : ℂ))
      = ∫ x in n..n + 1, (BP (k + 1) (x - n) : ℂ) * gf s (j + 1) x := by
    congr 1; ext x; ring
  rw [l1, key, l2]
  ring

/-- The order-8 boundary functional. -/
noncomputable def EMPhi (s : ℂ) (x : ℝ) : ℂ :=
  1 / 12 * gf s 0 x - 1 / 720 * gf s 2 x + 1 / 30240 * gf s 4 x - 1 / 1209600 * gf s 6 x

lemma interval_EM8 (s : ℂ) {n : ℝ} (hn : 0 < n) :
    ∫ x in n..n + 1, (BP 1 (x - n) : ℂ) * gf s 0 x =
      EMPhi s (n + 1) - EMPhi s n - ∫ x in n..n + 1, (BP 8 (x - n) : ℂ) * gf s 7 x := by
  rw [ibp_step s hn (k := 1) (k' := 2) rfl le_rfl (by norm_num) (j := 0) (j' := 1) rfl,
    ibp_step s hn (k := 2) (k' := 3) rfl (by norm_num) (by norm_num) (j := 1) (j' := 2) rfl,
    ibp_step s hn (k := 3) (k' := 4) rfl (by norm_num) (by norm_num) (j := 2) (j' := 3) rfl,
    ibp_step s hn (k := 4) (k' := 5) rfl (by norm_num) (by norm_num) (j := 3) (j' := 4) rfl,
    ibp_step s hn (k := 5) (k' := 6) rfl (by norm_num) (by norm_num) (j := 4) (j' := 5) rfl,
    ibp_step s hn (k := 6) (k' := 7) rfl (by norm_num) (by norm_num) (j := 5) (j' := 6) rfl,
    ibp_step s hn (k := 7) (k' := 8) rfl (by norm_num) (by norm_num) (j := 6) (j' := 7) rfl]
  simp only [EMPhi, BP]
  push_cast
  ring

lemma floor_piece (s : ℂ) (k : ℕ) :
    ∫ x in (k : ℝ)..(k : ℝ) + 1, (⌊x⌋ + 1 / 2 - x) * (x : ℂ) ^ (-(s + 1)) =
      -∫ x in (k : ℝ)..(k : ℝ) + 1, (BP 1 (x - k) : ℂ) * gf s 0 x := by
  rw [← intervalIntegral.integral_neg]
  refine intervalIntegral.integral_congr_ae ?_
  have hae : ∀ᵐ x : ℝ, x ≠ (k : ℝ) + 1 := by
    rw [ae_iff]; simp [measure_singleton]
  filter_upwards [hae] with x hx hmem
  rw [uIoc_of_le (by linarith)] at hmem
  have hfl : ⌊x⌋ = (k : ℤ) := by
    rw [Int.floor_eq_iff]; push_cast
    exact ⟨hmem.1.le, lt_of_le_of_ne hmem.2 hx⟩
  have hexp : (x : ℂ) ^ (-s - 1 - ((0 : ℕ) : ℂ)) = (x : ℂ) ^ (-(s + 1)) := by
    congr 1; push_cast; ring
  simp only [gf]
  rw [hexp]
  simp only [gc, hfl, BP]
  push_cast
  ring

lemma remainder_piece (s : ℂ) (hs : 0 < s.re) {n : ℝ} (hn : 0 < n) :
    ‖∫ x in n..n + 1, (BP 8 (x - n) : ℂ) * gf s 7 x‖ ≤
      ‖gc s 7‖ / 1000000 * ((n ^ (-(s.re + 7)) - (n + 1) ^ (-(s.re + 7))) / (s.re + 7)) := by
  set a := s.re + 7 with ha
  have hpos : ∀ x ∈ [[n, n + 1]], 0 < x := by
    intro x hx; rw [uIcc_of_le (by linarith)] at hx; linarith [hx.1]
  have hb : ∀ᵐ x ∂volume, x ∈ Ioc n (n + 1) →
      ‖(BP 8 (x - n) : ℂ) * gf s 7 x‖ ≤ ‖gc s 7‖ / 1000000 * x ^ (-a - 1) := by
    refine ae_of_all _ (fun x hx => ?_)
    have hx0 : 0 < x := by linarith [hx.1]
    rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, norm_gf hx0]
    have h8 := abs_BP8_le (u := x - n) (by linarith [hx.1]) (by linarith [hx.2])
    have he : -s.re - 1 - ((7 : ℕ) : ℝ) = -a - 1 := by rw [ha]; push_cast; ring
    rw [he]
    have hxp : 0 ≤ x ^ (-a - 1) := Real.rpow_nonneg hx0.le _
    calc |BP 8 (x - n)| * (‖gc s 7‖ * x ^ (-a - 1))
        ≤ 1 / 1000000 * (‖gc s 7‖ * x ^ (-a - 1)) :=
          mul_le_mul_of_nonneg_right h8 (by positivity)
      _ = _ := by ring
  have hint : IntervalIntegrable (fun x : ℝ => ‖gc s 7‖ / 1000000 * x ^ (-a - 1)) volume n
      (n + 1) := by
    refine ContinuousOn.intervalIntegrable ?_
    exact continuousOn_const.mul
      (continuousOn_id.rpow_const fun x hx => Or.inl (hpos x hx).ne')
  refine (intervalIntegral.norm_integral_le_of_norm_le (by linarith) hb hint).trans (le_of_eq ?_)
  rw [intervalIntegral.integral_const_mul, integral_rpow]
  · have : a ≠ 0 := by rw [ha]; linarith
    rw [show -a - 1 + 1 = -a by ring]
    field_simp
    ring
  · right
    refine ⟨by intro h; linarith, ?_⟩
    intro h0
    exact (lt_irrefl _ (hpos 0 h0))

lemma gf_tendsto (s : ℂ) (hs : 0 < s.re) (j : ℕ) :
    Tendsto (fun m : ℕ => gf s j ((m : ℝ))) atTop (𝓝 0) := by
  rw [tendsto_zero_iff_norm_tendsto_zero]
  have h1 : Tendsto (fun y : ℝ => ‖gc s j‖ * y ^ (-(s.re + 1 + j))) atTop (𝓝 (‖gc s j‖ * 0)) :=
    (tendsto_rpow_neg_atTop (by positivity)).const_mul _
  rw [mul_zero] at h1
  refine (h1.comp tendsto_natCast_atTop_atTop).congr' ?_
  filter_upwards [eventually_ge_atTop 1] with m hm
  have hm' : (0 : ℝ) < m := by exact_mod_cast hm
  simp only [Function.comp, norm_gf hm']
  congr 2; ring

/-- Order-8 Euler–Maclaurin for the sawtooth integral in PNT.riemannZeta0. -/
theorem integral_floor_EM8 {s : ℂ} (hs : 0 < s.re) {N : ℕ} (hN : 1 ≤ N) :
    ‖(∫ x in Ioi (N : ℝ), (⌊x⌋ + 1 / 2 - x) * (x : ℂ) ^ (-(s + 1))) - EMPhi s N‖ ≤
      ‖gc s 7‖ / 1000000 * ((N : ℝ) ^ (-(s.re + 7)) / (s.re + 7)) := by
  set f : ℝ → ℂ := fun x => (⌊x⌋ + 1 / 2 - x) * (x : ℂ) ^ (-(s + 1)) with hf
  have hfi : IntegrableOn f (Ioi (N : ℝ)) := PNT.integrableOn_of_Zeta0_fun (by omega) hs
  set R : ℕ → ℂ := fun i => ∫ x in ((N + i : ℕ) : ℝ)..((N + i : ℕ) : ℝ) + 1,
    (BP 8 (x - (N + i : ℕ)) : ℂ) * gf s 7 x with hR
  have hNpos : (0 : ℝ) < N := by exact_mod_cast hN
  have hii : ∀ a b : ℝ, (N : ℝ) ≤ a → a ≤ b → IntervalIntegrable f volume a b := by
    intro a b ha hab
    rw [intervalIntegrable_iff_integrableOn_Ioc_of_le hab]
    exact hfi.mono_set (Ioc_subset_Ioi_self.trans (Ioi_subset_Ioi ha))
  have key : ∀ m : ℕ, ∫ x in (N : ℝ)..((N + m : ℕ) : ℝ), f x =
      EMPhi s N - EMPhi s ((N + m : ℕ) : ℝ) + ∑ i ∈ Finset.range m, R i := by
    intro m
    induction m with
    | zero => simp
    | succ m ih =>
      have hc : ((N + (m + 1) : ℕ) : ℝ) = ((N + m : ℕ) : ℝ) + 1 := by push_cast; ring
      have hle : (N : ℝ) ≤ ((N + m : ℕ) : ℝ) := by exact_mod_cast Nat.le_add_right N m
      rw [hc, ← intervalIntegral.integral_add_adjacent_intervals (hii _ _ le_rfl hle)
        (hii _ _ hle (by linarith)), ih, Finset.sum_range_succ]
      have hp : (0 : ℝ) < ((N + m : ℕ) : ℝ) := by linarith
      have h1 := floor_piece s (N + m)
      have h2 := interval_EM8 s hp
      simp only [hf] at h1 ⊢
      rw [h1, h2]
      simp only [hR]
      ring
  have hbound : ∀ m : ℕ, ‖∑ i ∈ Finset.range m, R i‖ ≤
      ‖gc s 7‖ / 1000000 * ((N : ℝ) ^ (-(s.re + 7)) / (s.re + 7)) := by
    intro m
    set h : ℕ → ℝ := fun i => ((N + i : ℕ) : ℝ) ^ (-(s.re + 7)) / (s.re + 7) with hh
    have hpiece : ∀ i, ‖R i‖ ≤ ‖gc s 7‖ / 1000000 * (h i - h (i + 1)) := by
      intro i
      have hp : (0 : ℝ) < ((N + i : ℕ) : ℝ) := by
        have : (N : ℝ) ≤ ((N + i : ℕ) : ℝ) := by exact_mod_cast Nat.le_add_right N i
        linarith
      have := remainder_piece s hs hp
      have hc : ((N + (i + 1) : ℕ) : ℝ) = ((N + i : ℕ) : ℝ) + 1 := by push_cast; ring
      simp only [hh, hc]
      rw [← sub_div]
      exact this
    calc ‖∑ i ∈ Finset.range m, R i‖ ≤ ∑ i ∈ Finset.range m, ‖R i‖ := norm_sum_le _ _
      _ ≤ ∑ i ∈ Finset.range m, ‖gc s 7‖ / 1000000 * (h i - h (i + 1)) :=
          Finset.sum_le_sum fun i _ => hpiece i
      _ = ‖gc s 7‖ / 1000000 * (h 0 - h m) := by
          rw [← Finset.mul_sum, Finset.sum_range_sub']
      _ ≤ ‖gc s 7‖ / 1000000 * h 0 := by
          apply mul_le_mul_of_nonneg_left _ (by positivity)
          have : 0 ≤ h m := by
            simp only [hh]
            apply div_nonneg (Real.rpow_nonneg (by positivity) _) (by linarith)
          linarith
      _ = _ := by simp [hh]
  have hb : Tendsto (fun m : ℕ => ((N + m : ℕ) : ℝ)) atTop atTop :=
    tendsto_natCast_atTop_atTop.comp
      (tendsto_atTop_mono (fun m => Nat.le_add_left m N) tendsto_id)
  have hlim1 : Tendsto (fun m : ℕ => ∫ x in (N : ℝ)..((N + m : ℕ) : ℝ), f x) atTop
      (𝓝 (∫ x in Ioi (N : ℝ), f x)) :=
    intervalIntegral_tendsto_integral_Ioi _ hfi hb
  have hshift : ∀ j, Tendsto (fun m : ℕ => gf s j ((N + m : ℕ) : ℝ)) atTop (𝓝 0) := by
    intro j
    exact (gf_tendsto s hs j).comp (tendsto_atTop_mono (fun m => Nat.le_add_left m N) tendsto_id)
  have hlim2 : Tendsto (fun m : ℕ => EMPhi s ((N + m : ℕ) : ℝ)) atTop (𝓝 0) := by
    have := ((((hshift 0).const_mul (1 / 12 : ℂ)).sub ((hshift 2).const_mul (1 / 720 : ℂ))).add
      ((hshift 4).const_mul (1 / 30240 : ℂ))).sub ((hshift 6).const_mul (1 / 1209600 : ℂ))
    simp only [mul_zero, sub_zero, add_zero] at this
    exact this
  have hlim : Tendsto (fun m : ℕ => ∑ i ∈ Finset.range m, R i) atTop
      (𝓝 ((∫ x in Ioi (N : ℝ), f x) - EMPhi s N + 0)) := by
    refine ((hlim1.sub_const (EMPhi s N)).add hlim2).congr' ?_
    filter_upwards with m
    rw [key m]; ring
  rw [add_zero] at hlim
  exact le_of_tendsto' hlim.norm hbound

/-- **Order-8 Euler–Maclaurin for ζ** on Re s > 0. -/
theorem zeta_EM8 {s : ℂ} (hs : 0 < s.re) (hs1 : s ≠ 1) {N : ℕ} (hN : 1 ≤ N) :
    ‖riemannZeta s - ((∑ n ∈ Finset.range (N + 1), 1 / (n : ℂ) ^ s)
        + (-(N : ℂ) ^ (1 - s)) / (1 - s) + (-(N : ℂ) ^ (-s)) / 2 + s * EMPhi s N)‖ ≤
      ‖s‖ * (‖gc s 7‖ / 1000000 * ((N : ℝ) ^ (-(s.re + 7)) / (s.re + 7))) := by
  rw [← PNT.Zeta0EqZeta (N := N) (by omega) hs hs1, PNT.riemannZeta0_apply]
  have := integral_floor_EM8 hs hN
  set I := ∫ x in Ioi (N : ℝ), (⌊x⌋ + 1 / 2 - x) * (x : ℂ) ^ (-(s + 1))
  have e : (∑ n ∈ Finset.range (N + 1), 1 / (n : ℂ) ^ s) +
      ((-(N : ℂ) ^ (1 - s)) / (1 - s) + (-(N : ℂ) ^ (-s)) / 2 + s * I) -
      ((∑ n ∈ Finset.range (N + 1), 1 / (n : ℂ) ^ s) + (-(N : ℂ) ^ (1 - s)) / (1 - s)
        + (-(N : ℂ) ^ (-s)) / 2 + s * EMPhi s N) = s * (I - EMPhi s N) := by ring
  rw [e, norm_mul]
  exact mul_le_mul_of_nonneg_left this (norm_nonneg _)

end Fayed.FirstZero
