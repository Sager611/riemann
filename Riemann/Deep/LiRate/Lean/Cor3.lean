import Riemann.Deep.LiRate.Lean.Cor

/-!
# LiRate Corollaries 2–3: the weighted Arias de Reyna criterion (proofs.md Cor 2, Cor 3)

`y_m = R(m)/m`; `Σ_m m^δ y_m²` is `Summable (fun m => m^δ (R m / m)²)`.
-/

open Finset Real Filter

namespace LiRate

/-- **Corollary 2** (from Theorem 1): `Σ_m m^{1/2} y_m² = ∞`. -/
theorem cor2 (R : ℕ → ℝ) (h : Thm1Bound R) :
    ¬ Summable (fun m : ℕ => √(m : ℝ) * (R m / m) ^ 2) := by
  intro hs
  have hnn : ∀ m : ℕ, 0 ≤ √(m : ℝ) * (R m / m) ^ 2 :=
    fun m => mul_nonneg (Real.sqrt_nonneg _) (sq_nonneg _)
  set S := ∑' m : ℕ, √(m : ℝ) * (R m / m) ^ 2 with hS
  have hlog : Tendsto (fun N : ℕ => √(Real.log N)) atTop atTop :=
    Real.tendsto_sqrt_atTop.comp (Real.tendsto_log_atTop.comp tendsto_natCast_atTop_atTop)
  obtain ⟨N, hN, hN2, hN1⟩ := ((h (c₁ / 2) (by linarith [c₁_pos])).and
    ((hlog.eventually_gt_atTop (2 * S / c₁)).and (eventually_ge_atTop 1))).exists
  have hps : ∑ m ∈ Icc 1 N, √(m : ℝ) * (R m / m) ^ 2 ≤ S :=
    hs.sum_le_tsum _ (fun m _ => hnn m)
  have hNpos : (0 : ℝ) < N := by exact_mod_cast (by omega : 0 < N)
  have hNN : 0 < (N : ℝ) * √N := by positivity
  have hlow : (∑ m ∈ Icc 1 N, R m ^ 2) / (N * √N) ≤
      ∑ m ∈ Icc 1 N, √(m : ℝ) * (R m / m) ^ 2 := by
    rw [div_eq_mul_inv, sum_mul]
    apply sum_le_sum
    intro m hm
    rw [← div_eq_mul_inv]
    rw [mem_Icc] at hm
    have hm0 : (0 : ℝ) < m := by exact_mod_cast (by omega : 0 < m)
    have hmN : (m : ℝ) ≤ N := by exact_mod_cast hm.2
    have hsm : 0 < √(m : ℝ) := Real.sqrt_pos.2 hm0
    have hm' : (m : ℝ) ≠ 0 := hm0.ne'
    have h2 : √(m : ℝ) * √m = m := Real.mul_self_sqrt hm0.le
    have e : √(m : ℝ) * (R m / m) ^ 2 = R m ^ 2 / (m * √m) := by
      rw [eq_div_iff (by positivity)]
      calc √(m : ℝ) * (R m / m) ^ 2 * (m * √m) = (√m * √m) * m * (R m / m) ^ 2 := by ring
        _ = m * m * (R m / m) ^ 2 := by rw [h2]
        _ = R m ^ 2 := by field_simp
    rw [e]
    apply div_le_div_of_nonneg_left (sq_nonneg _) (by positivity)
    exact mul_le_mul hmN (Real.sqrt_le_sqrt hmN) hsm.le hNpos.le
  have h3 : c₁ / 2 * √(Real.log N) ≤ (∑ m ∈ Icc 1 N, R m ^ 2) / (N * √N) := by
    rw [le_div_iff₀ hNN]
    nlinarith [hN]
  rw [div_lt_iff₀ c₁_pos] at hN2
  nlinarith [c₁_pos]

lemma half_le_rpow (R : ℕ → ℝ) {δ : ℝ} (hδ : 1 / 2 ≤ δ) (m : ℕ) :
    √(m : ℝ) * (R m / m) ^ 2 ≤ (m : ℝ) ^ δ * (R m / m) ^ 2 := by
  rcases Nat.eq_zero_or_pos m with h | h
  · subst h; simp
  · apply mul_le_mul_of_nonneg_right _ (sq_nonneg _)
    rw [Real.sqrt_eq_rpow]
    exact Real.rpow_le_rpow_of_exponent_le (by exact_mod_cast h) hδ

/-- Corollary 3, `⇒`: under Theorem 1, `δ ≥ 1/2` is never summable. -/
theorem not_summable_of_half_le (R : ℕ → ℝ) (h : Thm1Bound R) {δ : ℝ} (hδ : 1 / 2 ≤ δ) :
    ¬ Summable (fun m : ℕ => (m : ℝ) ^ δ * (R m / m) ^ 2) := fun hs =>
  cor2 R h (Summable.of_nonneg_of_le (fun _ => mul_nonneg (Real.sqrt_nonneg _) (sq_nonneg _))
    (half_le_rpow R hδ) hs)

/-- Dyadic decomposition of `Σ_{0<m≤2^M}`. -/
lemma sum_Ioc_dyadic (f : ℕ → ℝ) (J₀ : ℕ) (M : ℕ) (hM : J₀ ≤ M) :
    ∑ m ∈ Ioc 0 (2 ^ M), f m =
      ∑ m ∈ Ioc 0 (2 ^ J₀), f m + ∑ j ∈ Ico J₀ M, ∑ m ∈ Ioc (2 ^ j) (2 * 2 ^ j), f m := by
  induction M, hM using Nat.le_induction with
  | base => simp
  | succ M hM ih =>
    rw [← sum_Ioc_consecutive f (Nat.zero_le (2 ^ M))
      (Nat.pow_le_pow_right two_pos (by omega : M ≤ M + 1)), ih, sum_Ico_succ_top hM, pow_succ']
    ring

/-- One dyadic block of `Σ m^δ y_m²`, from Theorem 2's bound at `N = 2^j`. -/
lemma block_bound (R : ℕ → ℝ) (C : ℝ) {δ : ℝ} (hδ : δ < 1 / 2) (j : ℕ)
    (hC : ∑ n ∈ Ioc (2 ^ j) (2 * 2 ^ j), R n ^ 2 ≤ C * (((2 ^ j : ℕ) : ℝ) * √((2 ^ j : ℕ) : ℝ)) *
      (Real.log ((2 ^ j : ℕ) : ℝ) * √(Real.log ((2 ^ j : ℕ) : ℝ)))) :
    ∑ m ∈ Ioc (2 ^ j : ℕ) (2 * 2 ^ j), (m : ℝ) ^ δ * (R m / m) ^ 2 ≤
      |C| * ((j : ℝ) ^ 2 * ((2 : ℝ) ^ (δ - 1 / 2)) ^ j) := by
  have hcast : (((2 ^ j : ℕ)) : ℝ) = (2 : ℝ) ^ j := by push_cast; ring
  rw [hcast] at hC
  set X : ℝ := (2 : ℝ) ^ j with hX
  have hXpos : 0 < X := by positivity
  have hX1 : 1 ≤ X := one_le_pow₀ one_le_two
  have hterm : ∀ m ∈ Ioc (2 ^ j : ℕ) (2 * 2 ^ j),
      (m : ℝ) ^ δ * (R m / m) ^ 2 ≤ X ^ (δ - 2) * R m ^ 2 := by
    intro m hm
    rw [mem_Ioc] at hm
    have hXm : X ≤ m := by rw [hX]; exact_mod_cast hm.1.le
    have hm0 : 0 < (m : ℝ) := lt_of_lt_of_le hXpos hXm
    have e : (m : ℝ) ^ δ * (R m / m) ^ 2 = (m : ℝ) ^ (δ - 2) * R m ^ 2 := by
      rw [Real.rpow_sub hm0, Real.rpow_two, div_pow]; ring
    rw [e]
    exact mul_le_mul_of_nonneg_right
      (Real.rpow_le_rpow_of_nonpos hXpos hXm (by linarith)) (sq_nonneg _)
  have h1 := sum_le_sum hterm
  rw [← mul_sum] at h1
  have hr : X ^ (δ - 2) * (X * √X) = ((2 : ℝ) ^ (δ - 1 / 2)) ^ j := by
    have e1 : X * √X = X ^ ((1 : ℝ) + 1 / 2) := by
      rw [Real.sqrt_eq_rpow, Real.rpow_add hXpos, Real.rpow_one]
    have e2 : X ^ (δ - 2) * X ^ ((1 : ℝ) + 1 / 2) = X ^ (δ - 1 / 2) := by
      rw [← Real.rpow_add hXpos]; congr 1; ring
    have e3 : X ^ (δ - 1 / 2) = ((2 : ℝ) ^ (δ - 1 / 2)) ^ j := by
      rw [hX, ← Real.rpow_natCast (2 : ℝ) j, ← Real.rpow_mul (by norm_num),
        ← Real.rpow_natCast ((2 : ℝ) ^ (δ - 1 / 2)) j, ← Real.rpow_mul (by norm_num)]
      congr 1; ring
    rw [e1, e2, e3]
  have hlog : Real.log X * √(Real.log X) ≤ (j : ℝ) ^ 2 := by
    have e : Real.log X = j * Real.log 2 := by rw [hX, Real.log_pow]
    have hl2 : 0 < Real.log 2 := Real.log_pos one_lt_two
    have hl2' : Real.log 2 < 1 := by linarith [Real.log_two_lt_d9]
    have hj0 : (0 : ℝ) ≤ j := Nat.cast_nonneg j
    have hjj : (j : ℝ) ≤ (j : ℝ) ^ 2 := by
      rcases Nat.eq_zero_or_pos j with h | h
      · subst h; simp
      · have : (1 : ℝ) ≤ j := by exact_mod_cast h
        nlinarith
    have ha : (j : ℝ) * Real.log 2 ≤ j := by nlinarith
    have hb : √((j : ℝ) * Real.log 2) ≤ j := by
      rw [Real.sqrt_le_left hj0]; linarith
    rw [e]
    calc (j : ℝ) * Real.log 2 * √((j : ℝ) * Real.log 2) ≤ j * j :=
          mul_le_mul ha hb (Real.sqrt_nonneg _) hj0
      _ = (j : ℝ) ^ 2 := by ring
  have hlX : 0 ≤ Real.log X := Real.log_nonneg hX1
  have hrj : 0 ≤ ((2 : ℝ) ^ (δ - 1 / 2)) ^ j := by positivity
  calc ∑ m ∈ Ioc (2 ^ j : ℕ) (2 * 2 ^ j), (m : ℝ) ^ δ * (R m / m) ^ 2
      ≤ X ^ (δ - 2) * ∑ m ∈ Ioc (2 ^ j : ℕ) (2 * 2 ^ j), R m ^ 2 := h1
    _ ≤ X ^ (δ - 2) * (C * (X * √X) * (Real.log X * √(Real.log X))) :=
        mul_le_mul_of_nonneg_left hC (Real.rpow_nonneg hXpos.le _)
    _ = C * ((X ^ (δ - 2) * (X * √X)) * (Real.log X * √(Real.log X))) := by ring
    _ = C * (((2 : ℝ) ^ (δ - 1 / 2)) ^ j * (Real.log X * √(Real.log X))) := by rw [hr]
    _ ≤ |C| * (((2 : ℝ) ^ (δ - 1 / 2)) ^ j * (Real.log X * √(Real.log X))) :=
        mul_le_mul_of_nonneg_right (le_abs_self C)
          (mul_nonneg hrj (mul_nonneg hlX (Real.sqrt_nonneg _)))
    _ ≤ |C| * (((2 : ℝ) ^ (δ - 1 / 2)) ^ j * (j : ℝ) ^ 2) :=
        mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_left hlog hrj) (abs_nonneg C)
    _ = |C| * ((j : ℝ) ^ 2 * ((2 : ℝ) ^ (δ - 1 / 2)) ^ j) := by ring

/-- Summability from dyadic block bounds. -/
lemma summable_of_dyadic (f g : ℕ → ℝ) (N₀ : ℕ) (hf0 : ∀ m, 0 ≤ f m) (hg : Summable g)
    (hg0 : ∀ j, 0 ≤ g j) (hblock : ∀ j, N₀ ≤ j → ∑ m ∈ Ioc (2 ^ j) (2 * 2 ^ j), f m ≤ g j) :
    Summable f := by
  apply summable_of_sum_range_le hf0 (c := ∑ m ∈ Icc 0 (2 ^ N₀), f m + ∑' j, g j)
  intro K
  have hsub : range K ⊆ Icc 0 (2 ^ max K N₀) := by
    intro m hm
    rw [mem_range] at hm
    rw [mem_Icc]
    have h1 : K ≤ max K N₀ := le_max_left _ _
    have h2 : max K N₀ < 2 ^ max K N₀ := Nat.lt_two_pow_self
    exact ⟨Nat.zero_le _, by omega⟩
  have hdec : ∑ m ∈ Icc 0 (2 ^ max K N₀), f m = ∑ m ∈ Icc 0 (2 ^ N₀), f m +
      ∑ j ∈ Ico N₀ (max K N₀), ∑ m ∈ Ioc (2 ^ j) (2 * 2 ^ j), f m := by
    rw [Icc_eq_cons_Ioc (Nat.zero_le _), sum_cons, Icc_eq_cons_Ioc (Nat.zero_le _), sum_cons,
      sum_Ioc_dyadic f N₀ (max K N₀) (le_max_right _ _)]
    ring
  have hb : ∑ j ∈ Ico N₀ (max K N₀), ∑ m ∈ Ioc (2 ^ j) (2 * 2 ^ j), f m ≤ ∑' j, g j :=
    (sum_le_sum fun j hj => hblock j (mem_Ico.1 hj).1).trans (hg.sum_le_tsum _ fun j _ => hg0 j)
  have h1 := sum_le_sum_of_subset_of_nonneg hsub (fun m _ _ => hf0 m)
  linarith

/-- **Corollary 3**, `⇐` (from Theorem 2's upper bound, dyadic blocks): `δ < 1/2` is summable. -/
theorem cor3_summable (R : ℕ → ℝ) (h2 : Thm2Upper R) {δ : ℝ} (hδ : δ < 1 / 2) :
    Summable (fun m : ℕ => (m : ℝ) ^ δ * (R m / m) ^ 2) := by
  obtain ⟨C, N₀, hC⟩ := h2
  have hr0 : 0 ≤ (2 : ℝ) ^ (δ - 1 / 2) := Real.rpow_nonneg zero_le_two _
  have hr1 : (2 : ℝ) ^ (δ - 1 / 2) < 1 :=
    Real.rpow_lt_one_of_one_lt_of_neg one_lt_two (by linarith)
  have hr2 : ‖(2 : ℝ) ^ (δ - 1 / 2)‖ < 1 := by rw [Real.norm_eq_abs, abs_of_nonneg hr0]; exact hr1
  have hg : Summable (fun j : ℕ => |C| * ((j : ℝ) ^ 2 * ((2 : ℝ) ^ (δ - 1 / 2)) ^ j)) :=
    (summable_pow_mul_geometric_of_norm_lt_one 2 hr2).mul_left _
  have hg0 : ∀ j : ℕ, 0 ≤ |C| * ((j : ℝ) ^ 2 * ((2 : ℝ) ^ (δ - 1 / 2)) ^ j) := fun j =>
    mul_nonneg (abs_nonneg C) (mul_nonneg (sq_nonneg _) (pow_nonneg hr0 j))
  have hf0 : ∀ m : ℕ, 0 ≤ (m : ℝ) ^ δ * (R m / m) ^ 2 := fun m =>
    mul_nonneg (Real.rpow_nonneg (Nat.cast_nonneg m) δ) (sq_nonneg _)
  have hblock : ∀ j : ℕ, N₀ ≤ j →
      ∑ m ∈ Ioc (2 ^ j : ℕ) (2 * 2 ^ j), (m : ℝ) ^ δ * (R m / m) ^ 2 ≤
        |C| * ((j : ℝ) ^ 2 * ((2 : ℝ) ^ (δ - 1 / 2)) ^ j) := fun j hj =>
    block_bound R C hδ j (hC (2 ^ j) (le_trans hj (Nat.lt_two_pow_self (n := j)).le))
  exact summable_of_dyadic (fun m : ℕ => (m : ℝ) ^ δ * (R m / m) ^ 2)
    (fun j : ℕ => |C| * ((j : ℝ) ^ 2 * ((2 : ℝ) ^ (δ - 1 / 2)) ^ j)) N₀ hf0 hg hg0 hblock

/-- Lemma 0 ⇒ `Σ m^δ y_m² = ∞` for every `δ ≥ 0` (the terms do not tend to 0). -/
theorem not_summable_of_lemma0 (R : ℕ → ℝ) (h0 : Lemma0 R) {δ : ℝ} (hδ : 0 ≤ δ) :
    ¬ Summable (fun m : ℕ => (m : ℝ) ^ δ * (R m / m) ^ 2) := by
  intro hs
  obtain ⟨c, hc, hf⟩ := h0
  have hev1 : ∀ᶠ m : ℕ in atTop, (m : ℝ) ^ δ * (R m / m) ^ 2 < 1 :=
    hs.tendsto_atTop_zero.eventually (gt_mem_nhds one_pos)
  have hev2 : ∀ᶠ m : ℕ in atTop, (m : ℝ) ≤ exp (c * m) ∧ (1 : ℝ) ≤ m := by
    filter_upwards [eventually_ge_atTop (⌈2 / c ^ 2⌉₊ + 1)] with m hm
    have hm1 : (1 : ℝ) ≤ m := by exact_mod_cast (by omega : 1 ≤ m)
    have hbig : 2 / c ^ 2 ≤ m := by
      have := Nat.le_ceil (2 / c ^ 2)
      have h2 : ((⌈2 / c ^ 2⌉₊ + 1 : ℕ) : ℝ) ≤ m := by exact_mod_cast hm
      push_cast at h2; linarith
    refine ⟨?_, hm1⟩
    have hq := Real.quadratic_le_exp_of_nonneg (by positivity : 0 ≤ c * m)
    have h3 : 2 ≤ c ^ 2 * m := by
      rw [div_le_iff₀ (by positivity)] at hbig; linarith
    nlinarith [mul_nonneg (by linarith : (0 : ℝ) ≤ m) (sub_nonneg.2 h3), mul_pos hc (by linarith : (0:ℝ) < m)]
  obtain ⟨m, hm1, hm2, hm3, hm4⟩ := (hf.and_eventually (hev1.and hev2)).exists
  have hm0 : (0 : ℝ) < m := by linarith
  have hRm : (m : ℝ) ≤ |R m| := hm3.trans hm1
  have h1 : 1 ≤ (R m / m) ^ 2 := by
    rw [← sq_abs, abs_div, abs_of_pos hm0]
    have : 1 ≤ |R m| / m := (one_le_div hm0).2 hRm
    nlinarith
  have h2 : 1 ≤ (m : ℝ) ^ δ := Real.one_le_rpow hm4 hδ
  nlinarith

/-- **Corollary 3** (RH: Theorems 1 and 2): `Σ_m m^δ y_m² < ∞ ⟺ δ < 1/2`. -/
theorem cor3 (R : ℕ → ℝ) (h1 : Thm1Bound R) (h2 : Thm2Upper R) (δ : ℝ) :
    Summable (fun m : ℕ => (m : ℝ) ^ δ * (R m / m) ^ 2) ↔ δ < 1 / 2 := by
  constructor
  · intro hs
    by_contra hδ
    exact not_summable_of_half_le R h1 (not_lt.1 hδ) hs
  · exact fun hδ => cor3_summable R h2 hδ

/-- **Corollary 3**, criterion form: `RH ⟺ Σ m^δ y_m² < ∞` for some `δ ∈ [0, 1/2)`. -/
theorem rh_iff_exists (R : ℕ → ℝ) (P : Prop) (h2 : P → Thm2Upper R) (h0 : ¬P → Lemma0 R) :
    P ↔ ∃ δ ∈ Set.Ico (0 : ℝ) (1 / 2), Summable (fun m : ℕ => (m : ℝ) ^ δ * (R m / m) ^ 2) := by
  constructor
  · intro hP
    exact ⟨0, ⟨le_refl _, by norm_num⟩, cor3_summable R (h2 hP) (by norm_num)⟩
  · rintro ⟨δ, hδ, hs⟩
    by_contra hP
    exact not_summable_of_lemma0 R (h0 hP) hδ.1 hs

/-- **Corollary 3**, criterion form: `RH ⟺ Σ m^δ y_m² < ∞` for every `δ ∈ [0, 1/2)`. -/
theorem rh_iff_forall (R : ℕ → ℝ) (P : Prop) (h2 : P → Thm2Upper R) (h0 : ¬P → Lemma0 R) :
    P ↔ ∀ δ ∈ Set.Ico (0 : ℝ) (1 / 2), Summable (fun m : ℕ => (m : ℝ) ^ δ * (R m / m) ^ 2) := by
  constructor
  · exact fun hP δ hδ => cor3_summable R (h2 hP) hδ.2
  · intro h
    exact (rh_iff_exists R P h2 h0).2 ⟨0, ⟨le_refl _, by norm_num⟩, h 0 ⟨le_refl _, by norm_num⟩⟩

end LiRate
