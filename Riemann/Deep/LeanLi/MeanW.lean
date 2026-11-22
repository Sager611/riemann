module

public import Riemann.Deep.LeanLi.Decomp
public import Mathlib.NumberTheory.Harmonic.Bounds

/-!
# LeanLi: `Σ_{N<n≤2N} W_T(n)² ≪ N · N(T) log T` (Dirichlet kernel + local zero counts)
-/

@[expose] public section

open Real Filter Set Finset

namespace LeanLi

/-- `2 sin(x/2) Σ_{a<n≤b} cos nx = sin((b + ½)x) - sin((a + ½)x)`. -/
lemma two_sin_mul_sum_cos (a b : ℕ) (hab : a ≤ b) (x : ℝ) :
    2 * Real.sin (x / 2) * ∑ n ∈ Ioc a b, Real.cos (n * x) =
      Real.sin ((b + 1 / 2) * x) - Real.sin ((a + 1 / 2) * x) := by
  induction b, hab using Nat.le_induction with
  | base => simp
  | succ b hab ih =>
    rw [Finset.sum_Ioc_succ_top (by omega), mul_add, ih]
    have e : 2 * Real.sin (x / 2) * Real.cos ((b + 1 : ℕ) * x) =
        Real.sin (((b + 1 : ℕ) + 1 / 2) * x) - Real.sin ((b + 1 / 2) * x) := by
      push_cast
      rw [show ((b : ℝ) + 1 + 1 / 2) * x = (b + 1) * x + x / 2 by ring,
        show ((b : ℝ) + 1 / 2) * x = (b + 1) * x - x / 2 by ring, Real.sin_add, Real.sin_sub]
      ring
    rw [e]; push_cast; ring

/-- The Dirichlet kernel bound. -/
lemma abs_sum_cos_le (N : ℕ) {x : ℝ} (hx : Real.sin (x / 2) ≠ 0) :
    |∑ n ∈ Ioc N (2 * N), Real.cos (n * x)| ≤ 1 / |Real.sin (x / 2)| := by
  have h := two_sin_mul_sum_cos N (2 * N) (by omega) x
  have hb : |2 * Real.sin (x / 2) * ∑ n ∈ Ioc N (2 * N), Real.cos (n * x)| ≤ 2 := by
    rw [h]
    have h1 := abs_sin_le_one (((2 * N : ℕ) : ℝ) + 1 / 2) |>.trans le_rfl
    calc |Real.sin ((((2 * N : ℕ) : ℝ) + 1 / 2) * x) - Real.sin (((N : ℝ) + 1 / 2) * x)| ≤
          |Real.sin ((((2 * N : ℕ) : ℝ) + 1 / 2) * x)| + |Real.sin (((N : ℝ) + 1 / 2) * x)| := abs_sub _ _
      _ ≤ 1 + 1 := add_le_add (abs_sin_le_one _) (abs_sin_le_one _)
      _ = 2 := by norm_num
  rw [abs_mul, abs_mul, abs_of_pos (by norm_num : (0 : ℝ) < 2)] at hb
  have hs : 0 < |Real.sin (x / 2)| := abs_pos.2 hx
  rw [le_div_iff₀ hs]
  nlinarith

lemma abs_sum_cos_le_N (N : ℕ) (x : ℝ) : |∑ n ∈ Ioc N (2 * N), Real.cos (n * x)| ≤ N := by
  refine (Finset.abs_sum_le_sum_abs _ _).trans ?_
  calc ∑ n ∈ Ioc N (2 * N), |Real.cos (n * x)| ≤ ∑ n ∈ Ioc N (2 * N), (1 : ℝ) :=
        Finset.sum_le_sum fun n _ => abs_cos_le_one _
    _ = N := by rw [Finset.sum_const, Nat.card_Ioc, show 2 * N - N = N by omega, nsmul_eq_mul, mul_one]

/-- `|D(x)| ≤ π/|x|` for `0 < |x| ≤ π`. -/
lemma abs_sum_cos_le_pi_div (N : ℕ) {x : ℝ} (hx0 : x ≠ 0) (hxπ : |x| ≤ π) :
    |∑ n ∈ Ioc N (2 * N), Real.cos (n * x)| ≤ π / |x| := by
  have hpos : 0 < |x| := abs_pos.2 hx0
  have hs : |x| / π ≤ |Real.sin (x / 2)| := by
    have h := sin_half_ge hpos.le hxπ
    rcases le_or_gt 0 x with h0 | h0
    · rw [abs_of_nonneg h0] at h ⊢
      rw [abs_of_nonneg (by linarith [sin_half_pos (lt_of_le_of_ne h0 (Ne.symm hx0)) (by rwa [abs_of_nonneg h0] at hxπ)])]
      exact h
    · rw [abs_of_neg h0] at h ⊢
      rw [show x / 2 = -(-x / 2) by ring, Real.sin_neg, abs_neg]
      rw [abs_of_pos (sin_half_pos (by linarith) (by rwa [abs_of_neg h0] at hxπ))]
      exact h
  have hsin : Real.sin (x / 2) ≠ 0 := by
    intro h; rw [h, abs_zero] at hs
    have : 0 < |x| / π := by positivity
    linarith
  refine (abs_sum_cos_le N hsin).trans ?_
  have : 0 < |x| / π := by positivity
  calc 1 / |Real.sin (x / 2)| ≤ 1 / (|x| / π) := one_div_le_one_div_of_le this hs
    _ = π / |x| := by field_simp

/-- `|θ(a) - θ(b)| ≥ |a - b|/(T² + 1/4)` on `[0, T]`. -/
lemma liTheta_sep {T a b : ℝ} (ha : 0 ≤ a) (haT : a ≤ T) (hb : 0 ≤ b) (hbT : b ≤ T) :
    |a - b| / (T ^ 2 + 1 / 4) ≤ |liTheta a - liTheta b| := by
  have key : ∀ x y : ℝ, 0 ≤ x → x ≤ y → y ≤ T → (y - x) / (T ^ 2 + 1 / 4) ≤ liTheta x - liTheta y := by
    intro x y hx hxy hyT
    have hD : Convex ℝ (Icc 0 T) := convex_Icc 0 T
    have hf : ∀ t, HasDerivAt (fun t => -liTheta t) (1 / (t ^ 2 + 1 / 4)) t := fun t => by
      have := (hasDerivAt_liTheta t).neg
      rw [neg_neg] at this
      exact this
    have h := hD.mul_sub_le_image_sub_of_le_deriv (f := fun t => -liTheta t)
      (fun t _ => (hf t).continuousAt.continuousWithinAt)
      (fun t _ => (hf t).differentiableAt.differentiableWithinAt)
      (C := 1 / (T ^ 2 + 1 / 4)) (fun t ht => by
        rw [interior_Icc] at ht
        rw [(hf t).deriv]
        apply one_div_le_one_div_of_le (by positivity)
        have : t ^ 2 ≤ T ^ 2 := pow_le_pow_left₀ ht.1.le ht.2.le 2
        linarith)
      x ⟨hx, le_trans hxy hyT⟩ y ⟨le_trans hx hxy, hyT⟩ hxy
    have e : (y - x) / (T ^ 2 + 1 / 4) = 1 / (T ^ 2 + 1 / 4) * (y - x) := by ring
    linarith
  rcases le_total a b with h | h
  · have := key a b ha h hbT
    rw [abs_sub_comm, abs_of_nonneg (by linarith : (0 : ℝ) ≤ b - a)]
    have hθ : liTheta b ≤ liTheta a := liTheta_anti h
    rw [abs_of_nonneg (by linarith : (0 : ℝ) ≤ liTheta a - liTheta b)]
    exact this
  · have := key b a hb h haT
    rw [abs_of_nonneg (by linarith : (0 : ℝ) ≤ a - b)]
    have hθ : liTheta a ≤ liTheta b := liTheta_anti h
    rw [abs_sub_comm, abs_of_nonneg (by linarith : (0 : ℝ) ≤ liTheta b - liTheta a)]
    exact this

section
variable (L : LitCore)
include L

/-- **Local zero count**: `#{ρ : u ≤ γ ≤ u + 1, γ ≤ T} ≤ C log T`. -/
lemma local_count : ∃ C T₀ : ℝ, 4 ≤ T₀ ∧ ∀ T, T₀ ≤ T → ∀ u : ℝ,
    ∑ ρ ∈ (zset T).filter (fun ρ => u ≤ ρ.im ∧ ρ.im ≤ u + 1), (mult ρ : ℝ) ≤ C * Real.log T := by
  obtain ⟨C₁, T₁, hT₁4, h₁⟩ := count_asymp L
  obtain ⟨C₂, T₂, hT₂4, h₂⟩ := jump_le L
  set T₃ := max T₁ T₂ with hT₃
  set K₀ := (zeroCount (T₃ + 1) : ℝ) with hK₀
  have hK0 : 0 ≤ K₀ := Nat.cast_nonneg _
  refine ⟨K₀ + 1 / Real.pi * 2 + 3 * |C₁| + |C₂|, T₃, le_trans hT₁4 (le_max_left _ _), fun T hT u => ?_⟩
  have hπ := Real.pi_gt_three
  have hπ0 := Real.pi_pos
  have hT4 : 4 ≤ T := le_trans (le_trans hT₁4 (le_max_left _ _)) hT
  have hlogT : 1 ≤ Real.log T := by
    rw [Real.le_log_iff_exp_le (by linarith)]
    have := Real.exp_one_lt_d9; linarith
  set F := (zset T).filter (fun ρ => u ≤ ρ.im ∧ ρ.im ≤ u + 1) with hF
  have hF0 : 0 ≤ ∑ ρ ∈ F, (mult ρ : ℝ) := Finset.sum_nonneg fun _ _ => Nat.cast_nonneg _
  have hRHS : K₀ ≤ (K₀ + 1 / Real.pi * 2 + 3 * |C₁| + |C₂|) * Real.log T := by
    have h0 : 0 ≤ 1 / Real.pi * 2 + 3 * |C₁| + |C₂| := by positivity
    nlinarith
  rcases lt_or_ge u T₃ with hu | hu
  · -- F ⊆ zset (T₃ + 1)
    have hsub : F ⊆ zset (T₃ + 1) := by
      intro ρ hρ
      rw [hF, Finset.mem_filter, mem_zset] at hρ
      rw [mem_zset]
      exact ⟨hρ.1.1, hρ.1.2.1, by linarith [hρ.2.2]⟩
    have := Finset.sum_le_sum_of_subset_of_nonneg (f := fun ρ => (mult ρ : ℝ)) hsub (fun _ _ _ => Nat.cast_nonneg _)
    rw [← zeroCount_eq] at this
    linarith
  rcases lt_or_ge T u with huT | huT
  · have : F = ∅ := by
      ext ρ; simp only [hF, Finset.mem_filter, mem_zset, Finset.notMem_empty, iff_false]
      rintro ⟨⟨_, _, h3⟩, h4, _⟩; linarith
    rw [this, Finset.sum_empty]
    have : 0 ≤ (K₀ + 1 / Real.pi * 2 + 3 * |C₁| + |C₂|) * Real.log T := by positivity
    linarith
  -- T₃ ≤ u ≤ T
  have hu1 : T₁ ≤ u := le_trans (le_max_left _ _) hu
  have hu2 : T₂ ≤ u := le_trans (le_max_right _ _) hu
  have hu0 : 0 < u := by linarith
  have hsub : F ⊆ zset (u + 1) \ zsetBelow u := by
    intro ρ hρ
    rw [hF, Finset.mem_filter, mem_zset] at hρ
    rw [Finset.mem_sdiff, mem_zset, mem_zsetBelow]
    refine ⟨⟨hρ.1.1, hρ.1.2.1, hρ.2.2⟩, fun h => ?_⟩
    linarith [h.2.2, hρ.2.1]
  have hbelow : zsetBelow u ⊆ zset (u + 1) := by
    intro ρ hρ; rw [mem_zsetBelow] at hρ; rw [mem_zset]; exact ⟨hρ.1, hρ.2.1, by linarith [hρ.2.2]⟩
  have hsum := Finset.sum_le_sum_of_subset_of_nonneg (f := fun ρ => (mult ρ : ℝ)) hsub (fun _ _ _ => Nat.cast_nonneg _)
  rw [Finset.sum_sdiff_eq_sub hbelow, ← zeroCount_eq] at hsum
  have hA1 := (abs_le.1 (h₁ (u + 1) (by linarith))).2
  have hA2 := (abs_le.1 (h₁ u hu1)).1
  have hJ := h₂ u hu2
  have hm := mainN_sub_le hu0 (show u ≤ u + 1 by linarith)
  have hlu : Real.log u ≤ Real.log T := Real.log_le_log hu0 huT
  have hlu1 : Real.log (u + 1) ≤ 2 * Real.log T := by
    have h1 : Real.log (u + 1) ≤ Real.log (2 * T) := Real.log_le_log (by linarith) (by linarith)
    rw [Real.log_mul (by norm_num) (by linarith)] at h1
    have : Real.log 2 ≤ Real.log T := Real.log_le_log (by norm_num) (by linarith)
    linarith
  have hlu0 : 0 ≤ Real.log u := Real.log_nonneg (by linarith)
  have hl2π : Real.log ((u + 1) / (2 * Real.pi)) ≤ Real.log (u + 1) :=
    Real.log_le_log (by positivity) (div_le_self (by linarith) (by linarith))
  have e1 : (u + 1 - u) / (2 * Real.pi) * Real.log ((u + 1) / (2 * Real.pi)) ≤ 1 / Real.pi * Real.log T := by
    have : (u + 1 - u) / (2 * Real.pi) = 1 / Real.pi / 2 := by field_simp; ring
    rw [this]
    have h0 : 0 ≤ 1 / Real.pi / 2 := by positivity
    have := mul_le_mul_of_nonneg_left (hl2π.trans hlu1) h0
    linarith
  have hC1a : C₁ * Real.log (u + 1) ≤ 2 * (|C₁| * Real.log T) := by
    have hl0 : 0 ≤ Real.log (u + 1) := Real.log_nonneg (by linarith)
    calc C₁ * Real.log (u + 1) ≤ |C₁| * Real.log (u + 1) := mul_le_mul_of_nonneg_right (le_abs_self _) hl0
      _ ≤ |C₁| * (2 * Real.log T) := mul_le_mul_of_nonneg_left hlu1 (abs_nonneg _)
      _ = 2 * (|C₁| * Real.log T) := by ring
  have hC1b : C₁ * Real.log u ≤ |C₁| * Real.log T :=
    (mul_le_mul_of_nonneg_right (le_abs_self _) hlu0).trans (mul_le_mul_of_nonneg_left hlu (abs_nonneg _))
  have hC2 : C₂ * Real.log u ≤ |C₂| * Real.log T :=
    (mul_le_mul_of_nonneg_right (le_abs_self _) hlu0).trans (mul_le_mul_of_nonneg_left hlu (abs_nonneg _))
  have hK0' : 0 ≤ K₀ * Real.log T := by positivity
  have hpi' : 0 ≤ 1 / Real.pi * Real.log T := by positivity
  have hexp : (K₀ + 1 / Real.pi * 2 + 3 * |C₁| + |C₂|) * Real.log T =
      K₀ * Real.log T + 2 * (1 / Real.pi * Real.log T) + 3 * (|C₁| * Real.log T) + |C₂| * Real.log T := by ring
  linarith

end

/-- `W² ≤ |Σ_k e^{inθ_k}|² = Σ_{j,k} cos n(θ_j - θ_k)`. -/
lemma Wv_sq_le (T : ℝ) (n : ℕ) :
    Wv T n ^ 2 ≤ ∑ j ∈ ZK T, ∑ k ∈ ZK T, Real.cos (n * (thK j - thK k)) := by
  have e : ∑ j ∈ ZK T, ∑ k ∈ ZK T, Real.cos (n * (thK j - thK k)) =
      (∑ j ∈ ZK T, Real.cos (n * thK j)) * (∑ k ∈ ZK T, Real.cos (n * thK k)) +
      (∑ j ∈ ZK T, Real.sin (n * thK j)) * (∑ k ∈ ZK T, Real.sin (n * thK k)) := by
    rw [Finset.sum_mul_sum, Finset.sum_mul_sum, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun j _ => ?_
    rw [← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun k _ => ?_
    rw [mul_sub, Real.cos_sub]
  rw [e, Wv, sq]
  nlinarith [mul_self_nonneg (∑ j ∈ ZK T, Real.sin (n * thK j))]

/-- The Dirichlet kernel `D_N(x) = Σ_{N<n≤2N} cos nx`. -/
noncomputable def Dk (N : ℕ) (x : ℝ) : ℝ := ∑ n ∈ Ioc N (2 * N), Real.cos (n * x)

lemma sum_Wv_sq_le (T : ℝ) (N : ℕ) :
    ∑ n ∈ Ioc N (2 * N), Wv T n ^ 2 ≤ ∑ j ∈ ZK T, ∑ k ∈ ZK T, |Dk N (thK j - thK k)| := by
  calc ∑ n ∈ Ioc N (2 * N), Wv T n ^ 2 ≤
        ∑ n ∈ Ioc N (2 * N), ∑ j ∈ ZK T, ∑ k ∈ ZK T, Real.cos (n * (thK j - thK k)) :=
        Finset.sum_le_sum fun n _ => Wv_sq_le T n
    _ = ∑ j ∈ ZK T, ∑ k ∈ ZK T, Dk N (thK j - thK k) := by
        rw [Finset.sum_comm]; refine Finset.sum_congr rfl fun j _ => ?_
        rw [Finset.sum_comm]; rfl
    _ ≤ _ := Finset.sum_le_sum fun j _ => Finset.sum_le_sum fun k _ => le_abs_self _

/-- Pair bound: `|D_N(θ_j - θ_k)| ≤ N·1_{|γ_j-γ_k|<1} + 2π(T²+¼)/(1 + |γ_j - γ_k|)`. -/
lemma Dk_pair_le {T : ℝ} (N : ℕ) {j k : ZeroIdx} (hj : j ∈ ZK T) (hk : k ∈ ZK T) :
    |Dk N (thK j - thK k)| ≤ (if |ordK j - ordK k| < 1 then (N : ℝ) else 0) +
      2 * (Real.pi * (T ^ 2 + 1 / 4)) / (1 + |ordK j - ordK k|) := by
  have hj0 : 0 < ordK j := j.1.2.2
  have hk0 : 0 < ordK k := k.1.2.2
  have hjT : ordK j ≤ T := mem_ZK.1 hj
  have hkT : ordK k ≤ T := mem_ZK.1 hk
  have hQ : 0 ≤ 2 * (Real.pi * (T ^ 2 + 1 / 4)) / (1 + |ordK j - ordK k|) := by positivity
  split_ifs with h
  · have := abs_sum_cos_le_N N (thK j - thK k)
    unfold Dk; linarith
  · push Not at h
    have hsep : |ordK j - ordK k| / (T ^ 2 + 1 / 4) ≤ |thK j - thK k| := liTheta_sep hj0.le hjT hk0.le hkT
    have hd : 0 < |ordK j - ordK k| := by linarith
    have hw : 0 < T ^ 2 + 1 / 4 := by positivity
    have hx : 0 < |thK j - thK k| := lt_of_lt_of_le (by positivity) hsep
    have hxπ : |thK j - thK k| ≤ Real.pi := by
      have := liTheta_pos hj0; have := liTheta_lt_pi hj0
      have := liTheta_pos hk0; have := liTheta_lt_pi hk0
      unfold thK; rw [abs_le]; constructor <;> linarith
    have h1 := abs_sum_cos_le_pi_div N (abs_pos.1 hx) hxπ
    have h2 : Real.pi / |thK j - thK k| ≤ Real.pi * (T ^ 2 + 1 / 4) / |ordK j - ordK k| := by
      rw [div_le_div_iff₀ hx hd]
      rw [div_le_iff₀ hw] at hsep
      have := mul_le_mul_of_nonneg_left hsep Real.pi_pos.le
      linarith
    have h3 : Real.pi * (T ^ 2 + 1 / 4) / |ordK j - ordK k| ≤
        2 * (Real.pi * (T ^ 2 + 1 / 4)) / (1 + |ordK j - ordK k|) := by
      rw [div_le_div_iff₀ hd (by positivity)]
      have : 0 ≤ Real.pi * (T ^ 2 + 1 / 4) := by positivity
      nlinarith
    unfold Dk at *
    linarith

/-- Fibre count: `Σ_i 1/(1+|y_i - c|) ≤ 2B H_{M+1}` when unit windows hold `≤ B` points. -/
lemma sum_inv_one_add_le {ι : Type*} (s : Finset ι) (y : ι → ℝ) (c : ℝ) (M : ℕ)
    (hM : ∀ i ∈ s, |y i - c| < M + 1) (B : ℝ)
    (hB : ∀ u : ℝ, ((s.filter (fun i => u ≤ y i ∧ y i ≤ u + 1)).card : ℝ) ≤ B) :
    ∑ i ∈ s, 1 / (1 + |y i - c|) ≤ 2 * B * ∑ m ∈ range (M + 1), 1 / ((m : ℝ) + 1) := by
  classical
  set f : ι → ℕ := fun i => ⌊|y i - c|⌋₊ with hf
  have hmaps : ∀ i ∈ s, f i ∈ range (M + 1) := by
    intro i hi
    rw [Finset.mem_range]
    have h1 := hM i hi
    have h0 : (⌊|y i - c|⌋₊ : ℝ) ≤ |y i - c| := Nat.floor_le (abs_nonneg _)
    have : (⌊|y i - c|⌋₊ : ℝ) < M + 1 := lt_of_le_of_lt h0 h1
    exact_mod_cast this
  have h1 : ∑ i ∈ s, 1 / (1 + |y i - c|) ≤ ∑ i ∈ s, 1 / ((f i : ℝ) + 1) := by
    apply Finset.sum_le_sum; intro i _
    apply one_div_le_one_div_of_le (by positivity)
    have : (⌊|y i - c|⌋₊ : ℝ) ≤ |y i - c| := Nat.floor_le (abs_nonneg _)
    simp only [hf]; linarith
  have h2 : ∑ i ∈ s, 1 / ((f i : ℝ) + 1) =
      ∑ m ∈ range (M + 1), ∑ i ∈ s.filter (fun i => f i = m), 1 / ((m : ℝ) + 1) :=
    (Finset.sum_fiberwise_of_maps_to' hmaps (fun m : ℕ => 1 / ((m : ℝ) + 1))).symm
  have h3 : ∀ m : ℕ, ((s.filter (fun i => f i = m)).card : ℝ) ≤ 2 * B := by
    intro m
    have hsub : s.filter (fun i => f i = m) ⊆ s.filter (fun i => c + m ≤ y i ∧ y i ≤ c + m + 1) ∪
        s.filter (fun i => c - m - 1 ≤ y i ∧ y i ≤ c - m - 1 + 1) := by
      intro i hi
      rw [Finset.mem_filter] at hi
      have hfl := hi.2
      simp only [hf] at hfl
      have hlo : (m : ℝ) ≤ |y i - c| := by rw [← hfl]; exact Nat.floor_le (abs_nonneg _)
      have hhi : |y i - c| < m + 1 := by rw [← hfl]; exact Nat.lt_floor_add_one _
      rw [Finset.mem_union, Finset.mem_filter, Finset.mem_filter]
      rcases le_or_gt 0 (y i - c) with h | h
      · rw [abs_of_nonneg h] at hlo hhi; left; exact ⟨hi.1, by linarith, by linarith⟩
      · rw [abs_of_neg h] at hlo hhi; right; exact ⟨hi.1, by linarith, by linarith⟩
    have hc1 := Finset.card_le_card hsub
    have hc2 := Finset.card_union_le (s.filter (fun i => c + m ≤ y i ∧ y i ≤ c + m + 1))
      (s.filter (fun i => c - m - 1 ≤ y i ∧ y i ≤ c - m - 1 + 1))
    have hb1 := hB (c + m)
    have hb2 := hB (c - m - 1)
    have : ((s.filter (fun i => f i = m)).card : ℝ) ≤
        ((s.filter (fun i => c + m ≤ y i ∧ y i ≤ c + m + 1)).card : ℝ) +
          ((s.filter (fun i => c - m - 1 ≤ y i ∧ y i ≤ c - m - 1 + 1)).card : ℝ) := by
      exact_mod_cast hc1.trans hc2
    linarith
  rw [h2] at h1
  refine h1.trans ?_
  rw [Finset.mul_sum]
  apply Finset.sum_le_sum; intro m _
  rw [Finset.sum_const, nsmul_eq_mul]
  have := h3 m
  have hm : 0 < 1 / ((m : ℝ) + 1) := by positivity
  nlinarith

/-- Points within distance `< 1` of `c` lie in two unit windows. -/
lemma card_near_le {ι : Type*} (s : Finset ι) (y : ι → ℝ) (c B : ℝ)
    (hB : ∀ u : ℝ, ((s.filter (fun i => u ≤ y i ∧ y i ≤ u + 1)).card : ℝ) ≤ B) :
    ((s.filter (fun i => |y i - c| < 1)).card : ℝ) ≤ 2 * B := by
  classical
  have hsub : s.filter (fun i => |y i - c| < 1) ⊆ s.filter (fun i => c ≤ y i ∧ y i ≤ c + 1) ∪
      s.filter (fun i => c - 1 ≤ y i ∧ y i ≤ c - 1 + 1) := by
    intro i hi
    rw [Finset.mem_filter] at hi
    rw [Finset.mem_union, Finset.mem_filter, Finset.mem_filter]
    have := abs_lt.1 hi.2
    rcases le_or_gt c (y i) with h | h
    · left; exact ⟨hi.1, h, by linarith⟩
    · right; exact ⟨hi.1, by linarith, by linarith⟩
  have hc1 := Finset.card_le_card hsub
  have hc2 := Finset.card_union_le (s.filter (fun i => c ≤ y i ∧ y i ≤ c + 1))
    (s.filter (fun i => c - 1 ≤ y i ∧ y i ≤ c - 1 + 1))
  have hb1 := hB c
  have hb2 := hB (c - 1)
  have : ((s.filter (fun i => |y i - c| < 1)).card : ℝ) ≤
      ((s.filter (fun i => c ≤ y i ∧ y i ≤ c + 1)).card : ℝ) +
        ((s.filter (fun i => c - 1 ≤ y i ∧ y i ≤ c - 1 + 1)).card : ℝ) := by
    exact_mod_cast hc1.trans hc2
  linarith

/-- Counting indexed zeros in a window. -/
lemma card_ZK_window (T u : ℝ) :
    (((ZK T).filter (fun k => u ≤ ordK k ∧ ordK k ≤ u + 1)).card : ℝ) =
      ∑ ρ ∈ (zset T).filter (fun ρ => u ≤ ρ.im ∧ ρ.im ≤ u + 1), (mult ρ : ℝ) := by
  rw [Finset.card_eq_sum_ones, Nat.cast_sum, Finset.sum_filter, Finset.sum_filter]
  have := sum_ZK T (fun x => if u ≤ x ∧ x ≤ u + 1 then (1 : ℝ) else 0)
  simp only [mul_ite, mul_one, mul_zero] at this
  push_cast
  exact this

lemma harmonic_sum_le (M : ℕ) {T : ℝ} (hM : (M : ℝ) ≤ T) :
    ∑ m ∈ range (M + 1), 1 / ((m : ℝ) + 1) ≤ 1 + Real.log (T + 1) := by
  have e : ∑ m ∈ range (M + 1), 1 / ((m : ℝ) + 1) = (harmonic (M + 1) : ℝ) := by
    simp only [harmonic]; push_cast; simp only [one_div]
  rw [e]
  have h := harmonic_le_one_add_log (M + 1)
  have h' : ((harmonic (M + 1) : ℚ) : ℝ) ≤ 1 + Real.log ((M + 1 : ℕ) : ℝ) := h
  have hl : Real.log ((M + 1 : ℕ) : ℝ) ≤ Real.log (T + 1) := by
    apply Real.log_le_log (by positivity); push_cast; linarith
  linarith

section
variable (L : LitCore)
include L

/-- **Unconditional mean square of `W_T`**: `Σ_{N<n≤2N} W_T(n)² ≤ N(T)·2B(N + 2Q(1 + log(T+1)))`. -/
theorem sumW_bound : ∃ C T₀ : ℝ, 4 ≤ T₀ ∧ 0 ≤ C ∧ ∀ T, T₀ ≤ T → ∀ N : ℕ,
    ∑ n ∈ Ioc N (2 * N), Wv T n ^ 2 ≤ (zeroCount T : ℝ) *
      (2 * (C * Real.log T) * N + 2 * (Real.pi * (T ^ 2 + 1 / 4)) * (2 * (C * Real.log T) * (1 + Real.log (T + 1)))) := by
  obtain ⟨C, T₀, hT₀, hC⟩ := local_count L
  refine ⟨|C|, T₀, hT₀, abs_nonneg _, fun T hT N => ?_⟩
  have hlog : 0 ≤ Real.log T := Real.log_nonneg (by linarith)
  set B := |C| * Real.log T with hBdef
  have hB : ∀ u : ℝ, (((ZK T).filter (fun k => u ≤ ordK k ∧ ordK k ≤ u + 1)).card : ℝ) ≤ B := by
    intro u
    rw [card_ZK_window]
    exact (hC T hT u).trans (mul_le_mul_of_nonneg_right (le_abs_self C) hlog)
  set Q := Real.pi * (T ^ 2 + 1 / 4) with hQ
  have hQ0 : 0 ≤ Q := by positivity
  have hB0 : 0 ≤ B := by positivity
  set M := ⌊T⌋₊ with hM
  have hMT : (M : ℝ) ≤ T := Nat.floor_le (by linarith)
  have hTM : T < M + 1 := Nat.lt_floor_add_one T
  have hrow : ∀ j ∈ ZK T, ∑ k ∈ ZK T, |Dk N (thK j - thK k)| ≤
      2 * B * N + 2 * Q * (2 * B * (1 + Real.log (T + 1))) := by
    intro j hj
    have hpt := fun k (hk : k ∈ ZK T) => Dk_pair_le N hj hk
    refine (Finset.sum_le_sum hpt).trans ?_
    rw [Finset.sum_add_distrib]
    have e1 : ∑ k ∈ ZK T, (if |ordK j - ordK k| < 1 then (N : ℝ) else 0) =
        N * (((ZK T).filter (fun k => |ordK k - ordK j| < 1)).card : ℝ) := by
      rw [Finset.sum_ite, Finset.sum_const_zero, add_zero, Finset.sum_const, nsmul_eq_mul, mul_comm]
      congr 3
      ext k; simp only [Finset.mem_filter, abs_sub_comm]
    have e2 : ∑ k ∈ ZK T, 2 * Q / (1 + |ordK j - ordK k|) =
        2 * Q * ∑ k ∈ ZK T, 1 / (1 + |ordK k - ordK j|) := by
      rw [Finset.mul_sum]
      refine Finset.sum_congr rfl fun k _ => ?_
      rw [abs_sub_comm]; ring
    rw [e1, e2]
    have hnear := card_near_le (ZK T) ordK (ordK j) B hB
    have hfar := sum_inv_one_add_le (ZK T) ordK (ordK j) M (fun k hk => by
      have h1 := mem_ZK.1 hk; have h2 := mem_ZK.1 hj
      have h3 : 0 < ordK k := k.1.2.2; have h4 : 0 < ordK j := j.1.2.2
      rw [abs_lt]; constructor <;> linarith) B hB
    have hH := harmonic_sum_le M hMT
    have hN0 : (0 : ℝ) ≤ N := Nat.cast_nonneg N
    have t1 : (N : ℝ) * (((ZK T).filter (fun k => |ordK k - ordK j| < 1)).card : ℝ) ≤ 2 * B * N := by
      nlinarith
    have t2 : 2 * Q * ∑ k ∈ ZK T, 1 / (1 + |ordK k - ordK j|) ≤ 2 * Q * (2 * B * (1 + Real.log (T + 1))) := by
      apply mul_le_mul_of_nonneg_left _ (by positivity)
      refine hfar.trans ?_
      exact mul_le_mul_of_nonneg_left hH (by positivity)
    linarith
  calc ∑ n ∈ Ioc N (2 * N), Wv T n ^ 2 ≤ ∑ j ∈ ZK T, ∑ k ∈ ZK T, |Dk N (thK j - thK k)| := sum_Wv_sq_le T N
    _ ≤ ∑ j ∈ ZK T, (2 * B * N + 2 * Q * (2 * B * (1 + Real.log (T + 1)))) := Finset.sum_le_sum hrow
    _ = (zeroCount T : ℝ) * (2 * B * N + 2 * Q * (2 * B * (1 + Real.log (T + 1)))) := by
        rw [Finset.sum_const, nsmul_eq_mul, card_ZK]

end

end LeanLi
